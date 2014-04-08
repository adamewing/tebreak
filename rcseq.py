#!/usr/bin/env python

import argparse
import subprocess
import gzip
import os
import sys
import pysam

from string import maketrans
from collections import Counter
from uuid import uuid4
from itertools import izip
from re import sub, search

'''
rcseq.py: Analyse RC-seq datasets.

Prereqs:
exonerate (pairwise alignments - not sure if we will need this)
bwa       (sr alignment workhorse)
pysam     (parsing SAM/BAM formatted alignments)
FLASH     (for assembling read pairs)

'''

class SplitRead:
    ''' gread = genome read, tread = TE read '''
    def __init__(self, gread, tread, tname, chrom, loc, tlen, rescue=False, breakloc=None):
        ''' if read is obtained through hard-clipped rescue, need to set breakloc manually '''
        self.gread     = gread
        self.tread     = tread  
        self.alignloc  = int(loc)
        self.breakloc  = None
        self.teside    = None # '5' or '3' (5' or 3')
        self.chrom     = chrom
        self.tlen      = int(tlen)
        self.hardclips = [] # keep hardclipped reads for later
        self.rescued   = False

        self.tclass, self.tname = tname.split(':')
        
        # find breakpoint
        if rescue:
            self.breakloc = breakloc
            self.rescued  = True
        else:
            if gread.qstart < gread.rlen - gread.qend:
                self.breakloc  = gread.positions[-1] # breakpoint on right
            else:
                self.breakloc  = gread.positions[0] # breakpoint on left

        self.dist_from_3p = self.tlen - self.tread.positions[-1]

        if self.dist_from_3p < 50:
            self.teside = '3'
        else:
            self.teside = '5'

        assert self.breakloc is not None
        assert self.teside is not None

    def get_tematch(self):
        ''' return number of mismatches / aligned length of TE sub-read '''
        nm = [value for (tag, value) in self.tread.tags if tag == 'NM'][0]
        return 1.0 - (float(nm)/float(self.tread.alen))

    def __gt__(self, other):
        ''' enables sorting of SplitRead objects '''
        if self.chrom == other.chrom:
            return self.breakloc > other.breakloc
        else:
            return self.chrom > other.chrom

    def __str__(self):
        return ','.join(map(str, ('SplitRead',self.chrom, self.breakloc, self.tname, self.tread.pos)))


class Cluster:
    ''' store and manipulate groups of SplitRead objects '''
    def __init__(self, firstread=None):
        self._splitreads = []
        self._start  = 0
        self._end    = 0
        self._median = 0
        self.chrom   = None

        if firstread is not None:
            self.add_splitread(firstread)

    def add_splitread(self, sr):
        ''' add a SplitRead and update '''
        self._splitreads.append(sr)
        if self.chrom is None:
            self.chrom = sr.chrom

        assert self.chrom == sr.chrom # clusters can't include > 1 chromosome

        if self._median > 0 and (self._median - sr.breakloc) > 1000:
            print "WARNING: Splitread", str(sr), "more than 1000 bases from median of cluster."

        ''' update statistics '''
        self._splitreads.sort()
        self._start  = self._splitreads[0].breakloc
        self._end    = self._splitreads[-1].breakloc
        self._median = self._splitreads[len(self)/2].breakloc

    def te_classes(self):
        ''' return distinct classes of TE in cluster '''
        return list(set([read.tclass for read in self._splitreads]))

    def te_names(self):
        ''' return distinct classes of TE in cluster '''
        return list(set([read.tname for read in self._splitreads]))

    def te_sides(self):
        ''' return distinct classes of TE in cluster '''
        return list(set([read.teside for read in self._splitreads]))

    def breakpoints(self):
        return Counter([read.breakloc for read in self._splitreads])

    def subcluster_by_class(self, teclass):
        ''' return a new cluster contaning only TEs of teclass '''
        new = Cluster()
        [new.add_splitread(read) for read in self._splitreads if read.tclass == teclass or read.tclass == 'POLYA']
        return new

    def get_hardclips(self):
        hclist = []
        for sr in self._splitreads:
            for hc in sr.hardclips:
                hclist.append(hc)
        return hclist

    def __str__(self):
        locstr  = self.chrom + ':' + str(self._start) + '-' + str(self._end)
        infostr = ' '.join(map(str, ('median:', self._median, 'len:', len(self._splitreads))))
        return locstr + ' ' + infostr

    def __len__(self):
        return len(self._splitreads)


def rc(dna):
    ''' reverse complement '''
    complements = maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
    return dna.translate(complements)[::-1]


def read_fasta(infa):
    ''' potato '''
    seqdict = {}

    with open(infa, 'r') as fa:
        seqid = ''
        seq   = ''
        for line in fa:
            if line.startswith('>'):
                if seq != '':
                    seqdict[seqid] = seq
                seqid = line.lstrip('>').strip()
                seq   = ''
            else:
                assert seqid != ''
                seq = seq + line.strip()

    if seqid not in seqdict and seq != '':
        seqdict[seqid] = seq

    return seqdict


def bwamem(fq, ref, threads=1, width=150, sortmem=2000000000):
    ''' FIXME: add parameters to commandline '''
    rg = '@RG\tID:RCSEQ\tSM:' + fq
    fqroot = sub('.extendedFrags.fastq.gz$', '', fq)
    fqroot = sub('.fq$', '', fqroot)

    print "DEBUG: fqroot:", fqroot

    sam_out  = '.'.join((fqroot, 'sam'))
    bam_out  = '.'.join((fqroot, 'bam'))
    sort_out = '.'.join((fqroot, 'sorted'))

    sam_cmd  = ['bwa', 'mem', '-t', str(threads), '-M', '-R', rg, ref, fq]
    bam_cmd  = ['samtools', 'view', '-bt', ref + '.bai', '-o', bam_out, sam_out]
    sort_cmd = ['samtools', 'sort', '-@', str(threads), '-m', str(sortmem), bam_out, sort_out]
    idx_cmd  = ['samtools', 'index', sort_out + '.bam']

    sys.stderr.write("running bwa-mem: " + ' '.join(sam_cmd) + "\n")

    with open(sam_out, 'w') as sam:
        p = subprocess.Popen(sam_cmd, stdout=subprocess.PIPE)
        for line in p.stdout:
            sam.write(line)

    sys.stderr.write("writing " + sam_out + " to BAM...\n")
    subprocess.call(bam_cmd)

    sys.stderr.write("sorting output: " + ' '.join(sort_cmd) + "\n")
    subprocess.call(sort_cmd)

    sys.stderr.write("indexing: " + ' '.join(idx_cmd) + "\n")
    subprocess.call(idx_cmd)

    return sort_out + '.bam'


def flash_wrapper(fq1, fq2, max_overlap, threads):
    ''' wrapper for FLASH (http://ccb.jhu.edu/software/FLASH/) '''
    out = 'RCTMP-' + str(uuid4())
    os.mkdir(out)
    args = ['flash', '-d', out, '-o', out, '-M', str(max_overlap), '-z', '-t', str(threads), fq1, fq2]
    print "calling FLASH:", args

    subprocess.call(args)

    return out + '/' + out + '.extendedFrags.fastq.gz'


def bamrec_to_fastq(read, diffseq=None, diffqual=None, diffname=None):
    ''' convert BAM record to FASTQ record: can substitute name/seq/qual using diffseq/diffqual/diffname '''
    name = read.qname
    if diffname is not None:
        name = diffname

    qual = read.qual
    if diffqual is not None:
        qual = diffqual

    seq  = read.seq
    if diffseq is not None:
        seq = diffseq

    if read.is_reverse:
        qual = qual[::-1]
        seq  = rc(seq)

    return '\n'.join(('@'+name,seq,"+",qual))


def filter_full_matches_dups(inbamfn, maxclip=10): # TODO params
    ''' return a FASTQ file of reads that didn't match the reference end-to-end '''
    outfqfn = sub('.bam$', '.imperfect.fq', inbamfn)
    inbam   = pysam.Samfile(inbamfn, 'rb')

    with open(outfqfn, 'w') as outfq:
        for read in inbam.fetch(until_eof=True):
            if not read.is_unmapped:
                if (not read.is_secondary) and read.rlen - read.alen < int(maxclip):
                    outfq.write(bamrec_to_fastq(read) + '\n')
            else:
                outfq.write(bamrec_to_fastq(read) + '\n')

    return outfqfn


def fetch_clipped_reads(inbamfn, minclip=50, maxaltclip=2): # TODO PARAMS
    ''' FIXME: add parameters to commandline '''
    assert minclip > maxaltclip

    outfqfn = sub('.bam$', '.clipped.fq', inbamfn)
    inbam   = pysam.Samfile(inbamfn, 'rb')

    # track tid:pos:length:is_reverse to remove potential PCR dups
    used  = {}

    with open(outfqfn, 'w') as outfq:
        for read in inbam.fetch():
            uid = ':'.join(map(str, (read.tid,read.pos,read.rlen,read.is_reverse)))
            unmapseq = None
            unmapqua = None

            if read.rlen - read.alen >= int(minclip): # 'soft' clipped?

                # length of 'minor' clip (want this to be small or zero - bad if just the middle part of read is aligned)
                altclip = min(read.qstart, read.rlen-read.qend)

                if altclip <= maxaltclip:
                    # get unmapped part
                    if read.qstart <= maxaltclip:
                        # (align) AAAAAAAAAAAAAAAAAA
                        # (read)  RRRRRRRRRRRRRRRRRRRRRRRRRRRR
                        unmapseq = read.seq[read.qend:]
                        unmapqua = read.qual[read.qend:]

                    else:
                        # (align)           AAAAAAAAAAAAAAAAAA
                        # (read)  RRRRRRRRRRRRRRRRRRRRRRRRRRRR
                        unmapseq = read.seq[:read.qstart]
                        unmapqua = read.qual[:read.qstart]

            if not read.is_unmapped and unmapseq is not None and uid not in used:
                infoname = ':'.join((read.qname, str(inbam.getrname(read.tid)), str(read.pos))) # help find read again
                outfq.write(bamrec_to_fastq(read, diffseq=unmapseq, diffqual=unmapqua, diffname=infoname) + '\n')
                used[uid] = True

        inbam.close()

    return outfqfn


def build_te_splitreads(refbamfn, tebamfn, teref, min_te_match = 0.9): # TODO PARAM
    ''' g* --> locations on genome; t* --> locations on TE '''
    refbam = pysam.Samfile(refbamfn, 'rb')
    tebam  = pysam.Samfile(tebamfn, 'rb')

    ''' get sorted list of split reads '''
    splitreads = []
    for tread in tebam.fetch():
        if not tread.is_unmapped:
            tname = tebam.getrname(tread.tid)
            gname = ':'.join(tread.qname.split(':')[:-2])
            gchrom, gloc = tread.qname.split(':')[-2:]

            gloc = int(gloc)
            sr = None
            hcreads = []
            for gread in refbam.fetch(reference=gchrom, start=gloc, end=gloc+1):
                if (gread.qname, gread.pos, refbam.getrname(gread.tid)) == (gname, gloc, gchrom):
                    tlen = len(teref[tname])
                    sr = SplitRead(gread, tread, tname, gchrom, gloc, tlen)
                    if sr.get_tematch() >= min_te_match:
                        splitreads.append(sr)
                if search('H', gread.cigarstring) and gread.is_secondary:
                    hcreads.append(gread)
            if sr is not None:
                [sr.hardclips.append(hc) for hc in hcreads]

    refbam.close()
    tebam.close()
    splitreads.sort()

    return splitreads


def build_te_clusters(splitreads, searchdist=100): # TODO PARAM
    ''' cluster SplitRead objects into Cluster objects and return a list of them '''
    clusters  = []

    for sr in splitreads:
        if len(clusters) == 0:
            clusters.append(Cluster(sr))

        elif clusters[-1].chrom != sr.chrom:
            clusters.append(Cluster(sr))

        else:
            #print "DEBUG: cluster median:", clusters[-1]._median, "breakloc:", sr.breakloc, "searchdist:", searchdist
            if abs(clusters[-1]._median - sr.breakloc) > searchdist:
                clusters.append(Cluster(sr))

            else:
                #print "DEBUG: adding sr", str(sr), "to cluster", str(clusters[-1])
                clusters[-1].add_splitread(sr)

    return clusters


def rescue_hardclips(clusters, refbamfn, telib, threads=1):
    ''' hard-clipped reads are often mostly TE-derived, hunt down the original read and find breakpoints '''
    hc_clusters = {}
    hc_reads    = {}
    hc_original = {}
    cn = 0

    for cluster in clusters:
        for hcread in cluster.get_hardclips():
            hc_clusters[hcread.qname] = cn
            hc_reads   [hcread.qname] = hcread
        cn += 1

    bam  = pysam.Samfile(refbamfn, 'rb')
    fqfn = sub('.bam$', '.rescue.fq', refbamfn)

    with open(fqfn, 'w') as fqout:
        for read in bam.fetch(until_eof=True):
            if not read.is_secondary:
                if read.qname in hc_clusters:
                    cliploc = read.seq.find(hc_reads[read.qname].seq)

                    if cliploc < 0:
                        cliploc = read.seq.find(rc(hc_reads[read.qname].seq))

                    assert cliploc >= 0
                    hc_original[read.qname] = read

                    fqrec = bamrec_to_fastq(read, diffseq=read.seq[:cliploc], diffqual=read.qual[:cliploc])
                    fqout.write(fqrec + '\n')

    print "Realigning hardclips to TE library..."
    bamout = bwamem(fqfn, telib, threads=threads)

    tebam = pysam.Samfile(bamout, 'rb')

    te_length = dict(zip(tebam.references, tebam.lengths))

    for teread in tebam.fetch():
        if not teread.is_secondary:
            breakloc = hc_reads[teread.qname].pos
            # break location is at right end of read
            if hc_reads[teread.qname].cigarstring.endswith('H'):
                breakloc = hc_reads[teread.qname].positions[-1]

            gread    = hc_reads[teread.qname]
            tname    = tebam.getrname(teread.tid)
            gchrom   = bam.getrname(gread.tid)
            gloc     = gread.pos
            tlen     = te_length[tname]

            sr = SplitRead(gread, teread, tname, gchrom, gloc, tlen, rescue=True, breakloc=breakloc)
            cn = hc_clusters[teread.qname]
            clusters[cn].add_splitread(sr)

    tebam.close()
    bam.close()

    return clusters


def filter_clusters(clusters, active_elts, refbamfn, minsize=4, bothends=False, mask=None): # TODO PARAM
    ''' return only clusters meeting cutoffs '''
    ''' active_elts is a list of active transposable element names (e.g. L1Hs) '''
    ''' refbamfn is the filename of the reference-aligned BAM '''
    ''' mask should be a tabix-indexed BED file of intervals to ignore (e.g. L1Hs reference locations) '''
    filtered = []
    active_elts = Counter(active_elts)

    for cluster in clusters:
        for teclass in cluster.te_classes(): # can consider peaks with more than one TE class
            subcluster = cluster.subcluster_by_class(teclass)
            reject = True
            
            # at least one alignment in a cluster has to be from an active element
            for tn in subcluster.te_names():
                if tn in active_elts:
                    reject = False

            if len(subcluster) < minsize:
                reject = True

            if bothends and len(cluster.te_sides()) < 2:
                reject = True

            if not reject:
                print subcluster, subcluster.te_names(), subcluster.breakpoints()
                filtered.append(subcluster)

    return filtered


def main(args):
    # merge paired ends
    mergefq = flash_wrapper(args.pair1, args.pair2, args.maxoverlap, args.threads)

    # filter reads unlikely to contain information about non-reference insertions (TODO: make optional)
    # optional map to reference followed by filter, not implemented yet

    # map merged pairs to reference
    refamfn = bwamem(mergefq, args.ref, threads=args.threads)

    # load references
    terefs = read_fasta(args.telib)

    # realign clipped sequences to TE library
    clipfastq = fetch_clipped_reads(refbamfn, minclip=50)
    tebamfn = bwamem(clipfastq, args.telib, threads=args.threads)
 
    # cluster split reads


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyse RC-seq data')
    parser.add_argument('-1', dest='pair1', required=True, help='fastq(.gz) containing first end reads')
    parser.add_argument('-2', dest='pair2', required=True, help='fastq(.gz) containint second end reads')
    parser.add_argument('-t', dest='threads', default=1, help='number of threads')
    parser.add_argument('-r', dest='ref', required=True, help='reference genome for bwa-mem, also expects .fai index (samtools faidx ref.fa)')
    parser.add_argument('--max-overlap', dest='maxoverlap', default=100, help='Maximum overlap used for joining paired reads with FLASH')
    parser.add_argument('--telib', dest='telib', required=True, help='TE library (BWA indexed FASTA), seq names must be CLASS:NAME')
    args = parser.parse_args()
    main(args)


