#!/usr/bin/env python

import argparse
import subprocess
import gzip
import os
import sys
import pysam

from collections import Counter
from uuid import uuid4
from itertools import izip
from re import sub

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
    def __init__(self, gread, tread, tname, chrom, loc):
        self.gread     = gread
        self.tread     = tread  
        self.alignloc  = loc
        self.breakloc  = 0
        self.breakside = None
        self.chrom     = chrom

        self.tclass, self.tname = tname.split(':')
        
        # find breakpoint
        
        if gread.qstart < gread.rlen - gread.qend:
            self.breakloc  = gread.positions[-1] # breakpoint on right
            self.breakside = 'R'
        else:
            self.breakloc  = gread.positions[0] # breakpoint on left
            self.breakside = 'L'

        assert self.breakloc >= 0
        assert self.breakside is not None

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

    def breakpoints(self):
        return Counter([read.breakloc for read in self._splitreads])

    def subcluster_by_class(self, teclass):
        ''' return a new cluster contaning only TEs of teclass '''
        new = Cluster()
        [new.add_splitread(read) for read in self._splitreads if read.tclass == teclass]
        return new

    def __str__(self):
        locstr  = self.chrom + ':' + str(self._start) + '-' + str(self._end)
        infostr = ' '.join(map(str, ('median:', self._median, 'len:', len(self._splitreads))))
        return locstr + ' ' + infostr

    def __len__(self):
        return len(self._splitreads)


def rc(dna):
    ''' reverse complement '''
    complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
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

    return seqdict


def align(qryseq, refseq):
    ''' find best alignment with exonerate '''
    rnd = str(uuid4()) 
    tgtfa = 'tmp.' + rnd + '.tgt.fa'
    qryfa = 'tmp.' + rnd + '.qry.fa'

    tgt = open(tgtfa, 'w')
    qry = open(qryfa, 'w')

    tgt.write('>ref' + '\n' + refseq + '\n')
    qry.write('>qry' + '\n' + qryseq + '\n')

    tgt.close()
    qry.close()

    cmd = ['exonerate', '--bestn', '1', '-m', 'ungapped', '--showalignment','0', '--ryo', 'SUMMARY\t%s\t%qab\t%qae\t%tab\t%tae\t%pi\n', qryfa, tgtfa]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    best = []
    topscore = 0

    for pline in p.stdout.readlines():
        if pline.startswith('SUMMARY'):
            c = pline.strip().split()
            if int(c[1]) > topscore:
                topscore = int(c[1])
                best = c

    os.remove(tgtfa)
    os.remove(qryfa)

    return best


def bwamem(fq, ref, threads=1, width=150, sortmem=2000000000):
    ''' FIXME: add parameters to commandline '''
    rg = '@RG\tID:RCSEQ\tSM:' + fq

    sam_out  = fq + '.realign.sam'
    bam_out  = fq + '.realign.bam'
    sort_out = fq + '.realign.sorted'

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


def fetch_clipped_reads(inbamfn, minclip=50, maxaltclip=2, minmapq=10): # TODO PARAMS
    ''' FIXME: add parameters to commandline '''
    assert minclip > maxaltclip

    outfqfn = sub('.bam$', '.clipped.fq', inbamfn)
    inbam   = pysam.Samfile(inbamfn, 'rb')

    # track tid:pos:length:is_reverse to remove potential PCR dups
    used = {}

    with open(outfqfn, 'w') as outfq:
        for read in inbam.fetch():
            uid = ':'.join(map(str, (read.tid,read.pos,read.rlen,read.is_reverse)))
            unmapseq = None
            unmapqua = None

            if read.rlen - read.alen >= int(minclip): # clipped?

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

            if not read.is_unmapped and unmapseq is not None and read.mapq > int(minmapq) and uid not in used:
                infoname = ':'.join((read.qname, str(inbam.getrname(read.tid)), str(read.pos))) # help find read again
                fqstr = "\n".join(("@" + infoname, unmapseq, '+', unmapqua))
                outfq.write(fqstr + "\n")
                used[uid] = True

        inbam.close()

    return outfqfn


def build_te_splitreads(inbamfn, tebamfn, min_te_match = 0.95): # TODO PARAM
    ''' g* --> locations on genome; t* --> locations on TE '''
    inbam = pysam.Samfile(inbamfn, 'rb')
    tebam = pysam.Samfile(tebamfn, 'rb')

    ''' get sorted list of split reads '''
    splitreads = []
    for tread in tebam.fetch():
        if not tread.is_unmapped:
            tname = tebam.getrname(tread.tid)
            gname = ':'.join(tread.qname.split(':')[:-2])
            gchrom, gloc = tread.qname.split(':')[-2:]

            gloc = int(gloc)
            for gread in inbam.fetch(reference=gchrom, start=gloc, end=gloc+1):
                if (gread.qname, gread.pos, inbam.getrname(gread.tid)) == (gname, gloc, gchrom):
                        sr = SplitRead(gread, tread, tname, gchrom, gloc)
                        if sr.get_tematch() >= min_te_match:
                            splitreads.append(sr)

    inbam.close()
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


def filter_clusters(clusters, minsize=4, bothends=True): # TODO PARAM
    ''' return only clusters meeting cutoffs '''
    filtered = []

    for cluster in clusters:
        for teclass in cluster.te_classes(): # can consider peaks with more than one TE class
            subcluster = cluster.subcluster_by_class(teclass)
            reject = False
            if len(subcluster) < minsize:
                reject = True

            if not reject:
                print subcluster, subcluster.te_names(), subcluster.breakpoints()
                filtered.append(subcluster)

    return filtered


def main(args):
    # merge paired ends
    mergefq = flash_wrapper(args.pair1, args.pair2, args.maxoverlap, args.threads)

    # map merged pairs to reference
    outbam = bwamem(mergefq, args.ref, threads=args.threads)

    # load references
    terefs = read_fasta(args.telib)

    # realign clipped sequences to TE library
    clipfastq = fetch_clipped_reads(inbamfn, minclip=50, maxaltclip=2)
    clipbam   = bwamem(clipfastq, args.telib, threads=args.threads)

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


