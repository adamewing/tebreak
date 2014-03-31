#!/usr/bin/env python

import argparse
import subprocess
import gzip
import os
import sys
import pysam

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
        self.gread = gread
        self.tread = tread

        self.tname = tname
        self.chrom = chrom
        self.alignloc   = loc
        self.breakloc   = 0

        # find breakpoint
        
        if gread.qstart < gread.rlen - gread.qend:
            # (align) AAAAAAAAAAAAAAAAAA
            # (read)  RRRRRRRRRRRRRRRRRRRRRRRRRRRR
            self.breakloc = gread.positions[-1]

        else:
            # (align)           AAAAAAAAAAAAAAAAAA
            # (read)  RRRRRRRRRRRRRRRRRRRRRRRRRRRR
            self.breakloc = gread.positions[0]

        assert self.breakloc >= 0

    def __gt__(self, other):
        ''' enables sorting of SplitRead objects '''
        if self.chrom == other.chrom:
            return self.breakloc > other.breakloc
        else:
            return self.chrom > other.chrom



def rc(dna):
    complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
    return dna.translate(complements)[::-1]


def read_fasta(infa):
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


def bestalign(qryseq, reflist):
    topscore = 0
    topref   = [] 
    for refid, refseq in reflist.iteritems():
        aln = align(qryseq, refseq)
        score = 0
        if aln:
            score = aln[1]
        if score > topscore:
            topscore = score
            topref   = aln
            topref.append(refid)

    return topref #FIXME ... testing currently


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


def fetch_clipped_reads(inbamfn, minclip=50, maxaltclip=2, minmapq=10):
    ''' FIXME: add parameters to commandline '''
    assert minclip > maxaltclip

    outfqfn = sub('.bam$', '.clipped.fq', inbamfn)
    inbam   = pysam.Samfile(inbamfn, 'rb')

    with open(outfqfn, 'w') as outfq:
        for read in inbam.fetch():
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

            if not read.is_unmapped and unmapseq is not None and read.mapq > int(minmapq):
                infoname = ':'.join((read.qname, str(inbam.getrname(read.tid)), str(read.pos))) # help find read again
                fqstr = "\n".join(("@" + infoname, unmapseq, '+', unmapqua))
                outfq.write(fqstr + "\n")

        inbam.close()

    return outfqfn


def build_te_clusters(inbamfn, tebamfn):
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
                    splitreads.append(SplitRead(gread, tread, tname, gchrom, gloc))

    inbam.close()
    tebam.close()
    splitreads.sort()

    return splitreads # FIXME 


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
    parser.add_argument('--telib', dest='telib', required=True, help='TE library (BWA indexed FASTA)')
    args = parser.parse_args()
    main(args)


