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


def rc(dna):
    complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
    return dna.translate(complements)[::-1]


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

    cmd = ['exonerate', '--bestn', '1', '-m', 'ungapped', '--showalignment','0', '--ryo', 'SUMMARY\t%s\t%qab\t%qae\t%tab\t%tae\n', qryfa, tgtfa]
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
    pass

def loadrefs(fa):
    pass

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


def fetch_clipped_reads(inbamfn, minclip=50, maxaltclip=2, refs=None):
    ''' FIXME: ass parameters to commandline '''
    assert minclip > maxaltclip

    outbamfn = sub('.bam$', '.clipped.bam', inbamfn)

    inbam  = pysam.Samfile(inbamfn, 'rb')
    outbam = pysam.Samfile(outbamfn, 'wb', template=inbam)

    for read in inbam.fetch():
        if read.rlen - read.alen >= int(minclip): # clipped?

            # length of 'minor' clip (want this to be small or zero - bad if just the midle part of read is aligned)
            altclip = min(read.qstart, read.rlen-read.qend)

            if altclip <= maxaltclip:
                # get unmapped part

                if read.qstart <= maxaltclip:
                    # (align) AAAAAAAAAAAAAAAAAA
                    # (read)  RRRRRRRRRRRRRRRRRRRRRRRRRRRR
                    unmapseq = read.seq[read.qend:]

                else:
                    # (align)           AAAAAAAAAAAAAAAAAA
                    # (read)  RRRRRRRRRRRRRRRRRRRRRRRRRRRR
                    unmapseq = read.seq[:read.qstart]

            if refs is not None:
                # TODO

                print "read.rlen, read.alen, read.qstart, read.qend, altclip, read.seq, unmapseq" 
                print read.rlen, read.alen, read.qstart, read.qend, altclip, read.seq, unmapseq 
                outbam.write(read)

    inbam.close()
    outbam.close()

    return outbamfn


def main(args):
    # merge paired ends
    mergefq = flash_wrapper(args.pair1, args.pair2, args.maxoverlap, args.threads)

    # map merged pairs to reference
    outbam  = bwamem(mergefq, args.ref, threads=args.threads)

    # load references
    # TODO

    # postprocess alignemnts
    # TODO

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyse RC-seq data')
    parser.add_argument('-1', dest='pair1', required=True, help='fastq(.gz) containing first end reads')
    parser.add_argument('-2', dest='pair2', required=True, help='fastq(.gz) containint second end reads')
    parser.add_argument('-t', dest='threads', default=1, help='number of threads')
    parser.add_argument('-r', dest='ref', required=True, help='reference genome for bwa-mem, also expects .fai index (samtools faidx ref.fa)')
    parser.add_argument('--max-overlap', dest='maxoverlap', default=100, help='Maximum overlap used for joining paired reads with FLASH')
    parser.add_argument('--telib', dest='telib', required=True, help='TE library (FASTA)')
    args = parser.parse_args()
    main(args)


