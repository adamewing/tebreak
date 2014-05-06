#!/usr/bin/env python

import sys, os
import argparse
import subprocess
import pysam

from re import sub
from uuid import uuid4


def bamtofastq(bam, samtofastq, threads=1):
    assert os.path.exists(samtofastq)
    assert bam.endswith('.bam')

    outfq = sub('bam$', 'fastq', bam)

    cmd = ['java', '-XX:ParallelGCThreads=' + str(threads), '-Xmx4g', '-jar', samtofastq + '/SamToFastq.jar', 'INPUT=' + bam, 'INTERLEAVE=true', 'FASTQ=' + outfq]
    sys.stderr.write("INFO\t" + now() + "\tconverting BAM " + bam + " to FASTQ\n")
    subprocess.call(cmd)

    return outfq


def cramtofastq(cram, cramjar, ref, threads=1):
    assert os.path.exists(cramjar)
    assert cram.endswith('.cram')

    outfq = sub('cram$', 'fastq', cram)

    cmd = ['java', '-XX:ParallelGCThreads=' + str(threads), '-Xmx4g', '-jar', cramjar, 'fastq', '-I', cram, '-R', ref, '--enumerate', '--skip-md5-check']
    sys.stderr.write("INFO\t" + now() + "\tconverting CRAM " + cram + " to FASTQ\n")

    with open(outfq, 'w') as fq:
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        for line in p.stdout:
            fq.write(line)

    return outfq


def bwamem(fq, ref, threads=1, width=150, uid=None, removefq=False):
    ''' FIXME: add parameters to commandline '''
    fqroot = sub('fastq$', '', fq)
    if uid is not None:
        fqroot = uid

    print "DEBUG: fqroot:", fqroot

    sam_out  = '.'.join((fqroot, 'sam'))
    bam_out  = '.'.join((fqroot, 'bam'))

    sam_cmd  = ['bwa', 'mem', '-t', str(threads), '-M', ref, fq]
    bam_cmd  = ['samtools', 'view', '-bt', ref + '.bai', '-o', bam_out, sam_out]

    sys.stderr.write("running bwa-mem: " + ' '.join(sam_cmd) + "\n")

    with open(sam_out, 'w') as sam:
        p = subprocess.Popen(sam_cmd, stdout=subprocess.PIPE)
        for line in p.stdout:
            sam.write(line)

    sys.stderr.write("writing " + sam_out + " to BAM...\n")
    subprocess.call(bam_cmd)
    os.remove(sam_out)

    if removefq:
    	os.remove(fq)

    return bam_out


def main(args):
    for seq in args.seqs:
        sys.stderr.write("INFO\t" + now() + "\tprocessing " + seq + "\n")
        if seq.endswith('.bam'):
            fastq = bamtofastq(seq, args.samtofastq, threads=int(args.threads))
        elif seq.endswith('.cram'):
            fastq = cramtofastq(seq, args.cramjar, threads=int(args.threads))
        elif seq.endswith('.fastq') or seq.endswith('.fastq.gz'):
            fastq = seq
        else:
            sys.stderr.write("ERROR\t" + now() + "\tunrecognized file format (extension != .bam or .cram or .fastq\n")
            sys.exit(1)


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='find relevant reads from sequence data')
	parser.add_argument(metavar='<input file bam/cram/fastq>', dest='seqs', nargs=1, help='input files')
	parser.add_argument('--samtofastq', default=None, help='path to picard SamToFastq.jar')
	parser.add_argument('--cramjar', default=None, help='path to cramtools .jar')
	parser.add_argument('-r', '--ref', required=True, help='reference fasta (needs bwa index and samtools faidx')
	parser.add_argument('-t', '--threads', default=1, help='threads for alignment')
	args = parser.parse_args()
	main(args)