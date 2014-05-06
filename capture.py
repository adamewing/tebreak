#!/usr/bin/env python

import sys, os
import argparse
import subprocess
import pysam

from re import sub
from uuid import uuid4


def bamtofastq(bam, picard, threads=1):
	assert os.path.exists(picard + '/SamToFastq.jar')
	assert bam.endswith('.bam')

	outfq = sub('bam$', 'fastq', bam)

	cmd = ['java', '-XX:ParallelGCThreads=' + str(threads), '-Xmx4g', '-jar', picard + '/SamToFastq.jar', 'INPUT=' + args.bam, 'INTERLEAVE=true', 'FASTQ=' + outfq]

    sys.stderr.write("DEBUG\t" + now() + "\tsplit cmd: " + ' '.join(cmd) + "\n")

    sys.stderr.write("INFO\t" + now() + "\tsplitting BAM " + args.bam + " into FASTQ files by readgroup\n")
    subprocess.call(cmd)

    return outfq


def crambam(cramjar, cramfn, chrom, start, end, ref, threads=1, outbam=str(uuid4())+'.bam'):
    ''' make local BAM '''
    assert os.path.exists(cramjar)

    region = str(chrom) + ':' + str(start) + '-' + str(end)

    args = ['java', '-XX:ParallelGCThreads=' + str(threads), '-Xmx4g', '-jar', cramjar, 'bam', '-I', cramfn, '-O', outbam, '-R', ref, '--skip-md5-check', region]
    subprocess.call(args)

    return outbam


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
	for seq in seqs:
		sys.stderr.write("INFO\t" + now() + "\tprocessing " + seq "\n")


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='find relevant reads from sequence data')
	parser.add_argument(metavar='<input file(s) (bam or fastq)>', dest='seqs', nargs='+', help='input files')
	parser.add_argument('--picard', required=True, help='path to picard JARs')
	parser.add_argument('--cramjar', required=True, help='path to cramtools .jar')
	parser.add_argument('-r', '--ref', required=True, help='reference fasta (needs bwa index and samtools faidx')
	parser.add_argument('-t', '--threads', default=1, help='threads for alignment')
	args = parser.parse_args()
	main(args)