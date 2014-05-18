#!/usr/bin/env python

import sys
import subprocess
import os

from uuid import uuid4
from collections import defaultdict as dd


'''
Preprocessing script to combine multiple samples from one individual into one merged sample
that can be fed to tebreak.py

input format should be individual ID followed by FASTQ paths for end1 and end2 (three columns) 
seperated by whitespace e.g.:

ID1  fastq1_1.fq fastq1_2.fq
ID1  fastq2_1.fq fastq2_2.fq
ID2  fastq3_1.fq fastq3_2.fq
ID2  fastq4_1.fq fastq4_2.fq
'''


def flash_wrapper(fq1, fq2, outdir, max_overlap, threads, uid=None):
    ''' wrapper for FLASH (http://ccb.jhu.edu/software/FLASH/) '''
    out = outdir + '/' + 'RCTMP-' + str(uuid4())
    if uid is not None:
         out = uid

    args = ['flash', '-d', os.path.dirname(out), '-o', os.path.basename(out), '-M', str(max_overlap), '-z', '-t', str(threads), fq1, fq2]
    print "INFO\t" + now() + "\tcalling FLASH:", args

    subprocess.call(args)

    return out + '.extendedFrags.fastq.gz'


def usage():
	print "usage:", sys.argv[0], "<individuals by samples list (ID[tab]FASTQ)>"


def main(infile):
	merged = dd(list)
	with open(infile, 'r') as samplemap
		for line in samplemap:
			indiv, fq1, fq2 = line.strip().split()
			

if __name__ == '__main__':
	if len(sys.argv) == 2:
		main(sys.argv[1])
	else:
		usage()