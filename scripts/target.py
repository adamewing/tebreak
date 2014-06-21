#!/usr/bin/env python

import os
import sys
import pysam
import datetime

from string import maketrans
from uuid import uuid4
from re import sub


## pulled in from tebreak.py ##

def now():
    return str(datetime.datetime.now())


def rc(dna):
    ''' reverse complement '''
    complements = maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
    return dna.translate(complements)[::-1]

###############################


def maketargetfastq(bamfn, bedfn):
    assert bamfn.endswith('.bam'), "not bam file: %r" % bamfn
    bam  = pysam.Samfile(bamfn, 'rb')
    fqfn = sub('bam$', str(uuid4()).split('-')[0] + '.fastq', os.path.basename(bamfn))
   
    print "INFO\t" + now() + "\tgetting readnames for target regions"
    rnames = {}

    with open(bedfn, 'r') as bed:
        for line in bed:
            chrom, start, end = line.strip().split()[:3]
            start, end = int(start), int(end)
            if chrom in bam.references:
                for read in bam.fetch(chrom, start, end):
                    rnames[read.qname] = True

    bam.close()
    bam = pysam.Samfile(bamfn, 'rb')

    print "INFO\t" + now() + "\tselected " + str(len(rnames)) + " read names based on coordinates"
    print "INFO\t" + now() + "\tdumping selected reads to fastq: " + fqfn

    with open(fqfn, 'w') as fq:
        for read in bam.fetch(until_eof=True):
            if not read.is_secondary:
                if read.qname in rnames:
                    fq.write('@' + read.qname + ':' + str(uuid4()).split('-')[0] + '\n')
                    if read.is_reverse:
                        fq.write(rc(read.seq) + '\n+\n' + read.qual[::-1] + '\n')
                    else:
                        fq.write(read.seq + '\n+\n' + read.qual + '\n')
    bam.close()
    return fqfn

if len(sys.argv) == 3:
    fastq = maketargetfastq(*sys.argv[1:])
else:
    print "usage:",sys.argv[0],"<bam> <bed> (script for testing)"
