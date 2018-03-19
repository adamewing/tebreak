#!/usr/bin/env python


import sys

from gzip import open
from pysam import Fastafile
from collections import defaultdict as dd
from string import maketrans


def usage():
    return "usage: %s <reference genome fasta> <refGenes.txt.gz>" % sys.argv[0]


def rc(dna):
    ''' reverse complement '''
    complements = maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
    return dna.translate(complements)[::-1]


if len(sys.argv) == 3:
    fa = Fastafile(sys.argv[1])

    assert sys.argv[2].endswith('.gz'), "refGenes.txt must be gzipped"

    genes = dd(list)

    with open(sys.argv[2], 'r') as ref:
        for line in ref:

            (bin,
            name,
            chrom,
            strand,
            txStart,
            txEnd,
            cdsStart,
            cdsEnd,
            exonCount,
            exonStarts,
            exonEnds,
            score,
            name2,
            cdsStartStat,
            cdsEndStat,
            exonFrames) = line.strip().split()

            exonStarts = map(int, exonStarts.split(',')[:-1])
            exonEnds   = map(int, exonEnds.split(',')[:-1])

            assert len(exonEnds) == len(exonStarts) == int(exonCount)

            seq = ''

            for start, end in zip(exonStarts, exonEnds):
                if chrom in fa.references:
                    seq += fa.fetch(chrom, start, end)

            if strand == '-':
                seq = rc(seq)

            #seq = seq + 'A'*100 # polyadenylate

            if seq:
                genes[name2].append(seq)


    for name in genes:
        for i, tx in enumerate(genes[name]):
            print ">%s.%d\n%s" % (name, i, tx)

else:
    sys.exit(usage())
