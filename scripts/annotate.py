#!/usr/bin/env python

import sys
import pysam

if len(sys.argv) == 4:
    tbx = pysam.Tabixfile(sys.argv[2])
    
    with open(sys.argv[1]) as l1seq:
        for line in l1seq:

            if line.startswith('UUID'):
                header = line.strip().split()
                header.append(sys.argv[3])
                print '\t'.join(header)
                continue

            c = line.strip().split()
            chrom = c[1]
            start = int(c[2])
            end   = int(c[3])

            annotations = []

            if chrom in tbx.contigs:
                for rec in tbx.fetch(chrom, start, end):
                    annotations.append('|'.join(rec.strip().split()))

            annotations = list(set(annotations)) # uniqify

            if len(annotations) == 0: annotations.append('NA')

            print line.strip() + '\t' + ','.join(annotations)

else:
    sys.exit("usage: %s <l1seq.py output> <tabix> <header name>" % sys.argv[0])
