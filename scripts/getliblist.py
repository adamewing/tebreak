#!/usr/bin/env python

import sys
import os

def getlibs(invcf):
    rgs = {}
    with open(invcf, 'r') as vcf:
        for line in vcf:
            if not line.startswith('#'):
                chrom, pos, id, ref, alt, qual, filter, info, format, sample = line.strip().split('\t')
                for rg in sample.split(':')[-1].split(','):
                    rgs[rg] = True

    return rgs.keys()


if len(sys.argv) == 2:
    rgs = []
    with open(sys.argv[1], 'r') as vcflist:
        for vcf in vcflist:
            vcf = vcf.strip()
            assert os.path.exists(vcf), "VCF not found: " + vcf
            for rg in getlibs(vcf):
                rgs.append(rg)

    print '\n'.join(sorted(list(set(rgs))))
            
else:
    print "usage:", sys.argv[0], "<tebreak output vcf list in a file>"
