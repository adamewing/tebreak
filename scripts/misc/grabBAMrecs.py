#!/usr/bin/env python

import sys
import pysam
import os
import re

from collections import defaultdict as dd

import logging
logger = logging.getLogger(__name__)
FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger.setLevel(logging.INFO)


def find_mate(read, bam):
    ''' AlignmentFile.mate() can return a non-primary alignment, so use this function instead '''
    chrom = read.next_reference_name
    for rec in bam.fetch(chrom, read.next_reference_start, read.next_reference_start+1):
        if rec.query_name == read.query_name and rec.reference_start == read.next_reference_start:
            if not rec.is_secondary and bin(rec.flag & 2048) != bin(2048):
                if rec.is_read1 != read.is_read1:
                    return rec
    return None


if len(sys.argv) == 3:
    inbam = pysam.AlignmentFile(sys.argv[1], 'rb')

    outfn = '.'.join(os.path.basename(sys.argv[1]).split('.')[:-1]) + '.' + re.sub(':', '_', sys.argv[2]) + '.bam'

    outbam = pysam.AlignmentFile(outfn, 'wb', template=inbam)

    seen = dd(list)

    for read in inbam.fetch(region=sys.argv[2]):
        if not read.is_supplementary and not read.is_secondary and not read.mate_is_unmapped:
            outbam.write(read)
            seen[read.qname].append(read.is_read1)

    seen_pairs = 0
    seen_alone = 0

    for qname, pair in seen.iteritems():
        assert len(set(pair)) <= 2
        
        if len(set(pair)) == 2:
            seen_pairs += 1
        if len(set(pair)) == 1:
            seen_alone += 1

    logger.info('%d pairs inside and %d mates outside region %s' % (seen_pairs, seen_alone, sys.argv[2]))

    matebam = pysam.AlignmentFile(sys.argv[1], 'rb') 

    for read in inbam.fetch(region=sys.argv[2]):
        if not read.is_supplementary and not read.is_secondary and not read.mate_is_unmapped:
            assert read.qname in seen
            if len(set(seen[read.qname])) == 1:
                mate = find_mate(read, matebam)
                if mate is not None:
                    outbam.write(mate)


else:
    sys.exit('usage: %s <BAM> <region chrom:start-end>' % sys.argv[0])
