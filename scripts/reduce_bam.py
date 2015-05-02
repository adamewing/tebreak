#!/usr/bin/env python

import pysam
import sys
import argparse
import logging


def parsereads(bamfn, outfn, maxdist=10000, minclip=5):
    bam = pysam.AlignmentFile(bamfn, 'rb')
    out = pysam.AlignmentFile(outfn, 'wb', template=bam)

    for read in bam.fetch(until_eof=True):
        output = False

        if read.is_unmapped:
            output = True

        else:
            if read.mate_is_unmapped:
                output = True

            else:
                if read.rlen - read.alen >= minclip: output = True # soft-clipped

                pair_dist = abs(read.reference_start - read.next_reference_start)
                if read.tid != read.next_reference_id or pair_dist > maxdist:
                    output = True # discordant

        if read.is_duplicate: output = False

        if output: out.write(read)

    bam.close()
    out.close()


def main(args):
    assert args.bam.endswith('.bam'), "not a BAM file: %s" % args.bam

    parsereads(args.bam, args.out, maxdist=int(args.dist), minclip=int(args.minclip))


if __name__ == '__main__':
    # set up logger
    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(format=FORMAT)
 
    parser = argparse.ArgumentParser(description='extract unmapped/partially mapped/discordant reads')
    parser.add_argument('-b', '--bam', required=True, help='input BAM')
    parser.add_argument('-o', '--out', required=True, help='output BAM')
    parser.add_argument('-d', '--dist', default=10000, help='threshold distance for discordant pairs')
    parser.add_argument('-m', '--minclip', default=5, help='minimum amount of soft-clipping to output')

    args = parser.parse_args()
    main(args)