#!/usr/bin/env python

import pysam
import sys
import os
import argparse
import logging

logger = logging.getLogger(__name__)

def rmtag(read):
    ''' remove non-essential tags '''
    oldtags = read.tags
    newtags = []
 
    for tag in oldtags:
        if tag[0] in ('NM', 'MC', 'MD', 'MQ', 'AS', 'XS', 'RG'):
            newtags.append(tag)
 
    read.tags = newtags
 
    return read


def parsereads(bamfn, outfn, maxdist=10000, minclip=5):
    bam = pysam.AlignmentFile(bamfn, 'rb')
    out = pysam.AlignmentFile(outfn, 'wb', template=bam)

    tick = 10000000
    try:
        tick = int((bam.mapped + bam.unmapped) * 0.01)
        if tick == 0: tick = 1
        logger.debug('outputting status every %d reads (1 pct)' % tick)
    except ValueError as e:
        logger.debug('no index found, outputting status every %d reads' % tick)

    for i, read in enumerate(bam.fetch(until_eof=True)):
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

        if output: out.write(rmtag(read))

        if i % tick == 0: logger.debug('parsed %d reads, last position: %s:%d' % (i, bam.getrname(read.tid), read.pos))

    bam.close()
    out.close()


def main(args):
    if args.verbose: logger.setLevel(logging.DEBUG)

    assert args.bam.endswith('.bam'), "not a BAM file: %s" % args.bam
    if args.out is None: args.out = '.'.join(os.path.basename(args.bam).split('.')[:-1]) + '.reduced.bam'

    parsereads(args.bam, args.out, maxdist=int(args.dist), minclip=int(args.minclip))


if __name__ == '__main__':
    # set up logger
    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(format=FORMAT)
 
    parser = argparse.ArgumentParser(description='extract unmapped/partially mapped/discordant reads')
    parser.add_argument('-b', '--bam', required=True, help='input BAM')
    parser.add_argument('-o', '--out', default=None, help='output BAM (default = <input>.reduced.bam')
    parser.add_argument('-d', '--dist', default=10000, help='threshold distance for discordant pairs')
    parser.add_argument('-m', '--minclip', default=5, help='minimum amount of soft-clipping to output')

    parser.add_argument('-v', '--verbose', action='store_true')

    args = parser.parse_args()
    main(args)