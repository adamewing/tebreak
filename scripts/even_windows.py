#!/usr/bin/env python

import pysam
import argparse
import logging


FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

def main(args):
    bam = pysam.AlignmentFile(args.bam, 'rb')

    read_count = 0
    last_chrom = None
    first_pos = 0
    prev_rec = None

    reads_per_window = bam.mapped / int(args.windows)

    assert reads_per_window > 0, "too few reads or too many windows"

    for rec in bam.fetch():
        if not rec.is_unmapped:
            read_count += 1

            if last_chrom is None:
                last_chrom = rec.reference_name

            if read_count >= reads_per_window or rec.reference_name != last_chrom:
                print '%s\t%d\t%d\t%d' % (last_chrom, first_pos, prev_rec.reference_end, read_count)

                read_count = 0
                last_chrom = rec.reference_name
                first_pos = rec.reference_start

                if int(args.padding) > 0:
                    if first_pos - int(args.padding) < 0:
                        first_pos = 0
                    else:
                        first_pos -= int(args.padding)

            prev_rec = rec

    print '%s\t%d\t%d\t%d' % (last_chrom, first_pos, prev_rec.reference_end, read_count)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='assign windows to BAMs such that reads are evenly distributed among the windows')
    parser.add_argument('-b', '--bam', required=True, help='coordinate-sorted BAM')
    parser.add_argument('-w', '--windows', required=True)
    parser.add_argument('-p', '--padding', default=0)
    args = parser.parse_args()
    main(args)