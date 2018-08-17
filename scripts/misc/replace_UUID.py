#!/usr/bin/env python

import pysam
import argparse

from collections import defaultdict as dd
from bx.intervals.intersection import Intersecter, Interval

def main(args):

    header = []

    forest = dd(Intersecter)


    with open(args.bedlike, 'r') as bed:
        for line in bed:
            chrom, start, end, uuid = line.strip().split()[:4]
            start = int(start)
            end   = int(end)

            assert start <= end, 'improperly formed BED-like file, start > end: %s' % line.strip()

            forest[chrom].add_interval(Interval(start, end, value=uuid))


    with open(args.table, 'r') as table:
        for i, line in enumerate(table):
            if i == 0:
                print line.strip()
                header = line.strip().split('\t')

            else:
                rec = {}

                for n, field in enumerate(line.strip().split('\t')):
                    rec[header[n]] = field

                chrom = rec['Chromosome']
                start = int(rec['Left_Extreme'])
                end   = int(rec['Right_Extreme'])


                uuid = rec['UUID']

                if chrom in forest:
                    for rec in forest[chrom].find(start, end):
                        uuid = rec.value
                        print '%s\t%s' % (uuid, '\t'.join(line.strip().split()[1:]))
                        continue


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Change UUIDs in tebreak table')
    parser.add_argument('-t', '--table', required=True, help='TEBreak table')
    parser.add_argument('-b', '--bedlike', required=True, help='BED-like file (1st three columns are chrom, start, end)')
    args = parser.parse_args()
    main(args)

