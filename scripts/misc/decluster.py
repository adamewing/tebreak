#!/usr/bin/env python

import pysam
import argparse

from collections import defaultdict as dd
from bx.intervals.intersection import Intersecter, Interval # pip install bx-python

import logging

FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def main(args):

    header = []

    support = {} 
    forest = dd(Intersecter)

    window = int(args.window)

    with open(args.table, 'r') as table:
        for i, line in enumerate(table):
            if i == 0:
                header = line.strip().split('\t')
                print line.strip()

            else:
                rec = {}

                for n, field in enumerate(line.strip().split('\t')):
                    rec[header[n]] = field

                uuid  = rec['UUID'] 
                chrom = rec['Chromosome']
                start = int(rec['Left_Extreme'])
                end   = int(rec['Right_Extreme'])

                support[uuid] = int(rec['Split_reads_5prime']) + int(rec['Split_reads_3prime'])

                forest[chrom].add_interval(Interval(start, end, value=uuid))


    with open(args.table, 'r') as table:
        for i, line in enumerate(table):
            if i == 0:
                header = line.strip().split('\t')

            else:
                rec = {}

                for n, field in enumerate(line.strip().split('\t')):
                    rec[header[n]] = field

                uuid  = rec['UUID'] 
                chrom = rec['Chromosome']
                start = int(rec['Left_Extreme'])-window
                end   = int(rec['Right_Extreme'])+window

                biggest = True

                for prox in forest[chrom].find(start, end):
                    if prox.value != uuid:
                        if support[prox.value] > int(rec['Split_reads_5prime']) + int(rec['Split_reads_3prime']):
                            biggest = False

                if biggest:
                    print line.strip()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Remove adjacent TE detections')
    parser.add_argument('-t', '--table', required=True, help='TEBreak table')
    parser.add_argument('-w', '--window', default=500)
    args = parser.parse_args()
    main(args)

