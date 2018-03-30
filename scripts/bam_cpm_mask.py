#!/usr/bin/env python

import pysam
import argparse
import multiprocessing as mp
import subprocess

import logging
logger = logging.getLogger(__name__)
FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger.setLevel(logging.INFO)

from math import log10

import pandas as pd
import scipy.stats as ss

class Segment:
    def __init__(self, chrom, start, end, cpm):
        self.chrom = chrom
        self.start = start
        self.end   = end
        self.cpm   = cpm
    
    def __lt__(self, other):
        if self.chrom != other.chrom:
            return self.chrom < other.chrom

        return self.start < other.start

    def __str__(self):
        return '%s\t%d\t%d\t%f' % (self.chrom, self.start, self.end, self.cpm)
        


def cpm(chrom, start, end, bamfn):
    bam = pysam.AlignmentFile(bamfn, 'rb')
    n = bam.mapped / float(1e6)

    count = 0

    for read in bam.fetch(chrom, start, end):
        if not read.is_secondary and read.mapq > 10: count += 1

    try:
        return (count / n) / (end-start)
    except ZeroDivisionError:
        return 0.0


def calc_seg(chrom, binstart, binend, bam):
    bin_cpm = cpm(chrom, binstart, binend, bam)

    return Segment(chrom, binstart, binend, bin_cpm)


def main(args):
    binsize = int(args.binsize)

    pool = mp.Pool(processes=int(args.procs))

    reslist = []

    with open(args.fai) as fai:
        for line in fai:
            chrom, chrlen = line.strip().split()[:2]
            chrlen = int(chrlen)

            for binstart in range(0, chrlen, binsize):
                binend = binstart + binsize
                if binend > chrlen: binend = chrlen

                res = pool.apply_async(calc_seg, [chrom, binstart, binend, args.bam])
                reslist.append(res)

    cn_segs = []
    for res in reslist:
        cn_segs.append(res.get())

    outfile = '.'.join(args.bam.split('.')[:-1]) + '.cpm.mask.txt'

    with open(outfile, 'w') as out:
        out.write('#Chrom\tStart\tEnd\tCPM\n')
        for s in sorted(cn_segs):
            out.write('%s\n' % str(s))

    data = pd.DataFrame.from_csv(outfile, sep='\t', header=0, index_col=None)

    data['z'] = ss.zscore(data['CPM'])

    highcov = data.loc[data['z'] > float(args.z)]

    highcov.to_csv(outfile, sep='\t', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Build mask from BAM based on CPM')
    parser.add_argument('--bam', required=True, help='indexed BAM')
    parser.add_argument('-f', '--fai', required=True, help='fasta index (.fai)')
    parser.add_argument('--binsize', required=True, help='bin size')
    parser.add_argument('-p', '--procs', default=1)
    parser.add_argument('-z', default=2.0, help='z-score cutoff (default = 2.0)')

    args = parser.parse_args()
    main(args)
