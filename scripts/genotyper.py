#!/usr/bin/env python


import sys
import argparse
import pysam


def break_count(bamfn, chrom, poslist, minpad=5, flex=1):
    bam = pysam.AlignmentFile(bamfn, 'rb')

    altcount = 0
    refcount = 0
    discards = 0

    tsd_start = min(poslist)
    tsd_end   = max(poslist)

    tsd_len = tsd_end - tsd_start

    for read in bam.fetch(chrom, tsd_start-minpad, tsd_end+minpad):
        if read.is_unmapped or read.is_duplicate or read.is_secondary:
            continue

        if read.mapq < 10:  # param
            continue

        rclip = len(read.seq) - read.query_alignment_end 
        lclip = read.query_alignment_start

        rbreak = 0
        lbreak = 0

        if rclip > max(tsd_len, minpad):
            rbreak = read.reference_end

        if lclip > max(tsd_len, minpad):
            lbreak = read.reference_start

        support_alt = False

        for pos in poslist: # does this read support a breakpoint in the list?
            if (rbreak >= pos-flex and rbreak <= pos+flex) or (lbreak >= pos-flex and lbreak <= pos+flex):
                support_alt = True

        if support_alt:
            altcount += 1

        else:
            #for pos in poslist: # does this read span a breakpoint in the list?
            if read.alen == len(read.seq):
                if read.reference_start < tsd_start and read.reference_end > tsd_end: # span TSD
                    refcount += 1


    return altcount, refcount
 

def teVAF(bamfn, chrom, poslist, fam):
    alt, ref = break_count(bamfn, chrom, poslist)
    vaf = 'NA'

    if float(ref+alt) > 0:
        vaf = float(alt)/float(alt+ref)

    print '\t'.join(map(str, (bamfn, fam, chrom, min(poslist), max(poslist), alt, ref, vaf)))


def main(args):
    with open(args.tebreak, 'r') as inlist:
        for line in inlist:
            if line.startswith('UUID'):
                continue

            c = line.strip().split()

            chrom = c[1]
            poslist = [int(c[4]), int(c[5])]

            fam = c[7]

            if poslist[0] == poslist[1]:
                continue

            teVAF(args.bam, chrom, poslist, fam)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='genotyper')
    parser.add_argument('-b', '--bam', required='True')
    parser.add_argument('-t,', '--tebreak', required='True')
    args = parser.parse_args()
    main(args)
