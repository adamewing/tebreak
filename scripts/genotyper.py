#!/usr/bin/env python


import os
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
 

def getVAF(bamfn, chrom, poslist):
    poslist = map(int, poslist)
    alt, ref = break_count(bamfn, chrom, poslist)
    vaf = 0.0 

    if float(ref+alt) > 0:
        vaf = float(alt)/float(alt+ref)

    return alt, ref, vaf


def genotype(runlist):

    with open(args.tabfile, 'r') as tab:
        for i, line in enumerate(tab):
            if i == 0: # header
                header = line.strip().split('\t')

                for name, _ in runlist:
                    header += [name+'_RefCount', name+'_Alt', name+'_VAF']

                print '\t'.join(header)
                continue

            rec = {}
            out = True
            for n, field in enumerate(line.strip().split('\t')):
                rec[header[n]] = field

            c = line.strip().split()

            for _, bamfn in runlist:
                alt, ref, vaf = getVAF(bamfn, rec['Chromosome'], (rec['5_Prime_End'], rec['3_Prime_End']))

                c.append(str(ref))
                c.append(str(alt))
                c.append(str(vaf))

            print '\t'.join(c)


def main(args):

    runlist = []

    with open(args.bamlist, 'r') as bamlist:
        for line in bamlist:
            name, bamfn = line.strip().split()

            runlist.append((name, bamfn))

    genotype(runlist)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='genotyper')
    parser.add_argument('-b', '--bamlist', required='True')
    parser.add_argument('-t', '--tabfile', required='True')
    args = parser.parse_args()
    main(args)
