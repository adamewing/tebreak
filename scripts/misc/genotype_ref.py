#!/usr/bin/env python

import argparse
import csv
import os

import pysam


def break_count(bam, chrom, poslist, minpad=5, flex=1, minmapq=10):
    ''' ref = number of reads spanning TSD, alt = number of reads clipped at breakpoint in poslist '''
    altcount = 0
    refcount = 0
    discards = 0

    poslist = list(poslist)

    tsd_start = min(poslist)
    tsd_end   = max(poslist)

    tsd_len = tsd_end - tsd_start

    if tsd_start < minpad: tsd_start = minpad

    for read in bam.fetch(chrom, tsd_start-minpad, tsd_end+minpad):
        if read.is_unmapped or read.is_duplicate:
            continue

        if read.mapq < minmapq:
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


def getVAF(bam, chrom, poslist):
    ''' return number of reads supporting alt (insertion), ref (reference) and vaf (variant allele fraction) '''
    poslist = map(int, poslist)
    alt, ref = break_count(bam, chrom, poslist)
    vaf = 0.0 

    if float(ref+alt) > 0:
        vaf = float(alt)/float(alt+ref)

    return alt, ref, vaf


def main(args):

    print('##fileformat=VCFv4.1')

    vcf_cols = ['#CHROM', 'POS ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']

    bams = []
    with open(args.bamlist) as _:
        for bam in _:
            fn, name = bam.strip().split()
            assert os.path.exists(fn.strip())
            bams.append(fn)
            vcf_cols.append(name)

    print('\t'.join(vcf_cols))

    with open(args.table) as table:
        csv_reader = csv.DictReader(table, delimiter='\t')

        for rec in csv_reader:
            tsd_5p = list(map(int, [rec['TSD_Start_5p'], rec['TSD_End_5p']]))
            tsd_3p = list(map(int, [rec['TSD_Start_3p'], rec['TSD_End_3p']]))

            info = 'ELT=%s;ORIENT=%s' % (rec['Name'], rec['Orientation'])

            vcf_line = [rec['Chromosome'], str(rec['TSD_Start_5p']), '.', rec['TSD_seq'][0], '<INS>', '100', 'PASS', info, 'GT:DS']

            for bam_fn in bams:
                bam = pysam.AlignmentFile(bam_fn)

                vaf_5p = getVAF(bam, rec['Chromosome'], tsd_5p)
                vaf_3p = getVAF(bam, rec['Chromosome'], tsd_3p)

                alt_total = vaf_5p[0] + vaf_3p[0]
                ref_total = vaf_5p[1] + vaf_3p[1]
                vaf_total = float(alt_total)/float(alt_total+ref_total)

                dose = 0.0
                gt = './.'

                if alt_total + ref_total >= int(args.mindepth):
                    dose = 2-(vaf_total*2) # ref-specific

                    gt = '1/1' # default to homz. reference for insertions in ref assembly

                    if dose > float(args.hetlow)*2 and dose < float(args.hethi)*2:
                        gt = '0/1'

                    if dose < float(args.hetlow)*2: # ref-specific
                        gt = '0/0'

                fmt = '%s:%.3f' % (gt, dose)

                vcf_line.append(fmt)
                
                bam.close()

            print('\t'.join(vcf_line))
                

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='genotype reference insertions')
    parser.add_argument('-b', '--bamlist', required=True)
    parser.add_argument('-t', '--table', required=True)
    parser.add_argument('--mindepth', default=10)
    parser.add_argument('--hetlow', default=0.15)
    parser.add_argument('--hethi', default=0.85)
    args = parser.parse_args()
    main(args)
