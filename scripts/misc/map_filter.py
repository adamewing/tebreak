#!/usr/bin/env python

import pysam
import argparse


def avgmap(maptabix, chrom, start, end):
    ''' return average mappability across chrom:start-end region; maptabix = pysam.Tabixfile'''
    scores = []

    if None in (start, end): return None

    if chrom in maptabix.contigs:
        for rec in maptabix.fetch(chrom, int(start), int(end)):
            mchrom, mstart, mend, mscore = rec.strip().split()
            mstart, mend = int(mstart), int(mend)
            mscore = float(mscore)

            while mstart < mend and mstart:
                mstart += 1
                if mstart >= int(start) and mstart <= int(end):
                    scores.append(mscore)

        if len(scores) > 0:
            return sum(scores) / float(len(scores))
        else:
            return 0.0
    else:
        return 0.0


def main(args):
    maptabix = pysam.Tabixfile(args.maptabix)
    cutoff = float(args.cutoff)

    header = []

    with open(args.table, 'r') as table:
        for i, line in enumerate(table):
            if i == 0:
                header = line.strip().split('\t')
                print line.strip()

            else:
                rec = {}

                for n, field in enumerate(line.strip().split('\t')):
                    rec[header[n]] = field

                mapscore = avgmap(maptabix, rec['Chromosome'], int(rec['Left_Extreme']), int(rec['Right_Extreme']))

                if mapscore > cutoff:
                    print line.strip()




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='foo')
    parser.add_argument('-t', '--table', required=True, help='TEBreak table')
    parser.add_argument('-m', '--maptabix', required=True, help='use mappability tabix')
    parser.add_argument('-c', '--cutoff', default=0.5, help='mappability cutoff (default = 0.5)')
    args = parser.parse_args()
    main(args)

