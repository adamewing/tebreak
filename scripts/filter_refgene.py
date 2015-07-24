#!/usr/bin/env python

import os
import sys
import pysam


def usage():
    return 'usage: %s </path/to/TEBreak directory> <tabular output from resolve.py>' % sys.argv[0]


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


if len(sys.argv) == 3:
    tebreak_dir = sys.argv[1]

    if not os.path.exists(tebreak_dir):
        sys.exit(usage())

    map_ref = tebreak_dir + '/lib/wgEncodeCrgMapabilityAlign50mer.bed.gz'
    pgo_ref = tebreak_dir + '/lib/PGO_Build74.coords.bed.gz'

    map_tbx = pysam.Tabixfile(map_ref)
    pgo_tbx = pysam.Tabixfile(pgo_ref)

    header = []
    with open(sys.argv[2], 'r') as tab:
        for i, line in enumerate(tab):

            if i == 0: # header
                header = line.strip().split('\t')
                print line.strip()

            else:
                rec = {}
                out = True
                for n, field in enumerate(line.strip().split('\t')):
                    rec[header[n]] = field

                if int(rec['3p_Cons_Len']) < 120 and int(rec['5p_Cons_Len']) < 120: out = False

                #if 'Y' not in (rec['5p_Improved'], rec['3p_Improved']): out = False

                if 'NA' in (rec['TE_Align_Start'], rec['TE_Align_End']):
                    out = False
                else:
                    if int(rec['TE_Align_End']) - int(rec['TE_Align_Start']) < 400:
                        out = False

                if rec['TSD_3prime'] != rec['TSD_5prime']: out = False

                if min(int(rec['Split_reads_5prime']), int(rec['Split_reads_3prime'])) < 2: out = False
                if float(rec['Remap_Disc_Fraction']) < 0.5 : out = False
                if max(float(rec['5p_Elt_Match']), float(rec['3p_Elt_Match'])) < 0.95: out = False
                if min(float(rec['5p_Elt_Match']), float(rec['3p_Elt_Match'])) < 0.90: out = False
                if max(float(rec['5p_Genome_Match']), float(rec['3p_Genome_Match'])) < 0.98: out = False
                if min(float(rec['5p_Genome_Match']), float(rec['3p_Genome_Match'])) < 0.95: out = False
                if int(rec['Remapped_Discordant']) < 8 : out = False

                if out:
                    if rec['Chromosome'] in pgo_tbx.contigs:
                        if len(list(pgo_tbx.fetch(rec['Chromosome'], int(rec['Left_Extreme']), int(rec['Right_Extreme'])))) > 0:
                            out = False

                if out:
                    if avgmap(map_tbx, rec['Chromosome'], rec['Left_Extreme'], rec['Right_Extreme']) < 0.9: out = False

                if out: print line.strip()


else:
    sys.exit(usage())
