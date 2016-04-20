#!/usr/bin/env python

''' script for filtering insertions vs. the human reference GRCh37/hg19 '''
''' may be useful as a template for extension to other species          '''


import pysam
import sys
import os
import logging
import argparse

verbose=True

FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
if verbose:
    logger.setLevel(logging.DEBUG)
else:
    logger.setLevel(logging.INFO)



def usage():
    return 'usage: %s </path/to/TEBreak directory> <tabular output from resolve.py>' % sys.argv[0]


def ref_filter(chrom, start, end, superfams, tbx):
    for sf in superfams.split(','):
        if sf == 'L1':
            for ins in tbx[sf].fetch(chrom, int(start), int(end)): return True

        if sf in ('ALU', 'SVA'):
            for ins in tbx[sf].fetch(chrom, int(start), int(end)): return True

        if sf == 'SVA':
            for ins in tbx[sf].fetch(chrom, int(start), int(end)): return True

    return False


def len_filter(rec):
    telen = int(rec['TE_Align_End']) - int(rec['TE_Align_Start'])
    if 'ALU' in rec['Superfamily'] and telen < 250: return True
    if 'SVA' in rec['Superfamily'] and telen < 1000: return True
    if 'L1' in rec['Superfamily'] and int(rec['TE_Align_End']) < 5950: return True

    return False


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
    tebreak_dir = args.tebreak

    if not os.path.exists(tebreak_dir):
        sys.exit(usage())

    l1_ref  = tebreak_dir + '/lib/mask.L1.hg19.bed.gz'
    alu_ref = tebreak_dir + '/lib/mask.Alu.hg19.bed.gz'
    sva_ref = tebreak_dir + '/lib/mask.SVA.hg19.bed.gz'
    map_ref = tebreak_dir + '/lib/wgEncodeCrgMapabilityAlign100mer.bed.gz'

    for fn in (l1_ref, alu_ref, sva_ref):
        if not os.path.exists(fn): sys.exit('reference %s not found' % fn)
        if not os.path.exists(fn + '.tbi'): sys.exit('index for reference %s not found' %fn)

    tbx = {}
    tbx['L1']  = pysam.Tabixfile(l1_ref)
    tbx['ALU'] = pysam.Tabixfile(alu_ref)
    tbx['SVA'] = pysam.Tabixfile(sva_ref)

    map_tbx = pysam.Tabixfile(map_ref)

    header = []
    with open(args.tabfile, 'r') as tab:
        for i, line in enumerate(tab):

            if i == 0: # header
                header = line.strip().split('\t')
                print line.strip()

            else:
                rec = {}
                out = True
                for n, field in enumerate(line.strip().split('\t')):
                    rec[header[n]] = field

                #logger.debug(rec['UUID'])

                if int(rec['3p_Cons_Len']) < 120 and int(rec['5p_Cons_Len']) < 120:
                    logger.debug('Filtered %s: consensus length < %d' % (rec['UUID'], 120))
                    out = False

                if 'NA' in (rec['TE_Align_Start'], rec['TE_Align_End']):
                    logger.debug('Filtered %s: TE_Align_Start or TE_Align_End is "NA"' % rec['UUID'])
                    out = False

                if ref_filter(rec['Chromosome'], rec['Left_Extreme'], rec['Right_Extreme'], rec['Superfamily'], tbx) and not args.ignore_ref_filter:
                    logger.debug('Filtered %s: proximity to reference TE of same superfamily' % rec['UUID']) 
                    out = False

                if max(float(rec['5p_Elt_Match']), float(rec['3p_Elt_Match'])) < 0.95:
                    logger.debug('Filtered %s: max(5p_Elt_Match, 3p_Elt_Match) < 0.95' % rec['UUID'])
                    out = False

                if max(float(rec['5p_Genome_Match']), float(rec['3p_Genome_Match'])) < 0.98:
                    logger.debug('Filtered %s: max(5p_Genome_Match, 3p_Genome_Match) < 0.98' % rec['UUID'])
                    out = False

                mapscore = avgmap(map_tbx, rec['Chromosome'], rec['Left_Extreme'], rec['Right_Extreme']) * (max(int(rec['3p_Cons_Len']), int(rec['5p_Cons_Len']))/100.)
                if mapscore < 0.1:
                    logger.debug('Filtered %s: mappability of %f < 0.1' % (rec['UUID'], mapscore))
                    out = False

                if float(rec['Remapped_Discordant']) < 4 or float(rec['Remap_Disc_Fraction']) < 0.5:
                    logger.debug('Filtered %s: low discordant evidence (< 4 reads or < 50pct supporting)' % rec['UUID'])
                    out = False

                if out and len_filter(rec):
                    logger.debug('Filtered %s: TE length filter' % rec['UUID'])
                    out = False

                if out: print line.strip()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='filter script for TEs on hg19')
    parser.add_argument('--tabfile', required=True, help='tabular output from resolve.py, requires header to be present')
    parser.add_argument('--tebreak', required=True, help='path to tebreak repo directory (needed to find annotation files)')
    parser.add_argument('--ignore_ref_filter', default=False, action='store_true', help='turn of filtering vs. reference elements')
    args = parser.parse_args()
    main(args)
