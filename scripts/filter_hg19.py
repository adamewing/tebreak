#!/usr/bin/env python

''' script for filtering insertions vs. the human reference GRCh37/hg19 '''
''' may be useful as a template for extension to other species          '''


import pysam
import sys
import os
import logging
import argparse
import align
import numpy as np

verbose=False

FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
if verbose:
    logger.setLevel(logging.DEBUG)
else:
    logger.setLevel(logging.INFO)

tebreak_dir = os.path.dirname(os.path.realpath(__file__))


def load_falib(infa):
    seqdict = {}

    with open(infa, 'r') as fa:
        seqid = ''
        seq   = ''
        for line in fa:
            if line.startswith('>'):
                if seq != '':
                    seqdict[seqid] = seq
                seqid = line.lstrip('>').strip().split()[0]
                seq   = ''
            else:
                assert seqid != ''
                seq = seq + line.strip()

    if seqid not in seqdict and seq != '':
        seqdict[seqid] = seq

    return seqdict


def ref_filter(chrom, start, end, superfams, tbx):
    for sf in superfams.split(','):
        if sf == 'L1':
            if chrom not in tbx[sf].contigs: return True
            for ins in tbx[sf].fetch(chrom, int(start), int(end)): return True

        if sf in ('ALU', 'SVA'):
            if chrom not in tbx[sf].contigs: return True
            for ins in tbx[sf].fetch(chrom, int(start), int(end)): return True

        if sf == 'SVA':
            if chrom not in tbx[sf].contigs: return True
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


def rc(dna):
    ''' reverse complement '''
    complements = maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
    return dna.translate(complements)[::-1]


def realign_filter(rec, inslib):
    S = -np.ones((256, 256)) + 2 * np.identity(256)
    S = S.astype(np.int16)

    seqn = rec['Superfamily'] + ':' + rec['Subfamily']
    if seqn not in inslib:
        return False

    seq_headers = ['Genomic_Consensus_5p', 'Genomic_Consensus_3p', 'Insert_Consensus_5p', 'Insert_Consensus_3p']

    for seqtype in seq_headers:
        s1 = align.string_to_alignment(rec[seqtype])
        s2 = align.string_to_alignment(inslib[seqn])

        (s, a1, a2) = align.align(s1, s2, -2, -2, S, local=True)
        a1 = align.alignment_to_string(a1)
        a2 = ''.join([b for b in list(align.alignment_to_string(a2)) if b != '-'])

        score = 0.0
        if len(a1) > 0:
            score = float(len(a1) - (len(a1)-s)) / float(len(a1))

        #print seqtype, score, len(a1)

        if score > 0.9 and len(a1) > 25:
            return False

        return True


def main(args):

    l1_ref  = tebreak_dir + '/../lib/mask.L1.hg19.bed.gz'
    alu_ref = tebreak_dir + '/../lib/mask.Alu.hg19.bed.gz'
    sva_ref = tebreak_dir + '/../lib/mask.SVA.hg19.bed.gz'
    map_ref = tebreak_dir + '/../lib/wgEncodeCrgMapabilityAlign100mer.bed.gz'

    inslib = None

    if args.insref:
        inslib = load_falib(args.insref)


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

                if rec['Insert_Consensus_5p'] == rec['Insert_Consensus_3p'] == 'NA':
                    logger.debug('Filtered %s: no insertion consensus mapped to insertion reference' % rec['UUID'])
                    out = False

                #if out and len_filter(rec):
                #    logger.debug('Filtered %s: TE length filter' % rec['UUID'])
                #    out = False

                if args.insref:
                    if realign_filter(rec, inslib): out = False

                if out: print line.strip()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='filter script for TEs on hg19')
    parser.add_argument('--tabfile', required=True, help='tabular output from resolve.py, requires header to be present')
    parser.add_argument('--insref', default=None, help='req. alignment to insertion sequence reference')
    parser.add_argument('--ignore_ref_filter', default=False, action='store_true', help='turn of filtering vs. reference elements')
    args = parser.parse_args()
    main(args)
