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
import pandas as pd

from collections import defaultdict as dd

from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC

import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

import seaborn as sns
sns.set()

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

    data = dd(dict)

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

                knowngood = 0
                if rec['NonRef'] != 'NA':
                    knowngood = 1

                uuid = rec['UUID']

                if 'NA' in (rec['TE_Align_Start'], rec['TE_Align_End']):
                    continue

                if min(float(rec['5p_Genome_Match']), float(rec['3p_Genome_Match'])) == 0.0:
                    continue

                if min(float(rec['5p_Elt_Match']), float(rec['3p_Elt_Match'])) == 0.0:
                    continue
                
                data[uuid]['Superfamily'] = rec['Superfamily']
                data[uuid]['KnownGood'] = knowngood

                
                data[uuid]['BigWindowSize']        = abs(int(rec['Left_Extreme']) - int(rec['Right_Extreme']))
                data[uuid]['SmallWindowSize']      = abs(int(rec['5_Prime_End']) - int(rec['3_Prime_End']))
                data[uuid]['TE_Align_Start']       = int(rec['TE_Align_Start'])
                data[uuid]['TE_Align_End']         = int(rec['TE_Align_End'])
                data[uuid]['5p_Elt_Match']         = float(rec['5p_Elt_Match'])
                data[uuid]['3p_Elt_Match']         = float(rec['3p_Elt_Match'])
                data[uuid]['5p_Genome_Match']      = float(rec['5p_Genome_Match'])
                data[uuid]['3p_Genome_Match']      = float(rec['3p_Genome_Match'])
                data[uuid]['Log_Split_reads']      = np.log(int(rec['Split_reads_5prime']) + int(rec['Split_reads_3prime']))
                data[uuid]['Remapped_Discordant']  = float(rec['Remapped_Discordant'])
                data[uuid]['Remap_Disc_Fraction']  = float(rec['Remap_Disc_Fraction'])
                data[uuid]['Remapped_Splitreads']  = float(rec['Remapped_Splitreads'])
                data[uuid]['Remap_Split_Fraction'] = float(rec['Remap_Split_Fraction'])
                data[uuid]['5p_Cons_Len']          = int(rec['5p_Cons_Len'])
                data[uuid]['3p_Cons_Len']          = int(rec['3p_Cons_Len'])
                data[uuid]['TSD_Len']              = len(rec['TSD_3prime'])
                data[uuid]['Varcount']             = len(rec['Variants'].split(','))

                data[uuid]['Mapscore'] = avgmap(map_tbx, rec['Chromosome'], rec['Left_Extreme'], rec['Right_Extreme'])

                if rec['TSD_3prime'] == 'NA': data[uuid]['TSD_Len'] = 0
                if rec['Variants'] == 'NA': data[uuid]['Varcount'] = 0




        data = pd.DataFrame.from_dict(data, orient='index')


        vplot_labels = ['BigWindowSize',
                        'SmallWindowSize',
                        '5p_Elt_Match',
                        '3p_Elt_Match',
                        '5p_Genome_Match',
                        '3p_Genome_Match',
                        'Log_Split_reads',
                        'Remapped_Discordant',
                        'Remap_Disc_Fraction',
                        'Remapped_Splitreads',
                        'Remap_Split_Fraction',
                        '5p_Cons_Len',
                        '3p_Cons_Len',
                        'TSD_Len',
                        'TE_Align_Start',
                        'TE_Align_End',
                        'Varcount',
                        'Mapscore']

        sns_plot = {}
        fig = {}

        # for elt in ['L1', 'ALU', 'SVA']:
        #     for vplot_label in vplot_labels:
        #         sns.plt.clf()
        #         logger.info('plotting %s' % vplot_label)
        #         sns_plot = sns.violinplot(x="KnownGood", y=vplot_label, data=data[data.Superfamily==elt], jitter=True)
        #         fig = sns_plot.figure
        #         fig.savefig(args.outbase+elt+'.'+vplot_label+'.png')


        for elt in ['L1', 'ALU', 'SVA']:
            clf = RandomForestClassifier()
            #clf = LinearDiscriminantAnalysis()
            #clf = AdaBoostClassifier()
            #clf = GaussianNB()
            #clf = KNeighborsClassifier()
            #clf = SVC(gamma=2, C=1)

            elt_data = data[data.Superfamily==elt]

            x = elt_data[vplot_labels]
            y = elt_data['KnownGood']

            clf.fit(x,y)
            score = clf.score(x,y)
            sys.stderr.write('%s score: %f\n' % (elt, score))

            pred = clf.predict(x)

            elt_data['Prediction'] = pred

            pred_data = elt_data[elt_data.Prediction==1]
            #print pred_data

            with open(args.tabfile, 'r') as tab:
                for i, line in enumerate(tab):

                    if i == 0: # header
                        header = line.strip().split('\t')
                        #print line.strip()

                    else:
                        line = line.strip()
                        uuid = line.split()[0]

                        if uuid in pred_data.index:
                            print line

            #data.to_csv('test.tsv', sep='\t')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='filter script for TEs on hg19')
    parser.add_argument('--tabfile', required=True, help='tabular output from resolve.py, requires header to be present')
    parser.add_argument('--insref', default=None, help='req. alignment to insertion sequence reference')
    parser.add_argument('--ignore_ref_filter', default=False, action='store_true', help='turn of filtering vs. reference elements')
    parser.add_argument('--outbase', default='')
    args = parser.parse_args()
    main(args)
