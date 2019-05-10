#!/usr/bin/env python3

import argparse

import pandas as pd
pd.options.mode.chained_assignment = None

import numpy as np
import scipy.stats as ss

from sklearn.ensemble import IsolationForest

import logging
FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


model_cols = [
'TE_Align_Start',
'TE_Align_End',
'5p_Elt_Match',
'3p_Elt_Match',
'5p_Genome_Match',
'3p_Genome_Match',
'Split_reads_5prime',
'Split_reads_3prime',
'Remapped_Discordant',
'Remap_Disc_Fraction',
'Remapped_Splitreads',
'Remap_Split_Fraction',
'5p_Cons_Len',
'3p_Cons_Len',
'TSD_Len',
'TSD_Dist',
'TSD_Bases',
'5p_Improved',
'3p_Improved',
'Max_VAF',
'Max_Depth',
'Elt_Len'
]


def levenshtein(s1, s2):
    ''' from https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#Python '''
    if len(s1) < len(s2):
        return levenshtein(s2, s1)

    # len(s1) >= len(s2)
    if len(s2) == 0:
        return len(s1)

    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1 # j+1 instead of j since previous_row and current_row are one character longer
            deletions = current_row[j] + 1       # than s2
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row
    
    return previous_row[-1]


def tsd_len(t):
    t = list(map(str, t))
    out = []
    for tsd in t:
        if tsd == 'nan':
            out.append(0)
        else:
            out.append(len(tsd))

    return out


def tsd_dist(data):
    d = []
    for uuid in data.index:
        d.append(levenshtein(str(data.loc[uuid]['TSD_5prime']), str(data.loc[uuid]['TSD_3prime'])))

    return d

def elt_len(data):
    l = []
    for uuid in data.index:
        l.append(data.loc[uuid]['TE_Align_End']-data.loc[uuid]['TE_Align_Start'])

    return l

def nr_count(t):
    t = list(map(str, t))
    out = []
    for nr in t:
        if nr == 'nan':
            out.append(0)
        else:
            out.append(len(nr.split('|')[-1].split(',')))

    return out

def tsd_bases(t):
    t = list(map(str, t))
    out = []
    for tsd in t:
        if tsd == 'nan':
            out.append(0)
        else:
            out.append(len(set(list(tsd))))

    return out


def max_gt_vaf(t):
    t = list(map(str, t))
    out = []
    for gt in t:
        if gt == 'nan':
            out.append(0.0)
        else:
            out.append(max([float(g.split('|')[-1]) for g in gt.split(',')]))

    return out

def max_gt_depth(t):
    t = list(map(str, t))
    out = []
    for gt in t:
        if gt == 'nan':
            out.append(0.0)
        else:
            out.append(max([int(g.split('|')[1])+int(g.split('|')[2]) for g in gt.split(',')]))

    return out


def main(args):
    data = pd.read_csv(args.table, sep='\t', header=0, index_col=0)
    orig = pd.read_csv(args.table, sep='\t', header=0, index_col=0, keep_default_na=False, na_values=['_'])

    orig['pred'] = 0

    # modify and generate additional columns
    data['Split_reads_5prime'] = data['Split_reads_5prime']/data['Sample_count']
    data['Split_reads_3prime'] = data['Split_reads_3prime']/data['Sample_count']
    data['Remapped_Discordant'] = data['Remapped_Discordant']/data['Sample_count']
    data['Remapped_Splitreads'] = data['Remapped_Splitreads']/data['Sample_count']

    data['5p_Improved'] = list(map(int, data['5p_Improved']=='Y'))
    data['3p_Improved'] = list(map(int, data['3p_Improved']=='Y'))

    data['TSD_Len'] = tsd_len(data['TSD_3prime'])
    data['TSD_Dist'] = tsd_dist(data)
    data['TSD_Bases'] = tsd_bases(data['TSD_3prime'])
    data['Max_VAF'] = max_gt_vaf(data['Genotypes'])
    data['Max_Depth'] = max_gt_depth(data['Genotypes'])
    data['Elt_Len'] = elt_len(data)

    data['NonRef'] = nr_count(data['NonRef'])

    subsets = list(set(data['Superfamily']))
    X_train = data.loc[data['NonRef'] > 1]


    for subset in subsets:
        
        X_train_subset = X_train[X_train['Superfamily'] == subset]
        X_test_subset  = data[data['Superfamily'] == subset]

        orig_subset = orig[orig['Superfamily'] == subset]

        X_train_subset = X_train_subset[model_cols]
        X_test_subset = X_test_subset[model_cols]

        logger.info("Training %s" % subset)

        clf = IsolationForest(behaviour='new', contamination='auto', random_state=42, max_samples='auto')
        clf.fit(X_train_subset)

        y_pred_train = clf.predict(X_train_subset)
        y_pred_test = clf.predict(X_test_subset)

        orig['pred'].loc[orig_subset.index] = y_pred_test

    # positive examples should be positive
    #orig['pred'].loc[orig['NonRef'] != 'NA'] = 1
    orig['pred'].loc[data['NonRef'] > 1] = 1


    orig.to_csv('%s.isoforest.txt' % args.table, sep='\t')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='giant bucket')
    parser.add_argument('-t', '--table', required=True, help='tebreak table')

    args = parser.parse_args()
    main(args)
