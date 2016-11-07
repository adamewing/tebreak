#!/usr/bin/env python

import argparse
import logging

import networkx as nx
import numpy as np

from uuid import uuid4
from collections import Counter
from collections import defaultdict as dd
from bx.intervals.intersection import Intersecter, Interval # pip install bx-python



verbose=False

FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
if verbose:
    logger.setLevel(logging.DEBUG)
else:
    logger.setLevel(logging.INFO)


def getrec(fn):
    ''' generator yields records '''
    header = []
    with open(fn, 'r') as table:
        for i, line in enumerate(table):
            if i == 0:
                header = line.strip().split('\t')

            else:
                rec = {}

                for n, field in enumerate(line.strip().split('\t')):
                    rec[header[n]] = field

                yield rec


def make_uuid_sets(tables):
    ''' UUID interval tree --> UUID intersection graph --> connected components --> UUID sets'''
    forest = dd(Intersecter)

    for table in tables:
        for rec in getrec(table):
            chrom = rec['Chromosome']
            start = int(rec['5_Prime_End'])-10
            end   = int(rec['3_Prime_End'])+10

            if start > end: start, end = end, start

            forest[chrom].add_interval(Interval(start, end, value=rec['UUID']))

    G = nx.Graph()

    for table in tables:
        for rec in getrec(table):
            chrom = rec['Chromosome']
            start = int(rec['5_Prime_End'])-10
            end   = int(rec['3_Prime_End'])+10

            if start > end: start, end = end, start

            for i in forest[chrom].find(start, end):
                G.add_edge(rec['UUID'], i.value)

    return list(nx.connected_components(G))


def make_rec_table(tables):
    ''' make two level dict: master --> UUID --> column '''
    master = {}

    for table in tables:
        for rec in getrec(table):
            assert rec['UUID'] not in master, 'UUID collision! %s' % rec['UUID']

            master[rec['UUID']] = rec

    return master


def combine(uuids, recs):

    combo = {}

    uncertainty_warning = []

    use_any = ['Chromosome']
    use_max = ['Right_Extreme', '5p_Cons_Len', '3p_Cons_Len']
    use_min = ['Left_Extreme']

    use_maj = ['5_Prime_End',
            '3_Prime_End',
            'Superfamily',
            'Subfamily',
            'Orient_5p',
            'Orient_3p',
            'Inversion',
            'TSD_3prime',
            'TSD_5prime']

    use_sum = ['Split_reads_5prime',
            'Split_reads_3prime',
            'Remapped_Discordant',
            'Remapped_Splitreads',
            'Sample_count']

    use_med_int = ['TE_Align_Start',
            'TE_Align_End']

    use_med_flt = ['5p_Elt_Match',
            '3p_Elt_Match',
            '5p_Genome_Match',
            '3p_Genome_Match',
            'Remap_Disc_Fraction',
            'Remap_Split_Fraction']

    use_YN = ['5p_Improved', '3p_Improved']

    use_union = ['Sample_support', 'Variants']

    use_len = ['Genomic_Consensus_5p',
            'Genomic_Consensus_3p',
            'Insert_Consensus_5p',
            'Insert_Consensus_3p']


    combo['UUID'] = str(uuid4())

    for col in use_any:
        combo[col] = recs[list(uuids)[0]][col]

    for col in use_max:
        combo[col] = max([int(recs[uuid][col]) for uuid in uuids])

    for col in use_min:
        combo[col] = min([int(recs[uuid][col]) for uuid in uuids])

    for col in use_sum:
        combo[col] = sum([int(recs[uuid][col]) for uuid in uuids])

    for col in use_med_int:
        combo[col] = int(np.median([int(recs[uuid][col]) for uuid in uuids]))

    for col in use_med_flt:
        combo[col] = np.median([float(recs[uuid][col]) for uuid in uuids])

    for col in use_len:
        seqs = [recs[uuid][col] for uuid in uuids]
        seqs.sort(key=len)
        combo[col] = seqs[-1]

    for col in use_YN:
        yn = [recs[uuid][col] for uuid in uuids]
        if 'Y' in yn:
            combo[col] = 'Y'
        else:
            combo[col] = 'N'

    for col in use_union:
        union = []
        for uuid in uuids:
            union += recs[uuid][col].split(',')

        combo[col] = ','.join(union)

    for col in use_maj:
        maj = Counter([recs[uuid][col] for uuid in uuids]).most_common(1)[0]
        combo[col] = maj[0]

        if maj[1] < len(uuids):
            uncertainty_warning.append(col)

    combo['Old_UUIDs'] = ','.join(uuids)

    if uncertainty_warning:
        combo['CombineWarnings'] = ','.join(uncertainty_warning)

    else:
        combo['CombineWarnings'] = 'NA'

    return combo

    # for c in combo:
    #     print c, ':', combo[c]
    # print '*'*100


def main(args):
    logger.info('building UUID sets...')

    uuid_sets = make_uuid_sets(args.tables)

    logger.info('built %d UUID sets' % len(uuid_sets))

    logger.info('building master record table...')

    recs = make_rec_table(args.tables)

    logger.info('master record table has %d entires' % len(recs))

    combo_recs = [combine(u, recs) for u in uuid_sets]

    header = ['UUID',
    'Chromosome',
    'Left_Extreme',
    'Right_Extreme',
    '5_Prime_End',
    '3_Prime_End',
    'Superfamily',
    'Subfamily',
    'TE_Align_Start',
    'TE_Align_End',
    'Orient_5p',
    'Orient_3p',
    'Inversion',
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
    '5p_Improved',
    '3p_Improved',
    'TSD_3prime',
    'TSD_5prime',
    'Sample_count',
    'Sample_support',
    'Genomic_Consensus_5p',
    'Genomic_Consensus_3p',
    'Insert_Consensus_5p',
    'Insert_Consensus_3p',
    'Variants',
    'Old_UUIDs',
    'CombineWarnings']

    print '\t'.join(header)

    for rec in combo_recs:
        out = [str(rec[h]) for h in header]

        print '\t'.join(out)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='merge tables in a relatively sane way')
    parser.add_argument('tables', nargs='+')
    args = parser.parse_args()
    main(args)