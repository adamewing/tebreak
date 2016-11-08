#!/usr/bin/env python

import sys
import os
import subprocess
import logging
import argparse

#import align
import pysam
import numpy as np

from uuid import uuid4


FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def rc(dna):
    ''' reverse complement '''
    complements = maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
    return dna.translate(complements)[::-1]


def align(qryseq, refseq, elt='PAIR', minmatch=90.0):
    rnd = str(uuid4())
    tgtfa = 'tmp.' + rnd + '.tgt.fa'
    qryfa = 'tmp.' + rnd + '.qry.fa'

    tgt = open(tgtfa, 'w')
    qry = open(qryfa, 'w')

    tgt.write('>ref' + '\n' + refseq + '\n')
    qry.write('>qry' + '\n' + qryseq + '\n')

    tgt.close()
    qry.close()

    cmd = ['exonerate', '--bestn', '1', '-m', 'ungapped', '--showalignment','0', '--ryo', elt + '\t%s\t%qab\t%qae\t%tab\t%tae\t%pi\t%qS\t%tS\n', qryfa, tgtfa]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    best = []
    topscore = 0

    for pline in p.stdout.readlines():
        if pline.startswith(elt):
            c = pline.strip().split()
            if int(c[1]) > topscore and float(c[6]) >= minmatch:
                topscore = int(c[1])
                best = c

    os.remove(tgtfa)
    os.remove(qryfa)

    return best


def flip_ends(rec):
    rec['5_Prime_End'], rec['3_Prime_End'] = rec['3_Prime_End'], rec['5_Prime_End']
    rec['Orient_5p'], rec['Orient_3p'] = rec['Orient_3p'], rec['Orient_5p']
    rec['5p_Elt_Match'], rec['3p_Elt_Match'] = rec['3p_Elt_Match'], rec['5p_Elt_Match']
    rec['5p_Genome_Match'], rec['3p_Genome_Match'] = rec['3p_Genome_Match'], rec['5p_Genome_Match']
    rec['Split_reads_5prime'], rec['Split_reads_3prime'] = rec['Split_reads_3prime'], rec['Split_reads_5prime']
    rec['5p_Cons_Len'], rec['3p_Cons_Len'] = rec['3p_Cons_Len'], rec['5p_Cons_Len']
    rec['5p_Improved'], rec['3p_Improved'] = rec['3p_Improved'], rec['5p_Improved']
    rec['TSD_3prime'], rec['TSD_5prime'] = rec['TSD_5prime'], rec['TSD_3prime']
    rec['Genomic_Consensus_5p'], rec['Genomic_Consensus_3p'] = rec['Genomic_Consensus_3p'], rec['Genomic_Consensus_5p']
    rec['Insert_Consensus_5p'], rec['Insert_Consensus_3p'] = rec['Insert_Consensus_3p'], rec['Insert_Consensus_5p']

    return rec


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


def main(args):

    inslib = None

    if args.insref:
        inslib = load_falib(args.insref)

    ref = pysam.Fastafile(args.refgenome)

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

                ins_id = '%s:%s' % (rec['Superfamily'], rec['Subfamily'])
                refseq = ref.fetch(rec['Chromosome'], int(rec['Left_Extreme']), int(rec['Right_Extreme']))

                #print rec['Genomic_Consensus_5p'], inslib[ins_id]
                elt_5p_align = align(rec['Genomic_Consensus_5p'], inslib[ins_id])
                elt_3p_align = align(rec['Genomic_Consensus_3p'], inslib[ins_id])
                gen_5p_align = align(rec['Genomic_Consensus_5p'], refseq)
                gen_3p_align = align(rec['Genomic_Consensus_3p'], refseq)

                elt_5p_orient = 'NA'
                elt_3p_orient = 'NA'
                gen_5p_orient = 'NA'
                gen_3p_orient = 'NA'

                if elt_5p_align:
                    elt_5p_orient = elt_5p_align[-1]

                if elt_3p_align:
                    elt_3p_orient = elt_3p_align[-1]

                if gen_5p_align:
                    gen_5p_orient = gen_5p_align[-1]

                if gen_3p_align:
                    gen_3p_orient = gen_3p_align[-1]

                new_5p_orient = 'NA'
                new_3p_orient = 'NA'

                if 'NA' not in (elt_5p_orient, gen_5p_orient):
                    if elt_5p_orient == gen_5p_orient:
                        new_5p_orient = '+'
                    else:
                        new_5p_orient = '-'

                if 'NA' not in (elt_3p_orient, gen_3p_orient):
                    if elt_3p_orient == gen_3p_orient:
                        new_3p_orient = '+'
                    else:
                        new_3p_orient = '-'

                coords_5p = []
                coords_3p = []

                if elt_5p_align:
                    coords_5p = sorted(map(int, (elt_5p_align[4], elt_5p_align[5])))

                if elt_3p_align:
                    coords_3p = sorted(map(int, (elt_3p_align[4], elt_3p_align[5])))

                #print elt_5p_align, coords_5p
                #print elt_3p_align, coords_3p

                flip = False
                if coords_5p and coords_3p and coords_5p[1] > coords_3p[1]:
                    flip = True

                #print flip

                rec['Orient_5p'] = new_5p_orient
                rec['Orient_3p'] = new_3p_orient

                if 'NA' not in (new_5p_orient, new_3p_orient) and 'None' not in (rec['Orient_5p'], rec['Orient_3p']):
                    if rec['Orient_5p'] != rec['Orient_3p']:
                        rec['Inversion'] = 'Y'
                    else:
                        rec['Inversion'] = 'N'

                else:
                    rec['Inversion'] = 'N'


                if flip:
                    rec = flip_ends(rec)

                out = [rec[h] for h in header]

                print '\t'.join(out)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='foo')
    parser.add_argument('-t', '--table', required=True, help='tabular output from resolve.py, requires header to be present')
    parser.add_argument('-i', '--insref', required=True, help='insertion sequence reference')
    parser.add_argument('-r', '--refgenome', required=True, help='reference genome')
    args = parser.parse_args()
    main(args)