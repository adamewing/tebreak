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


def hompol_scan(seq, base):
    ''' return start and end coords of longest homopolymer run'''
    max_hompol = []
    cur_hompol = []

    for i, b in enumerate(seq.upper()):
        if b == base:
            cur_hompol.append(i)
        else:
            if len(cur_hompol) > len(max_hompol):
                max_hompol = cur_hompol
            cur_hompol = []

    if len(cur_hompol) > len(max_hompol):
        max_hompol = cur_hompol

    if len(max_hompol) > 0:
        return min(max_hompol), max(max_hompol)-min(max_hompol)
    else:
        return 0,0


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
    #print ' '.join(cmd)

    return best


def flip_ends(rec):
    rec['5_Prime_End'], rec['3_Prime_End'] = rec['3_Prime_End'], rec['5_Prime_End']
    rec['Orient_5p'], rec['Orient_3p'] = rec['Orient_3p'], rec['Orient_5p']
    rec['5p_Elt_Match'], rec['3p_Elt_Match'] = rec['3p_Elt_Match'], rec['5p_Elt_Match']
    rec['5p_Genome_Match'], rec['3p_Genome_Match'] = rec['3p_Genome_Match'], rec['5p_Genome_Match']
    rec['Split_reads_5prime'], rec['Split_reads_3prime'] = rec['Split_reads_3prime'], rec['Split_reads_5prime']
    if 'Sample_support_5p' in rec:
        rec['Sample_support_5p'], rec['Sample_support_3p'] = rec['Sample_support_3p'], rec['Sample_support_5p']
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


def fix_ins_id(ins_id, inslib):
    superfam, subfam = ins_id.split(':')

    for i in inslib.keys():
        if i.split(':')[-1] == subfam:
            superfam = i.split(':')[0]

    return '%s:%s' % (superfam, subfam)


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

    inslib = None

    if args.insref:
        inslib = load_falib(args.insref)

    ref = pysam.Fastafile(args.refgenome)

    header = []

    count_5p_diff = 0
    count_3p_diff = 0
    count_5p_switchcons = 0
    count_3p_switchcons = 0

    if args.maptabix:
        maptabix = pysam.Tabixfile(args.maptabix)

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

                if rec['Superfamily'] == 'NA':
                    ins_id = rec['Subfamily']

                if rec['Subfamily'] == 'NA':
                    ins_id = rec['Superfamily']

                if ins_id not in inslib:
                    if ':' not in ins_id:
                        logger.warn('No insertion identification for %s (ins_id %s)' % (rec['UUID'], ins_id))
                        continue
                    else:
                        ins_id = fix_ins_id(ins_id, inslib)

                    if ins_id not in inslib:
                        logger.warn('No insertion identification for %s (ins_id %s)' % (rec['UUID'], ins_id))
                        continue

                # filtering here

                out = True

                if rec['Insert_Consensus_5p'] == rec['Insert_Consensus_3p'] == 'NA':
                    logger.debug('Filtered %s: no insertion consensus mapped to insertion reference' % rec['UUID'])
                    out = False

                if int(rec['3p_Cons_Len']) + int(rec['5p_Cons_Len']) < int(args.conslen):
                    logger.info('Filtered %s: total consensus length < %d' % (rec['UUID'], int(args.conslen)))
                    out = False

                if max(float(rec['5p_Elt_Match']), float(rec['3p_Elt_Match'])) < float(args.eltmatch):
                    logger.info('Filtered %s: max(5p_Elt_Match, 3p_Elt_Match) < %f' % (rec['UUID'], float(args.eltmatch)))
                    out = False

                if max(float(rec['5p_Genome_Match']), float(rec['3p_Genome_Match'])) < float(args.refmatch):
                    logger.info('Filtered %s: max(5p_Genome_Match, 3p_Genome_Match) < %f' % (rec['UUID'], float(args.refmatch)))
                    out = False

                if float(rec['Remapped_Discordant']) < int(args.numdiscord):
                    logger.info('Filtered %s: low discordant evidence (< %d reads)' % (rec['UUID'], int(args.numdiscord)))
                    out = False

                if int(rec['Split_reads_5prime']) + int(rec['Split_reads_3prime']) < int(args.numsplit):
                    logger.info('Filtered %s: low split read evidence (< %d reads)' % (rec['UUID'], int(args.numsplit)))
                    out = False

                if levenshtein(rec['TSD_3prime'], rec['TSD_5prime']) > 1:
                    logger.info('Filtered %s: TSD mismatch: %s vs %s' % (rec['UUID'], rec['TSD_5prime'], rec['TSD_3prime']))
                    out = False

                elif (len(list(set(list(rec['TSD_3prime'])))) == 1 or len(list(set(list(rec['TSD_5prime'])))) == 1) and len(rec['TSD_3prime']) > 4:
                    logger.info('Filtered %s: TSD is a long homopolymer: %s' % (rec['UUID'], rec['TSD_3prime']))
                    out = False

                if args.fracend is not None and rec['Genomic_Consensus_3p'] != 'NA' and rec['TE_Align_End'] != 'NA':
                    fracend = float(len(inslib[ins_id]) - int(rec['TE_Align_End'])) / float(len(inslib[ins_id]))
                    if fracend > float(args.fracend):
                        out = False

                if args.maxvars is not None and 'Variants' in rec:
                    numvars = len(rec['Variants'].split(','))
                    if numvars > int(args.maxvars):
                        out = False

                if out:
                    if args.tabix is not None:
                        for posfilter in [pysam.Tabixfile(fn) for fn in args.tabix.split(',')]:
                            if rec['Chromosome'] in posfilter.contigs:
                                if len(list(posfilter.fetch(rec['Chromosome'], int(rec['Left_Extreme']), int(rec['Right_Extreme'])))) > 0:
                                    out = False

                    if args.refgene is not None:
                        gene_tbx = pysam.Tabixfile(args.refgene)
                        if rec['Chromosome'] in gene_tbx.contigs:
                            gene_overlaps = []
                            for overlap in gene_tbx.fetch(rec['Chromosome'], int(rec['Left_Extreme']), int(rec['Right_Extreme'])):
                                gene_overlaps.append(overlap.split()[3])

                            gene = '.'.join(rec['Superfamily'].split('.')[:-1])
                            if gene in gene_overlaps:
                                out = False

                
                if out and args.maptabix:
                    mapscore = avgmap(maptabix, rec['Chromosome'], int(rec['Left_Extreme']), int(rec['Right_Extreme']))
                    if mapscore < 0.5:
                        logger.info('Filtered %s: mappability low' % rec['UUID'])
                        out = False


                if out:
                    left = min(int(rec['5_Prime_End']), int(rec['3_Prime_End'])) - 20
                    right = max(int(rec['5_Prime_End']), int(rec['3_Prime_End'])) + 20
                    refseq = ref.fetch(rec['Chromosome'], left, right)

                    for b in ('A', 'T', 'C', 'G'):
                        hp = hompol_scan(refseq, b)
                        #print hp, refseq
                        if hp[1] > 20:
                            logger.info('Filtered %s: homopolymer in insertion site' % rec['UUID'])
                            out = False

                if out:
                    refseq = ref.fetch(rec['Chromosome'], int(rec['Left_Extreme'])-100, int(rec['Right_Extreme'])+100)

                    # check that the reference genome region doesn't have a good match to the reference element sequence
                    # inslib[ins_id] vs refseq

                    self_align = align(inslib[ins_id], refseq, minmatch=85)

                    if self_align:
                        logger.info('Filtered %s: self-match between genome and refelt: %s' % (rec['UUID'], str(self_align)))
                        out = False

                if out:
                    # realignment

                    refseq = ref.fetch(rec['Chromosome'], int(rec['Left_Extreme']), int(rec['Right_Extreme']))

                    elt_5p_align = align(rec['Genomic_Consensus_5p'], inslib[ins_id])
                    elt_3p_align = align(rec['Genomic_Consensus_3p'], inslib[ins_id])
                    gen_5p_align = align(rec['Genomic_Consensus_5p'], refseq)
                    gen_3p_align = align(rec['Genomic_Consensus_3p'], refseq)

                    if args.realign_all_isoforms:
                        align_store = {'elt_5p': elt_5p_align, 'elt_3p': elt_3p_align, 'gen_5p': gen_5p_align, 'gen_3p': gen_3p_align}

                        for align_result in align_store:

                            basename = '.'.join(ins_id.split('.')[:-1])
                            tried_isoforms = [ins_id]
                            isoform_results = []

                            if align_store[align_result]:
                                isoform_results.append(align_store[align_result])

                            isoform_names = [tx for tx in inslib.keys() if '.'.join(tx.split('.')[:-1]) == basename]

                            for isoname in isoform_names:
                                if isoname != ins_id:

                                    targetseq = refseq
                                    if align_result.startswith('elt'):
                                        targetseq = inslib[isoname]

                                    isores = align(rec['Genomic_Consensus_' + align_result.split('_')[1]], targetseq)

                                    if isores:
                                        isoform_results.append(isores)

                            if isoform_results:
                                align_store[align_result] = sorted(isoform_results, key=lambda isores: isores[1])[0]

                        elt_5p_align = align_store['elt_5p']
                        elt_3p_align = align_store['elt_3p']
                        gen_5p_align = align_store['gen_5p']
                        gen_3p_align = align_store['gen_3p']

                    # try using the insertion-based consensus if no luck with the genomic one

                    if not elt_5p_align or not gen_5p_align:
                        retry_elt_5p_align = align(rec['Insert_Consensus_5p'], inslib[ins_id])
                        retry_gen_5p_align = align(rec['Insert_Consensus_5p'], refseq)

                        if retry_gen_5p_align and retry_elt_5p_align:
                            elt_5p_align = retry_elt_5p_align
                            gen_5p_align = retry_gen_5p_align
                            count_5p_switchcons += 1


                    if not elt_3p_align or not gen_3p_align:
                        retry_elt_3p_align = align(rec['Insert_Consensus_3p'], inslib[ins_id])
                        retry_gen_3p_align = align(rec['Insert_Consensus_3p'], refseq)

                        if retry_gen_3p_align and retry_elt_3p_align:
                            elt_3p_align = retry_elt_3p_align
                            gen_3p_align = retry_gen_3p_align
                            count_3p_switchcons += 1

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

                    flip = False
                    if coords_5p and coords_3p and coords_5p[1] > coords_3p[1]:
                        flip = True

                    if rec['Orient_5p'] != new_5p_orient:
                        logger.info('Changed 5p orientation for %s' % rec['UUID'])
                        count_5p_diff += 1

                    if rec['Orient_3p'] != new_3p_orient:
                        logger.info('Changed 3p orientation for %s' % rec['UUID'])
                        count_3p_diff += 1

                    rec['Orient_5p'] = new_5p_orient
                    rec['Orient_3p'] = new_3p_orient

                    if new_3p_orient == new_5p_orient == 'NA':
                        out = False

                    if args.require_realign:
                        if 'NA' in (rec['Orient_5p'], rec['Orient_3p']):
                            out = False

                    if 'NA' not in (new_5p_orient, new_3p_orient) and 'None' not in (rec['Orient_5p'], rec['Orient_3p']):
                        if rec['Orient_5p'] != rec['Orient_3p']:
                            rec['Inversion'] = 'Y'
                        else:
                            rec['Inversion'] = 'N'

                    else:
                        rec['Inversion'] = 'N'


                    if flip:
                        rec = flip_ends(rec)

                    if out:
                        print '\t'.join([rec[h] for h in header])


    logger.info('Changed orientation on %d 5p ends' % count_5p_diff)
    logger.info('Changed orientation on %d 3p ends' % count_3p_diff)
    logger.info('Used insertion consensus for %d 5p ends' % count_5p_switchcons)
    logger.info('Used insertion consensus for %d 3p ends' % count_3p_switchcons)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='foo')
    parser.add_argument('-t', '--table', required=True, help='tabular output from resolve.py, requires header to be present')
    parser.add_argument('-i', '--insref', required=True, help='insertion sequence reference')
    parser.add_argument('-r', '--refgenome', required=True, help='reference genome')

    parser.add_argument('--conslen', default=300, help='min total consensus length (default=300)')
    parser.add_argument('--eltmatch', default=0.95, help='min element match (default=0.95)')
    parser.add_argument('--refmatch', default=0.98, help='min reference match (default=0.98)')
    parser.add_argument('--numsplit', default=4, help='min number of supporting split read mappings (default=4)')
    parser.add_argument('--numdiscord', default=4, help='min number of supporting discordant read mappings (default=4)')
    parser.add_argument('--fracend', default=None, help='filter insertions by distance to reference 3p end by fraction of length')
    parser.add_argument('--maxvars', default=None, help='maximum number of variants in Variants column')

    parser.add_argument('--tabix', default=None, help='tabix files to be used as positional filters (can be comma-delimited list of tabix-indexed files)')
    parser.add_argument('--refgene', default=None, help='filter out insertions with superfamily overlapping matching refgene (important for GRIP searches)')
    parser.add_argument('--maptabix', default=None, help='use mappability tabix')

    parser.add_argument('--realign_all_isoforms', action='store_true', help='enable for potentially better realignment when working on GRIPs')
    parser.add_argument('--require_realign', action='store_true', help='require that all breakends are realignable')
    args = parser.parse_args()
    main(args)
