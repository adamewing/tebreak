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
import subprocess

from uuid import uuid4

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


def ref_filter(chrom, start, end, superfams, tbx, extend=0):
    start = int(start)
    end   = int(end)

    if extend > 0:
        start -= extend
        end   += extend

        if start < 0: start = 0


    for sf in superfams.split(','):
        if sf == 'L1':
            if chrom not in tbx[sf].contigs: return True
            for ins in tbx[sf].fetch(chrom, start, end): return True

        if sf in ('ALU', 'SVA'):
            if chrom not in tbx[sf].contigs: return True
            for ins in tbx[sf].fetch(chrom, start, end): return True

        if sf == 'SVA':
            if chrom not in tbx[sf].contigs: return True
            for ins in tbx[sf].fetch(chrom, start, end): return True

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
    seqn = rec['Superfamily'] + ':' + rec['Subfamily']
    if seqn not in inslib:
        return False

    seq_headers = ['Genomic_Consensus_5p', 'Genomic_Consensus_3p', 'Insert_Consensus_5p', 'Insert_Consensus_3p']

    matches = []

    for seqtype in seq_headers:
        if rec[seqtype] == 'NA':
            continue

        #print seqtype, rec[seqtype]

        alignment = align(rec[seqtype], inslib[seqn], rec['Subfamily'])

        if alignment:
            matches.append([seqtype] + alignment)

    return matches


def align(qryseq, refseq, elt):
    rnd = str(uuid4())
    tgtfa = 'tmp.' + rnd + '.tgt.fa'
    qryfa = 'tmp.' + rnd + '.qry.fa'

    tgt = open(tgtfa, 'w')
    qry = open(qryfa, 'w')

    tgt.write('>ref' + '\n' + refseq + '\n')
    qry.write('>qry' + '\n' + qryseq + '\n')

    tgt.close()
    qry.close()

    cmd = ['exonerate', '--bestn', '1', '-m', 'ungapped', '--showalignment','0', '--ryo', elt + '\t%s\t%qab\t%qae\t%tab\t%tae\t%pi\n', qryfa, tgtfa]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    best = []
    topscore = 0

    for pline in p.stdout.readlines():
        if pline.startswith(elt):
            c = pline.strip().split()
            if int(c[1]) > topscore:
                topscore = int(c[1])
                best = c

    os.remove(tgtfa)
    os.remove(qryfa)

    return best


def overlap(iv1, iv2): 
    if min(iv1[1], iv2[1]) - max(iv1[0], iv2[0]) > 0: # is there overlap?
        return [max(iv1[0], iv2[0]), min(iv1[1], iv2[1])]
 
    return None


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

                if args.insref:
                    header += ['Mappability', 'ExonerateRealign']

                if args.chimera:
                    header += ['ChimeraBaseCount', 'InsSiteHomology']

                print '\t'.join(header)

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

                ref_present = False

                if args.wideref:
                    ref_present = ref_filter(rec['Chromosome'], rec['Left_Extreme'], rec['Right_Extreme'], rec['Superfamily'], tbx, extend=10000)
                else:
                    ref_present = ref_filter(rec['Chromosome'], rec['Left_Extreme'], rec['Right_Extreme'], rec['Superfamily'], tbx)

                if ref_present and not args.ignore_ref_filter:
                    logger.debug('Filtered %s: proximity to reference TE of same superfamily' % rec['UUID']) 
                    out = False

                if max(float(rec['5p_Elt_Match']), float(rec['3p_Elt_Match'])) < 0.95:
                    logger.debug('Filtered %s: max(5p_Elt_Match, 3p_Elt_Match) < 0.95' % rec['UUID'])
                    out = False

                if max(float(rec['5p_Genome_Match']), float(rec['3p_Genome_Match'])) < 0.98:
                    logger.debug('Filtered %s: max(5p_Genome_Match, 3p_Genome_Match) < 0.98' % rec['UUID'])
                    out = False

                mapscore = avgmap(map_tbx, rec['Chromosome'], rec['Left_Extreme'], rec['Right_Extreme'])# * (max(int(rec['3p_Cons_Len']), int(rec['5p_Cons_Len']))/100.)
                if mapscore < 0.1:
                    logger.debug('Filtered %s: mappability of %f < 0.1' % (rec['UUID'], mapscore))
                    out = False

                if float(rec['Remapped_Discordant']) < 4:
                    logger.debug('Filtered %s: low discordant evidence (< 4 reads)' % rec['UUID'])
                    out = False

                if float(rec['Remap_Disc_Fraction']) < 0.5:
                    logger.debug('Filtered %s: low discordant evidence (< 50pct supporting)' % rec['UUID'])
                    out = False

                if rec['Insert_Consensus_5p'] == rec['Insert_Consensus_3p'] == 'NA':
                    logger.debug('Filtered %s: no insertion consensus mapped to insertion reference' % rec['UUID'])
                    out = False

                if args.lenfilter and out and len_filter(rec):
                    logger.debug('Filtered %s: TE length filter' % rec['UUID'])
                    out = False

                align_info = 'NA'

                if out and args.insref and 'ExonerateRealign' not in header:
                    align_info = realign_filter(rec, inslib)

                    if len(align_info) == 0:
                        out = False

                    well_aligned = False
                    for alignment in align_info:
                        seqtype, _, score, qstart, qend, tstart, tend, pi = alignment
                        tstart = int(tstart)
                        tend   = int(tend)
                        pi     = float(pi)

                        if pi >= 95.0 and abs(tend-tstart) >= 100:
                            well_aligned = True

                    if not well_aligned: out = False


                ins_site_homlen = 0 # insertion site homology length
                ins_site_homseq = 'NA' # sequence of overlapped region

                if out and args.chimera:
                    if not args.refgenome:
                        sys.exit('--refgenome required in conjunction with --chimera')

                    ref = pysam.Fastafile(args.refgenome)

                    left  = int(rec['Left_Extreme']) - 1000
                    right = int(rec['Right_Extreme']) + 1000

                    if left < 0: left = 0

                    ref_seq = ref.fetch(rec['Chromosome'], left, right)

                    seqn = rec['Superfamily'] + ':' + rec['Subfamily']

                    ins_seq = inslib[seqn]

                    alignside = ''

                    ins_align = []
                    gen_align = []

                    if rec['Genomic_Consensus_3p'] != 'NA':
                        ins_align = align(rec['Genomic_Consensus_3p'], ins_seq, rec['Subfamily'])
                        gen_align = align(rec['Genomic_Consensus_3p'], ref_seq, 'Genomic')
                        alignside = 'Genomic_Consensus_3p'

                    else:
                        ins_align = align(rec['Genomic_Consensus_5p'], ins_seq, rec['Subfamily'])
                        gen_align = align(rec['Genomic_Consensus_5p'], ref_seq, 'Genomic')
                        alignside = 'Genomic_Consensus_5p'

                    ins_subcoords = None

                    if ins_align:
                        ins_subcoords = map(int, ins_align[2:4])

                    gen_subcoords = None

                    if gen_align:
                        gen_subcoords = map(int, gen_align[2:4])
                    else:
                        out = False

                    ol = None

                    if gen_subcoords is not None and ins_subcoords is not None:
                        ol = overlap(ins_subcoords, gen_subcoords)

                    if ol is not None:
                        ins_site_homlen = ol[1]-ol[0]
                        ins_site_homseq = rec[alignside][ol[0]:ol[1]]

                if out:
                    fields = line.strip().split()
                    fields.append(str(mapscore))
                    fields.append(','.join([';'.join(alignment) for alignment in align_info]))

                    if args.chimera:
                        fields.append(str(ins_site_homlen))
                        fields.append(ins_site_homseq)
                    
                    print '\t'.join(fields)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='filter script for TEs on hg19')
    parser.add_argument('--tabfile', required=True, help='tabular output from resolve.py, requires header to be present')
    parser.add_argument('--insref', default=None, help='req. alignment to insertion sequence reference')
    parser.add_argument('--ignore_ref_filter', default=False, action='store_true', help='turn of filtering vs. reference elements')
    parser.add_argument('--lenfilter', default=False, action='store_true', help='turn on filter by insertion length')
    parser.add_argument('--refgenome', default=None)
    parser.add_argument('--chimera', default=False, action='store_true')
    parser.add_argument('--wideref', default=False, action='store_true')
    args = parser.parse_args()
    main(args)
