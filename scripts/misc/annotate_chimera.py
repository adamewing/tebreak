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


    header = []

    with open(args.tabfile, 'r') as tab:
        for i, line in enumerate(tab):

            if i == 0: # header
                header = line.strip().split('\t')
                header += ['ChimeraBaseCount', 'ChimeraMatchIns', 'ChimeraMatchRef', 'InsSiteHomology', 'PossibleRefEltChimera']

                print '\t'.join(header)

            else:
                rec = {}
                for n, field in enumerate(line.strip().split('\t')):
                    rec[header[n]] = field

                ins_site_homlen = 0 # insertion site homology length
                ins_site_homseq = 'NA' # sequence of overlapped region
                ch_ref_present = False
                ins_pct_match = 0.0
                ref_pct_match = 0.0


                ref = pysam.Fastafile(args.refgenome)

                left  = int(rec['Left_Extreme']) - 1000
                right = int(rec['Right_Extreme']) + 1000

                if left < 0: left = 0

                ref_seq = ref.fetch(rec['Chromosome'], left, right)

                seqn = rec['Superfamily'] + ':' + rec['Subfamily']

                if 'NA' in (rec['Superfamily'], rec['Subfamily']):
                    continue

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

                    ch_align_ins = align(ins_site_homseq, ins_seq, 'Ins')
                    ch_align_ref = align(ins_site_homseq, ref_seq, 'Ref')

                    if ch_align_ins:
                        ins_pct_match = ch_align_ins[-1]
                    if ch_align_ref:
                        ref_pct_match = ch_align_ref[-1]

                # chimera with adjacent ref element check
                ch_ref_present = ref_filter(rec['Chromosome'], rec['Left_Extreme'], rec['Right_Extreme'], rec['Superfamily'], tbx, extend=10000)

                # output
                fields = line.strip().split()
                fields.append(str(ins_site_homlen))
                fields.append(str(ins_pct_match))
                fields.append(str(ref_pct_match))
                fields.append(ins_site_homseq)
                fields.append(str(ch_ref_present))
                
                print '\t'.join(fields)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='foo')
    parser.add_argument('-t', '--tabfile', required=True, help='tabular output from resolve.py, requires header to be present')
    parser.add_argument('-i', '--insref', required=True, help='insertion sequence reference')
    parser.add_argument('-r', '--refgenome', required=True, help='reference genome')
    args = parser.parse_args()
    main(args)
