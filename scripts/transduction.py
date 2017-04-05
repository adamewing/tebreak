#!/usr/bin/env python

import os
import argparse
import subprocess
import logging
import sys

from uuid import uuid4
from string import maketrans
from collections import namedtuple
from collections import defaultdict as dd

import pysam
import numpy as np
from bx.intervals.intersection import Intersecter, Interval


FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

class Loc:
    def __init__(self, chrom, pos, strand):
        self.chrom = chrom
        self.pos = int(pos)
        self.strand = strand
        self.annot = []

    def __str__(self):
        return '%s:%d:%s:%s' % (self.chrom, self.pos, self.strand, ';'.join(self.annot))


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


def align(qryseq, refseq):
    rnd = str(uuid4())
    tgtfa = 'tmp.' + rnd + '.tgt.fa'
    qryfa = 'tmp.' + rnd + '.qry.fa'

    tgt = open(tgtfa, 'w')
    qry = open(qryfa, 'w')

    tgt.write('>ref' + '\n' + refseq + '\n')
    qry.write('>qry' + '\n' + qryseq + '\n')

    tgt.close()
    qry.close()

    cmd = ['exonerate', '--bestn', '1', '-m', 'ungapped', '--showalignment', '1', '--ryo', 'ALN' + '\t%s\t%qab\t%qae\t%tab\t%tae\t%tS\t%pi\n', qryfa, tgtfa]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    best = []
    topscore = 0

    for pline in p.stdout.readlines():
        #print "***|" + pline.strip()
        if pline.startswith('ALN'):
            c = pline.strip().split()
            if int(c[1]) > topscore:
                topscore = int(c[1])
                best = c
        #else:
        #    print pline.rstrip()

    os.remove(tgtfa)
    os.remove(qryfa)

    return best


def rc(dna):
    ''' reverse complement '''
    complements = maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
    return dna.translate(complements)[::-1]


def polya_scan(seq):
    max_polya = []
    cur_polya = []

    for i, b in enumerate(seq.upper()):
        if b == 'A':
            cur_polya.append(i)
        else:
            if len(cur_polya) > len(max_polya):
                max_polya = cur_polya
            cur_polya = []

    if len(cur_polya) > len(max_polya):
        max_polya = cur_polya

    if len(max_polya) > 0:
        return min(max_polya), max(max_polya)
    else:
        return 0,0


def guess_forward(seq):
    nA = len([b for b in list(seq.upper()) if b == 'A'])
    nT = len([b for b in list(seq.upper()) if b == 'T'])

    if nA > nT:
        return True

    return False


def get_trn_seq(cons_seq, ins_seq, ref_seq, minlen=20):
    if not guess_forward(cons_seq):
        cons_seq = rc(cons_seq)

    ins_align = align(cons_seq, ins_seq)
    gen_align = align(cons_seq, ref_seq)

    #print 'INS ' + str(ins_align)
    #print 'GEN ' + str(gen_align)

    if ins_align:
        ins_start, ins_end = sorted(map(int, (ins_align[2], ins_align[3])))

    trn_seq = None

    gen_start = len(cons_seq)
    if gen_align:
        gen_start, gen_end = sorted(map(int, (gen_align[2], gen_align[3])))

    if ins_align:
        trd_start = ins_end
        trd_end   = gen_start

        if trd_end - trd_start > minlen:
            trn_seq = cons_seq[trd_start:trd_end]

    if trn_seq is None:
        cons_cover = np.zeros(len(cons_seq))
        if ins_align:
            for i in range(ins_start, ins_end):
                cons_cover[i] = 1

        if gen_align:
            for i in range(gen_start, gen_end):
                cons_cover[i] = 1

        bases = [cons_seq[i] for i in range(len(cons_seq)) if cons_cover[i] == 0]
    
        if len(bases) >= minlen:
            trn_seq = ''.join(bases)
    
    return trn_seq


def locate(seq, refgenome):
    rnd = str(uuid4())
    qryfq = 'tmp.' + rnd + '.qry.fq'

    with open(qryfq, 'w') as qry:
        qry.write('@query' + '\n' + seq + '\n+\n' + 'B'*len(seq) + '\n')

    FNULL = open(os.devnull, 'w')

    sam_cmd  = ['bwa', 'mem', refgenome, qryfq]
    view_cmd = ['samtools', 'view', '-F', '0x4', '-'] # view mapped reads only
    aln  = subprocess.Popen(sam_cmd, stdout=subprocess.PIPE, stderr=FNULL)
    view = subprocess.Popen(view_cmd, stdin=aln.stdout, stdout=subprocess.PIPE, stderr=FNULL)

    locs = []

    for line in view.stdout:
        c = line.strip().split()

        bits  = int(c[1])
        chrom = c[2]
        pos   = int(c[3])

        strand = '-'
        if bits & 0x10 == 0:
            strand = '+'

        locs.append(Loc(chrom, pos, strand))


    os.remove(qryfq)

    return locs


def main(args):

    inslib = load_falib(args.insref)
    ref = pysam.Fastafile(args.refgenome)

    logger.info("loading bwa index %s into shared memory ..." % args.refgenome)
    p = subprocess.Popen(['bwa', 'shm', args.refgenome], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    for line in p.stdout: pass # wait for bwa to load
    logger.info("loaded.")

    uuid_forest = dd(Intersecter)

    with open(args.tabfile, 'r') as tab:
        for line in tab:
            if line.startswith('UUID'):
                continue

            c = line.strip().split()
            uuid, chrom, start, end = c[:4]
            start = int(start)
            end = int(end)

            uuid_forest[chrom].add_interval(Interval(start, end, value=uuid))


    header = []
    with open(args.tabfile, 'r') as tab:
        for i, line in enumerate(tab):
            if i == 0: # header
                header = line.strip().split('\t')
                header += ['Truth', 'Transduction']

            else:
                rec = {}
                for n, field in enumerate(line.strip().split('\t')):
                    rec[header[n]] = field

                left  = int(rec['Left_Extreme'])-1000
                right = int(rec['Right_Extreme'])+1000

                if left < 0: left = 0

                ref_seq = ref.fetch(rec['Chromosome'], left, right)

                seqn = rec['Superfamily'] + ':' + rec['Subfamily']

                ins_seq = inslib[seqn].rstrip('Aa') # trim ins ref polyA if present

                locs = []

                if rec['Insert_Consensus_3p'] != 'NA':
                    trn_seq = get_trn_seq(rec['Insert_Consensus_3p'], ins_seq, ref_seq)
                    if trn_seq:
                        #print trn_seq
                        locs += locate(trn_seq, args.refgenome)

                if rec['Genomic_Consensus_3p'] != 'NA':
                    trn_seq = get_trn_seq(rec['Genomic_Consensus_3p'], ins_seq, ref_seq)
                    if trn_seq:
                        #print trn_seq
                        locs += locate(trn_seq, args.refgenome)
                
                # annotate reference transductions
                if args.reftranslib is not None:
                    annot_locs = []
                    reflib = pysam.Tabixfile(args.reftranslib)
                    for loc in locs:
                        if loc.chrom in reflib.contigs:
                            for refelt in reflib.fetch(loc.chrom, loc.pos-100, loc.pos+100):
                                c = refelt.split()
                                refstart = int(c[1])
                                refend = int(c[2])
                                refstr = c[3]
                                #print c, loc.pos, loc.strand
                                if refstr == '+' and loc.pos > refend-20:
                                    loc.annot.append('REF:' + ':'.join(refelt.split()))

                                elif refstr == '-' and loc.pos < refstart+20:
                                    loc.annot.append('REF:' + ':'.join(refelt.split()))

                # annotate non-reference transductions
                for loc in locs:
                    if loc.chrom in uuid_forest:
                        for nr_elt in uuid_forest[loc.chrom].find(loc.pos-100, loc.pos+100):
                            if nr_elt.value != rec['UUID']:
                                loc.annot.append('NR:' + nr_elt.value)

                locstr = []
                if locs:
                    for loc in locs:
                        locstr.append(str(loc))
                else:
                    locstr.append('NA')

                tru_string = []
                if args.truth is not None:
                    tru = pysam.Tabixfile(args.truth)
                    
                    if rec['Chromosome'] in tru.contigs:
                        for vcfrec in tru.fetch(rec['Chromosome'], int(rec['Left_Extreme']), int(rec['Right_Extreme'])):
                            tru_string.append(vcfrec.strip().split()[2])
                else:
                    tru_string.append('NA')

                print line.strip() + '\t' + ','.join(tru_string) + '\t' + ','.join(locstr)




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='find transductions')
    parser.add_argument('-i', '--insref', required=True)
    parser.add_argument('-r', '--refgenome', required=True)
    parser.add_argument('-t', '--tabfile', required=True)
    parser.add_argument('--reftranslib', default=None)
    parser.add_argument('--truth', default=None, help='debug')
    args = parser.parse_args()
    main(args)