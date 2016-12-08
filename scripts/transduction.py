#!/usr/bin/env python

import os
import argparse
import subprocess
import logging

from uuid import uuid4
from string import maketrans
from collections import namedtuple

import pysam

FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


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


def get_trn_seq(cons_seq, ins_seq, ref_seq, minlen=20):
    ins_align = align(cons_seq, ins_seq)
    gen_align = align(cons_seq, ref_seq)

    #print 'INS ' + str(ins_align)
    #print 'GEN ' + str(gen_align)

    if ins_align and gen_align:
        ins_start, ins_end = sorted(map(int, (ins_align[2], ins_align[3])))
        gen_start, gen_end = sorted(map(int, (gen_align[2], gen_align[3])))

        trd_start = min(ins_end, gen_end)
        trd_end   = max(ins_start, gen_start)

        if trd_end - trd_start > minlen:
            return cons_seq[trd_start:trd_end]

    return None


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

    Loc = namedtuple('Loc', 'chrom pos strand')
    locs = []

    for line in view.stdout:
        c = line.strip().split()

        bits  = int(c[1])
        chrom = c[2]
        pos   = int(c[3])

        strand = '-'
        if bits & 0x10 == 0:
            strand = '+'

        locs.append(Loc(chrom=chrom, pos=pos, strand=strand))


    os.remove(qryfq)

    return locs


def main(args):

    inslib = load_falib(args.insref)
    ref = pysam.Fastafile(args.refgenome)

    logger.info("loading bwa index %s into shared memory ..." % args.refgenome)
    p = subprocess.Popen(['bwa', 'shm', args.refgenome], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    for line in p.stdout: pass # wait for bwa to load
    logger.info("loaded.")

    header = []
    with open(args.tabfile, 'r') as tab:
        for i, line in enumerate(tab):
            if i == 0: # header
                header = line.strip().split('\t')
                header.append('Transduction')

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
                
                locstr = []
                if locs:
                    for loc in locs:
                        locstr.append(':'.join(map(str, loc)))
                else:
                    locstr.append('NA')

                print line.strip() + '\t' + ','.join(locstr)




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='find transductions')
    parser.add_argument('-i', '--insref', required=True)
    parser.add_argument('-r', '--refgenome', required=True)
    parser.add_argument('-t', '--tabfile', required=True)
    args = parser.parse_args()
    main(args)