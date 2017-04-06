#!/usr/bin/env python

import sys
import re
import pysam
import argparse

from string import maketrans
from random import betavariate, uniform, choice, randint


## 1691   3.0  0.0  0.0  chr10      495536  495734 (135039013) +  L1HS           LINE/L1               5957 6155    (0) 359265


class RefElt:
    __slots__ = ['chrom', 'start', 'end', 'top_strand', 'diverge', 'elt_type']

    def __init__(self, chrom, start, end, top_strand, diverge, elt_type, ref):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.top_strand = top_strand
        self.diverge = diverge
        self.elt_type = elt_type

        self.dna = self._seq(ref)


    def _seq(self, ref):
        if self.top_strand:
            return ref.fetch(self.chrom, self.start, self.end)
        else:
            return rc(ref.fetch(self.chrom, self.start, self.end))


    def notes(self):
        orig_strand = ('-', '+')[int(self.top_strand)]
        return 'type=%s,chrom=%s,orig_start=%d,orig_end=%d,orig_strand=%s,orig_div=%f' % (self.elt_type,self.chrom, self.start, self.end, orig_strand, self.diverge)


    def truncate(self, pct_fl=0.25):
        ''' simulate TPRT "U-shaped" distribution, return header note, seq '''

        if uniform(0,1) < pct_fl:
            trunc_start = randint(0,10)

        else:
            trunc_start = int(len(self.dna) * betavariate(6,1))

        note = 'trunc_start=%d,ins_len=%d' % (trunc_start,len(self.dna))

        return note, self.dna[trunc_start:len(self.dna)]


    def inversion(self):
        ''' simulate twin-priming inversion, return header note, seq '''
        offset = int(.3*len(self.dna)*betavariate(1.5,2))
        inv_pt = len(self.dna) - offset                        # inversion point
        t_start = int((len(self.dna)-offset)*betavariate(4,1)) # start relative to insertion ref
        delta = int(uniform(-10,10))                           # internal dup/del
        
        note = 'trunc_start=%d,inv_loc=%d,inv_slop=%d,ins_len=%d' % (t_start, inv_pt, delta, len(self.dna))

        end5p = self.dna[t_start:inv_pt]
        end3p = self.dna[inv_pt-delta:len(self.dna)]

        return note, rc(end5p) + end3p

    
    def transduction(self, ref, search=100):
        ''' simulate 3' transductions, return transduction only incl. polya signal '''
        signals = ['AATAAA', 'ATAAA', 'GTAAA']
        flank = ''

        if self.top_strand:
            flank = ref.fetch(self.chrom, self.end, self.end+search)

        else:
            flank = rc(ref.fetch(self.chrom, self.start - search, self.start))

        tr_seq = ''

        for signal in signals:
            if re.search(signal, flank):
                tr_seq = flank.split(signal)[0]+signal

        note = 'trd=%d' % len(tr_seq)

        return note, tr_seq


def rc(dna):
    ''' reverse complement '''
    complements = maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
    return dna.translate(complements)[::-1]


def pa(seq):
    ''' polyadenylate '''
    tail = 'A'*randint(10,100)
    return 'pa=%d' % len(tail), seq+tail


def main(args):

    ref = pysam.Fastafile(args.fasta)
    chr_prefix = ref.references[0].startswith('chr')

    elts = args.elts.split(',')

    ref_elts = []

    with open(args.rmsk, 'r') as rmsk:
        for line in rmsk:
            c = line.strip().split()

            if len(c) == 0 or c[0] in ('SW', 'score'):
                continue

            diverge = float(c[1])

            chrom = c[4]

            if not chr_prefix and chrom.startswith('chr'):
                chrom = chrom.replace('chr', '')

            if chr_prefix and not chrom.startswith('chr'):
                chrom = 'chr' + chrom

            start = int(c[5])
            end   = int(c[6])

            elt_start = int(c[11].strip('()'))
            elt_left = int(c[13].strip('()'))

            full_length = elt_left + elt_start < 150

            top_strand = c[8] == '+'

            elt_type = c[9]

            if full_length and diverge <= float(args.maxdiv) and elt_type in elts and chrom in ref.references:
                ref_elts.append(RefElt(chrom, start, end, top_strand, diverge, elt_type, ref))

    assert len(ref_elts) > 0, "no elements found!"

    for _ in range(int(args.num)):
        elt = choice(ref_elts)

        elt_seq = ''
        notes = elt.notes() # track insertion information

        # either invert+truncate or truncate
        if not args.notrunc:
            if uniform(0,1) <= float(args.inv):
                inv_note, elt_seq = elt.inversion()
                notes += ','+inv_note

            else:
                trunc_note, elt_seq = elt.truncate()
                notes += ','+trunc_note
        else:
            elt_seq = elt.dna


        # maybe add a transduction
        tr_seq = ''
        if uniform(0,1) <= float(args.tr): 
            tr_note, tr_seq = elt.transduction(ref, 1000)
            if tr_seq:
                elt_seq += tr_seq
                notes += ','+tr_note

        assert elt_seq, 'fail on elt: %s' % notes

        # polyadenylate
        if not args.nopolya:
            pa_note, elt_seq = pa(elt_seq)
            notes += ','+pa_note

        # insertion is minus strand 50% of the time
        if uniform(0,1) < .5:
            elt_seq = rc(elt_seq)
            notes += ','+'ins_strand=minus'

        else:
            notes += ','+'ins_strand=plus'

        notes += ','+'len=%d' % len(elt_seq)

        print '>%s\n%s' % (notes, elt_seq.upper())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='make insertion reference input for bamsurgeon/addsv')
    parser.add_argument('-f', '--fasta', required=True, help='reference genome FASTA (indexed with samtools faidx)')
    parser.add_argument('-r', '--rmsk', required=True, help='repeatmasker .out file')
    parser.add_argument('-e', '--elts', required=True, help='element type to grab from rmsk output')
    parser.add_argument('--maxdiv', default=0.5, help='maximum consensus divergence, default=0.5')
    parser.add_argument('-n', '--num', default=1000, help='number of elements to simulate')
    parser.add_argument('--tr', default=0.15, help='transduction chance, default=0.15')
    parser.add_argument('--inv', default=0.2, help='inversion chance, default=0.2')
    parser.add_argument('--notrunc', action='store_true', default=False, help='do not truncate insertions')
    parser.add_argument('--nopolya', action='store_true', default=False, help='do not polyadenylate')


    args = parser.parse_args()
    main(args)
