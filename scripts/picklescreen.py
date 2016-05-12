#!/usr/bin/env python

import os
import re
import shutil
import cPickle as pickle
import argparse
import logging
import subprocess

from uuid import uuid4
from collections import defaultdict as dd


logger = logging.getLogger(__name__)


def count_alen(samrec):
    ''' return number of "M" bases in CIGAR string '''
    cigar = samrec[5]
    return float(sum([int(m[0]) for m in re.findall(r'(\d+)([M]{1})', cigar)]))


def count_mm(samrec):
    ''' parse NM flag for mismatches in SAM record, return -1 if no NM tag is present '''
    for field in samrec:
        if field.startswith('NM'):
            return float(field.split(':')[-1])

    return -1.


def matchpct(samrec):
    return 1.0-(count_mm(samrec)/count_alen(samrec))


def mapfilter(fq, ref, minscore=20, minmatch=0.95, threads=4):
    logger.debug('map %s to %s ...' % (fq, ref))

    bwa = ['bwa', 'mem', '-t', str(threads), '-M', '-Y', '-k', '10', '-P', '-S', '-T', str(minscore), ref, fq]

    keep = {}

    p = subprocess.Popen(bwa, stdout=subprocess.PIPE)

    for rec in p.stdout:
        if not rec.startswith('@'):
            rec = rec.split('\t')
            name = '-'.join(rec[0].split('-')[:-1])
            flag = int(rec[1])
            if bin(flag & 4) != bin(4): # not unmapped
                if matchpct(rec) > minmatch:
                    keep[name] = True

    return keep


def prepare_ref(fasta, refoutdir='tebreak_refs'):
    if not os.path.exists(refoutdir):
        os.mkdir(refoutdir)
        assert os.path.exists(refoutdir), 'could not create ref output directory: %s' % refoutdir

    ref_fa = refoutdir + '/' + os.path.basename(fasta)

    if not os.path.exists(ref_fa):
        logger.debug('Copying %s to %s ...' % (fasta, ref_fa))
        shutil.copy(fasta, ref_fa)

    logger.debug('Samtools indexing %s ...' % ref_fa)
    subprocess.call(['samtools', 'faidx', ref_fa])
    assert os.path.exists(ref_fa + '.fai'), 'could not samtools faidx %s' % ref_fa


    logger.debug('Create BWA db for %s ...' % ref_fa)
    p = subprocess.Popen(['bwa', 'index', ref_fa], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    for line in p.stdout: pass
    assert os.path.exists(ref_fa + '.bwt'), 'could not bwa index %s' % ref_fa

    return ref_fa


def makefq(insertions, tmpdir='/tmp', use_distal=False):
    ''' write FASTQ file from consensus sequences in insertions '''
    insfq = tmpdir + '/' + str(uuid4()) + '.fq'

    with open(insfq, 'w') as fq:
        for ins in insertions:
            for be in ('be1', 'be2'):
                if use_distal:
                    if be+'_dist_seq' in ins['INFO'] and ins['INFO'][be+'_dist_seq']:
                        fq.write('@%s\n%s\n+\n%s\n' % (ins['INFO']['ins_uuid']+'-'+be, ins['INFO'][be+'_dist_seq'], len(ins['INFO'][be+'_dist_seq'])*'B'))
                else:
                    if be+'_cons_seq' in ins['INFO'] and ins['INFO'][be+'_cons_seq']:
                        fq.write('@%s\n%s\n+\n%s\n' % (ins['INFO']['ins_uuid']+'-'+be, ins['INFO'][be+'_cons_seq'], len(ins['INFO'][be+'_cons_seq'])*'B'))

    return insfq


def main(args):
    logger.debug('loading pickle: %s' % args.pickle)

    with open(args.pickle, 'r') as pickin:
        insertions = pickle.load(pickin)

    logger.debug('finished loading %s' % args.pickle)
    logger.debug('raw candidate count: %d' % len(insertions))

    ref = prepare_ref(args.ref)

    fq = makefq(insertions, use_distal=args.use_distal)

    mapped = mapfilter(fq, ref, minscore=args.minscore, minmatch=float(args.minmatch), threads=int(args.threads))

    filtered = []

    for ins in insertions:
        keep = True
        if args.invert:
            if ins['INFO']['ins_uuid'] in mapped: keep = False
        else:
            if ins['INFO']['ins_uuid'] not in mapped: keep = False

        if keep:
            rec = dd(dict)
            rec['INFO'] = ins['INFO']
            rec['READSTORE'] = ins['READSTORE']
            filtered.append(rec)

    logger.debug('kept %d records' % len(filtered))

    with open(args.out, 'w') as pickout:
        pickle.dump(filtered, pickout)


if __name__ == '__main__':
    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(format=FORMAT)
    logger.setLevel(logging.DEBUG)

    parser = argparse.ArgumentParser(description='filter pickle based on alignments')
    parser.add_argument('-p', '--pickle', required=True, help='input filename (tebreak.py pickle)')
    parser.add_argument('-r', '--ref', required=True, help='reference to align consensus sequences against')
    parser.add_argument('-t', '--threads', default=4, help='alignment threads')
    parser.add_argument('-o', '--out', default='filtered.pickle', help='output filename (tebreak.py pickle)')
    parser.add_argument('-i', '--invert', default=False, action='store_true', help='retain insertions that do not match library')
    parser.add_argument('-s', '--minscore', default=20, help='minimum alignment score (-T option to bwa mem) default=20')
    parser.add_argument('-m', '--minmatch', default=0.95, help='minimum match pct/100 (default 0.95)')
    parser.add_argument('-d', '--use_distal', default=False, action='store_true', help='use only distal part of consensus for alignment')
    args = parser.parse_args()
    main(args)
