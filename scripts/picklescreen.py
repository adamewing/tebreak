#!/usr/bin/env python

import os
import shutil
import cPickle as pickle
import argparse
import logging
import subprocess

from uuid import uuid4
from collections import defaultdict as dd


logger = logging.getLogger(__name__)


def map(fq, ref, threads=4):
    logger.debug('map %s to %s ...' % (fq, ref))

    bwa = ['bwa', 'mem', '-t', str(threads), '-M', '-Y', '-k', '10', '-P', '-S', '-T', '20', ref, fq]

    keep = {}

    p = subprocess.Popen(bwa, stdout=subprocess.PIPE)

    for rec in p.stdout:
        if not rec.startswith('@'):
            rec = rec.split('\t')
            name = '-'.join(rec[0].split('-')[:-1])
            flag = int(rec[1])
            if bin(flag & 4) != bin(4): # not unmapped
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


def makefq(insertions, tmpdir='/tmp'):
    ''' write FASTQ file from consensus sequences in insertions '''
    insfq = tmpdir + '/' + str(uuid4()) + '.fq'

    with open(insfq, 'w') as fq:
        for ins in insertions:
            for be in ('be1', 'be2'):
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

    fq = makefq(insertions)

    mapped = map(fq, ref, threads=int(args.threads))

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
    parser.add_argument('--invert', default=False, action='store_true', help='retain insertions that do not match library')
    args = parser.parse_args()
    main(args)
