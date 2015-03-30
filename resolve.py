#!/usr/bin/env python

import os
import shutil
import cPickle as pickle
import argparse
import logging
import tebreak
import subprocess

from uuid import uuid4

logger = logging.getLogger(__name__)


def prepare_ref(fasta, refoutdir='tebreak_refs'):
    if not os.path.exists(refoutdir):
        os.mkdir(refoutdir)
        assert os.path.exists(refoutdir), 'could not create ref output directory: %s' % refoutdir

    ref_fa = refoutdir + '/' + fasta

    if not os.path.exists(ref_fa):
        logger.debug('Copying %s to %s ...' % (fasta, ref_fa))
        shutil.copy(fasta, ref_fa)

    if not os.path.exists(ref_fa + '.fai'):
        logger.debug('Samtools indexing %s ...' % ref_fa)
        subprocess.call(['samtools', 'faidx', ref_fa])
        assert os.path.exists(ref_fa + '.fai'), 'could not samtools faidx %s' % ref_fa

    if not os.path.exists(ref_fa + '.tis'):
        logger.debug('Create LAST db for %s ...' % ref_fa)
        subprocess.call(['lastdb', '-s', '4G', ref_fa, ref_fa])
        assert os.path.exists(ref_fa + '.tis'), 'could not lastdb -4G %s %s' % (ref_fa, ref_fa)

    return ref_fa


def lastal_cons(ins, ref_fa, tmpdir='/tmp'):
    tmpfa = tmpdir + '/' + 'tebreak.resolve.%s.fa' % str(uuid4())
    with open(tmpfa, 'w') as fa:
        fa.write('>%s\n%s\n' % (ins['SR']['be1_obj_uuid'], ins['SR']['be1_dist_seq']))
        if 'be2_cons_seq' in ins['SR']:
            fa.write('>%s\n%s\n' % (ins['SR']['be2_obj_uuid'], ins['SR']['be2_dist_seq']))

    cmd = ['lastal', '-e 10', ref_fa, tmpfa]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    for line in p.stdout:
        if not line.startswith('#'):
            print line.strip()

    os.remove(tmpfa)



def main(args):
    if args.verbose: logger.setLevel(logging.DEBUG)
    insertions = []

    with open(args.pickle, 'r') as pickin:
        insertions = pickle.load(pickin)

    tebreak.text_summary(insertions) # debug

    ref_fa = prepare_ref(args.fasta, refoutdir=args.refoutdir)

    for ins in insertions:
        lastal_cons(ins, ref_fa)


if __name__ == '__main__':
    # set up logger
    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(format=FORMAT)
 
    parser = argparse.ArgumentParser(description='Resolve TE insertions from TEbreak data')
    parser.add_argument('-p', '--pickle', required=True)
    parser.add_argument('-f', '--fasta', required=True)
    parser.add_argument('--refoutdir', default='tebreak_refs')

    parser.add_argument('--tmpdir', default='/tmp')
    parser.add_argument('-v', '--verbose', action='store_true', default=False)
    args = parser.parse_args()
    main(args)