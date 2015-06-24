#!/usr/bin/env python

import cPickle as pickle
import argparse
import logging
import subprocess

from uuid import uuid4

logger = logging.getLogger(__name__)


def makefq(insertions):
    ''' write FASTA file from consensus sequences in insertions '''
    insfa = tmpdir + '/' + str(uuid4()) + '.fa'

    with open(insfa, 'w') as fa:
        for ins in insertions:
            for be in ('be1', 'be2'):
                if be+'_cons_seq' in ins['INFO'] and ins['INFO'][be+'_cons_seq']:
                    fa.write('>%s\n%s\n' % (ins['INFO']['ins_uuid'], ins['INFO'][be+'_cons_seq']))

    return insfa


def main(args):
    logger.debug('loading pickle: %s' % args.pickle)

    with open(args.pickle, 'r') as pickin:
        insertions = pickle.load(pickin)

    logger.debug('finished loading %s' % args.pickle)
    logger.debug('raw candidate count: %d' % len(insertions))

    print makefa(insertions)


if __name__ == '__main__':
    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(format=FORMAT)
    logger.setLevel(logging.DEBUG)

    parser = argparse.ArgumentParser(description='trim pickle based on alignments')
    parser.add_argument('-p', '--pickle', required=True)
    parser.add_argument('-r', '--ref', required=True)
    args = parser.parse_args()
    main(args)
