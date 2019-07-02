#!/usr/bin/env python

import os
import pickle
import argparse
import logging

from uuid import uuid4
from collections import defaultdict as dd


logger = logging.getLogger(__name__)

def output_fastq(ins, pickle, uuid):
    out_sr_fn = '.'.join(pickle.strip().split('.')[:-1]) + '.' + uuid + '.SR.fastq'
    out_dr_fn = '.'.join(pickle.strip().split('.')[:-1]) + '.' + uuid + '.DR.fastq'

    sr_count = 0
    dr_count = 0

    out_sr = open(out_sr_fn, 'w')
    out_dr = open(out_dr_fn, 'w')

    for read in ins['READSTORE']:
        if read.find('.SR/') > 0:
            out_sr.write(read)
            sr_count += 1

        if read.find('.DR/') > 0:
            out_dr.write(read)
            dr_count += 1

    out_sr.close()
    out_dr.close()
    
    return out_sr_fn, out_dr_fn, sr_count, dr_count


def main(args):
    logger.debug('loading pickle: %s' % args.pickle)

    with open(args.pickle, 'rb') as pickin:
        insertions = pickle.load(pickin)

    logger.debug('finished loading %s' % args.pickle)
    logger.debug('raw candidate count: %d' % len(insertions))

    uuids = {}
    with open(args.uuids) as _:
        for line in _:
            if not line.startswith('UUID') and not line.startswith ('#'):
                uuids[line.strip().split()[0]] = True

    for ins in insertions:
        if ins['INFO']['ins_uuid'] in uuids:
            if len(ins['READSTORE']) == 0:
                logger.warning('no reads for insertion: %s' % ins['INFO']['ins_uuid'])
                continue

            sr_fq, dr_fq, sr_count, dr_count = output_fastq(ins, args.pickle, ins['INFO']['ins_uuid'])

            logger.info('wrote %d split reads to %s' % (sr_count, sr_fq))
            logger.info('wrote %d discordant reads to %s' % (dr_count, dr_fq))


if __name__ == '__main__':
    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(format=FORMAT)
    logger.setLevel(logging.DEBUG)

    parser = argparse.ArgumentParser(description='output reads supporting insertions')
    parser.add_argument('-p', '--pickle', required=True, help='input filename (tebreak.py pickle)')
    parser.add_argument('-u', '--uuids', required=True, help='list of UUIDS in a .txt file - can use a tebreak table')
    args = parser.parse_args()
    main(args)
