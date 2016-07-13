#!/usr/bin/env python

import os
import cPickle as pickle
import argparse
import logging

from collections import OrderedDict as od

logger = logging.getLogger(__name__)


def main(args):
    logger.info('loading pickle: %s' % args.pickle)

    uuids = od() 

    with open(args.uuids) as uuid_in:
        for line in uuid_in:
            uuid = line.strip().split()[0]
            uuids[uuid] = True

    insertions = []

    with open(args.pickle, 'r') as pickin:
        insertions = pickle.load(pickin)

    logger.info('finished loading %s' % args.pickle)
    logger.info('raw candidate count: %d' % len(insertions))

    filtered = []

    for ins in insertions:
        if ins['INFO']['ins_uuid'] in uuids:
            filtered.append(ins)

    logger.info('kept %d records' % len(filtered))

    with open(args.out, 'w') as pickout:
        pickle.dump(filtered, pickout)


if __name__ == '__main__':
    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(format=FORMAT)
    logger.setLevel(logging.INFO)

    parser = argparse.ArgumentParser(description='filter pickle based on UUIDs')
    parser.add_argument('-p', '--pickle', required=True, help='input filename (tebreak.py pickle)')
    parser.add_argument('-u', '--uuids', required=True, help='list of UUIDs to keep (if >1 column, use the first column)')
    args = parser.parse_args()
    main(args)
