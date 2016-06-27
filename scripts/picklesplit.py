#!/usr/bin/env python

import sys
import cPickle as pickle
import logging

from collections import defaultdict as dd

logger = logging.getLogger(__name__)
FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger.setLevel(logging.DEBUG)


if len(sys.argv) == 2:
    logger.debug('loading pickle: %s' % sys.argv[1])

    with open(sys.argv[1], 'r') as pickin:
        insertions = pickle.load(pickin)

    logger.debug('raw candidate count: %d' % len(insertions))

    ins_by_chrom = dd(list) 

    for ins in insertions:
        ins_by_chrom[ins['INFO']['chrom']].append(ins)

    for chrom in ins_by_chrom:
        outfn = '%s.%s.pickle' % ('.'.join(sys.argv[1].split('.')[:-1]), chrom)
        with open(outfn, 'w') as out:
            pickle.dump(ins_by_chrom[chrom], out)


else:
    sys.exit('usage: %s <pickle>' % sys.argv[0])

