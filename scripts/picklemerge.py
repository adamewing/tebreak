#!/usr/bin/env python


import os
import sys
import cPickle as pickle
import logging

logger = logging.getLogger(__name__)

from collections import OrderedDict as od

def resolve_duplicates(insertions):
    ''' resolve instances where breakpoints occur > 1x in the insertion list '''
    ''' this can happen if intervals overlap, e.g. in  genome chunking '''

    insdict = od() # --> index in insertions

    for n, ins in enumerate(insertions):
        be1 = ins['INFO']['chrom'] + ':' + str(ins['INFO']['be1_breakpos'])
        be2 = ins['INFO']['chrom'] + ':' + str(ins['INFO']['be2_breakpos'])
        
        if be1 not in insdict:
            insdict[be1] = n 
            insdict[be2] = n
        
        else:
            if prefer_insertion(ins, insertions[insdict[be1]]):
                insdict[be1] = n
                insdict[be2] = n

    return [insertions[n] for n in list(set(insdict.values()))]


def prefer_insertion(ins1, ins2):
    ''' return true if ins1 has more evidence than ins2, false otherwise '''
    # prefer two-end support over one end
    if ins1['INFO']['be1_breakpos'] != ins1['INFO']['be2_breakpos'] and ins2['INFO']['be1_breakpos'] == ins2['INFO']['be2_breakpos']:
        return True

    # prefer higher split read count
    if ins1['INFO']['be1_sr_count'] + ins1['INFO']['be2_sr_count'] > ins2['INFO']['be1_sr_count'] + ins2['INFO']['be2_sr_count']:
        return True

    # prefer higher discordant read count
    if ins1['INFO']['dr_count'] > ins2['INFO']['dr_count']:
        return True

    return False


if __name__ == '__main__':
    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(format=FORMAT)
    logger.setLevel(logging.INFO)


if len(sys.argv) > 2:
    assert not os.path.exists(sys.argv[1]), "file exists: %s, refusing to overwrite." % sys.argv[1]

    insertions = []

    for pfn in sys.argv[2:]:
        with open(pfn, 'r') as pickin:
            insertions += pickle.load(pickin)

        logger.info('loaded %s' % pfn)

    logger.info('finished loading pickles, found %d records, merging results...' % len(insertions))

    insertions = resolve_duplicates(insertions)

    logger.info('%d records remain after merge.' % len(insertions))

    with open(sys.argv[1], 'w') as pickout:
        pickle.dump(insertions, pickout)

    logger.info('wrote %d records to %s' % (len(insertions), sys.argv[1]))

else:
    sys.exit("usage: %s <output.pickle> <input1.pickle> <input2.pickle> ... ")