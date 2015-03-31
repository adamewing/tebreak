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


#######################################
## Classes                           ##
#######################################


class LASTResult:
    def __init__(self, res):
        self.raw = res
        self.score = int(res[0].split()[1].replace('score=', ''))

        self.target_id      = res[1].split()[1]
        self.target_start   = int(res[1].split()[2])
        self.target_alnsize = int(res[1].split()[3])
        self.target_strand  = res[1].split()[4]
        self.target_seqsize = int(res[1].split()[5])
        self.target_align   = res[1].split()[6]

        self.query_id      = res[2].split()[1]
        self.query_start   = int(res[2].split()[2])
        self.query_alnsize = int(res[2].split()[3])
        self.query_strand  = res[2].split()[4]
        self.query_seqsize = int(res[2].split()[5])
        self.query_align   = res[2].split()[6]

    def __lt__(self, other):
        return self.score > other.score

    def __gt__(self, other):
        return self.score < other.score

    def __str__(self):
        return "\n".join(self.raw)


#######################################
## Functions                         ##
#######################################


def best_match(last_results, query_name):
    qres = [res for res in sorted(last_results) if res.query_id == query_name]
    if len(qres) > 0:
        return qres[0]

    return None


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

    cmd = ['lastal', '-e 20', ref_fa, tmpfa]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)

    la_lines   = []
    la_results = []

    for line in p.stdout:
        if not line.startswith('#'):
            if line.strip() != '':
                la_lines.append(line.strip())

            else:
                la_results.append(LASTResult(la_lines))
                la_lines = []

    os.remove(tmpfa)

    return la_results


def add_insdata(ins, last_res):
    be1_bestmatch = best_match(last_res, ins['be1_obj_uuid'])
    be2_bestmatch = None
    if 'be2_obj_uuid' in ins and ins['be2_obj_uuid'] != ins['be1_obj_uuid']:
        be2_bestmatch = best_match(last_res, ins['be2_obj_uuid'])

    print be1_bestmatch

    if be1_bestmatch is not None: ins['SR']['be1_bestmatch'] = be1_bestmatch
    if be2_bestmatch is not None: ins['SR']['be2_bestmatch'] = be2_bestmatch

    return ins


def resolve_insertions(insertions, ref_fa):
    for ins in insertions:
        last_res = lastal_cons(ins, ref_fa)
        ins = add_insdata(ins, last_res)

    tebreak.text_summary(insertions) # debug


def main(args):
    if args.verbose: logger.setLevel(logging.DEBUG)
    insertions = []

    with open(args.pickle, 'r') as pickin:
        insertions = pickle.load(pickin)

    ref_fa = prepare_ref(args.fasta, refoutdir=args.refoutdir)

    resolve_insertions(insertions, ref_fa)


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