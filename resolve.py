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

    def pct_match(self):
        return float(sum([a==b for a,b in zip(list(self.query_align), list(self.target_align))])) / float(self.query_alnsize)

    def __lt__(self, other):
        return self.score > other.score

    def __gt__(self, other):
        return self.score < other.score

    def __str__(self):
        return "\n".join(self.raw)


#######################################
## Functions                         ##
#######################################


def best_match(last_results, query_name, req_target=None):
    qres = [res for res in sorted(last_results) if res.query_id == query_name and (req_target is None or req_target == res.target_id)]
    if len(qres) > 0:
        return qres[0]

    return None


def load_inslib(infa):
    seqdict = {}

    with open(infa, 'r') as fa:
        seqid = ''
        seq   = ''
        for line in fa:
            if line.startswith('>'):
                if seq != '':
                    seqdict[seqid] = seq
                seqid = line.lstrip('>').strip()
                seq   = ''
            else:
                assert seqid != ''
                seq = seq + line.strip()

    if seqid not in seqdict and seq != '':
        seqdict[seqid] = seq

    return seqdict


def prepare_ref(fasta, refoutdir='tebreak_refs', makeFAI=True, makeBWA=True, makeLAST=True):
    if not os.path.exists(refoutdir):
        os.mkdir(refoutdir)
        assert os.path.exists(refoutdir), 'could not create ref output directory: %s' % refoutdir

    ref_fa = refoutdir + '/' + os.path.basename(fasta)

    if not os.path.exists(ref_fa):
        logger.debug('Copying %s to %s ...' % (fasta, ref_fa))
        shutil.copy(fasta, ref_fa)

    if not os.path.exists(ref_fa + '.fai') and makeFAI:
        logger.debug('Samtools indexing %s ...' % ref_fa)
        subprocess.call(['samtools', 'faidx', ref_fa])
        assert os.path.exists(ref_fa + '.fai'), 'could not samtools faidx %s' % ref_fa

    if not os.path.exists(ref_fa + '.bwt') and makeBWA:
        logger.debug('Create BWA db for %s ...' % ref_fa)
        subprocess.call(['bwa', 'index', ref_fa])
        assert os.path.exists(ref_fa + '.bwt'), 'could not bwa index %s' % ref_fa

    if not os.path.exists(ref_fa + '.tis') and makeLAST:
        logger.debug('Create LAST db for %s ...' % ref_fa)
        subprocess.call(['lastdb', '-s', '4G', ref_fa, ref_fa])
        assert os.path.exists(ref_fa + '.tis'), 'could not lastdb -4G %s %s' % (ref_fa, ref_fa)

    return ref_fa


def lastal_cons(ins, ref_fa, tmpdir='/tmp'):
    tmpfa = tmpdir + '/' + 'tebreak.resolve.%s.fa' % str(uuid4())
    with open(tmpfa, 'w') as fa:
        fa.write('>%s\n%s\n' % (ins['SR']['be1_obj_uuid'], ins['SR']['be1_dist_seq'].split(',')[0]))
        if 'be2_dist_seq' in ins['SR'] and ins['SR']['be2_dist_seq'] is not None:
            fa.write('>%s\n%s\n' % (ins['SR']['be2_obj_uuid'], ins['SR']['be2_dist_seq'].split(',')[0]))

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
    be1_bestmatch = best_match(last_res, ins['SR']['be1_obj_uuid'])
    be2_bestmatch = None
    if 'be2_obj_uuid' in ins['SR'] and ins['SR']['be2_obj_uuid'] != ins['SR']['be1_obj_uuid']:
        be2_bestmatch = best_match(last_res, ins['SR']['be2_obj_uuid'])

    if None not in (be1_bestmatch, be2_bestmatch):
        if be1_bestmatch.target_id != be2_bestmatch.target_id:
            if be1_bestmatch.score > be2_bestmatch.score:
                # try to get be2 target to match be1 target
                newmatch = best_match(last_res, ins['SR']['be2_obj_uuid'], req_target=be1_bestmatch.target_id)
                if newmatch is not None and be2_bestmatch.score - newmatch.score < 5:
                    be2_bestmatch = newmatch

            else:
                # try to get be1 target to match be2 target
                newmatch = best_match(last_res, ins['SR']['be1_obj_uuid'], req_target=be2_bestmatch.target_id)
                if newmatch is not None and be1_bestmatch.score - newmatch.score < 5:
                    be1_bestmatch = newmatch

    if be1_bestmatch is not None: ins['SR']['be1_bestmatch'] = be1_bestmatch
    if be2_bestmatch is not None: ins['SR']['be2_bestmatch'] = be2_bestmatch

    ins = assign_insertion_ends(ins)

    ins['SR']['ins_length'] = infer_length(ins)

    ins = infer_orientation(ins)

    if be1_bestmatch is not None: ins['SR']['best_ins_matchpct'] = be1_bestmatch.pct_match()

    if be2_bestmatch is not None:
        if 'best_ins_matchpct' in ins['SR'] and be2_bestmatch.pct_match() > be1_bestmatch.pct_match():
            ins['SR']['best_ins_matchpct'] = be2_bestmatch.pct_match()
        elif 'best_ins_matchpct' not in ins['SR']:
            ins['SR']['best_ins_matchpct'] = be2_bestmatch.pct_match()

    return ins


def assign_insertion_ends(ins):
    be1 = None
    be2 = None

    if 'be1_bestmatch' in ins['SR']: be1 = ins['SR']['be1_bestmatch']
    if 'be2_bestmatch' in ins['SR']: be2 = ins['SR']['be2_bestmatch']

    if None not in (be1, be2): # compare to decide 3' end
        ins['SR']['be1_is_3prime'] = be1.target_start > be2.target_start
        ins['SR']['be2_is_3prime'] = not ins['SR']['be1_is_3prime']

    else: # for single end evidence, use coordinates to decide 3' end
        if be1 is not None:
            # 10bp to end of elt reference or has >= 10bp polyA we'll call it 3 prime
            ins['SR']['be1_is_3prime'] = be1.target_seqsize - (be1.target_start+be1.target_alnsize) < 10
            if not ins['SR']['be1_is_3prime']: ins['SR']['be1_is_3prime'] = be1.query_align.endswith('A'*10)
            if not ins['SR']['be1_is_3prime']: ins['SR']['be1_is_3prime'] = be1.target_align.endswith('A'*10)
            ins['SR']['be2_is_3prime'] = not ins['SR']['be1_is_3prime']

        elif be2 is not None:
            ins['SR']['be2_is_3prime'] = be2.target_seqsize - (be2.target_start+be2.target_alnsize) < 10
            if not ins['SR']['be2_is_3prime']: ins['SR']['be2_is_3prime'] = be2.query_align.endswith('A'*10)
            if not ins['SR']['be2_is_3prime']: ins['SR']['be2_is_3prime'] = be2.target_align.endswith('A'*10)
            ins['SR']['be1_is_3prime'] = not ins['SR']['be2_is_3prime']

    return ins


def infer_orientation(ins):
    # not sure if the prox. alignment must always be +, but it seems that way... testing assumption
    if 'be1_prox_str' in ins['SR']:
        assert ins['SR']['be1_prox_str'] == '+', "be1 negative strand alignment in proximal sequence"

    if 'be2_prox_str' in ins['SR']:
        assert ins['SR']['be2_prox_str'] == '+', "be2 negative strand alignment in proximal sequence"

    # defaults
    ins['SR']['be1_orient'] = None
    ins['SR']['be2_orient'] = None
    
    for be in ('be1', 'be2'):
        if be+'_prox_loc' in ins['SR'] and len(ins['SR'][be+'_prox_loc']) == 1:
            # work out which end of the consensus belongs to the distal sequence: _dist to proximal seq.
            right_dist = len(ins['SR'][be+'_prox_seq']) - max(ins['SR'][be+'_prox_loc'][0])
            left_dist  = min(ins['SR'][be+'_prox_loc'][0])

            distal_is_left = False

            if left_dist > right_dist:
                distal_is_left = True

            if be+'_is_3prime' in ins['SR']:
                if ins['SR'][be+'_is_3prime']:
                    if distal_is_left:
                        ins['SR'][be+'_orient'] = '+'
                    else:
                        ins['SR'][be+'_orient'] = '-'

                else:
                    if distal_is_left:
                        ins['SR'][be+'_orient'] = '-'
                    else:
                        ins['SR'][be+'_orient'] = '+'

    ins['SR']['inversion'] = ins['SR']['be1_orient'] != ins['SR']['be2_orient']

    return ins


def infer_length(ins):
    be1 = None
    be2 = None

    if 'be1_bestmatch' in ins['SR']: be1 = ins['SR']['be1_bestmatch']
    if 'be2_bestmatch' in ins['SR']: be2 = ins['SR']['be2_bestmatch']

    # length is known if both ends are known
    if None not in (be1, be2):
        coords = (be1.target_start+be1.target_alnsize, be1.target_start, be2.target_start+be2.target_alnsize, be2.target_start)
        return max(coords) - min(coords)

    # can reasonably approximate length if only 5' end is known
    if be1 is not None and not ins['SR']['be1_is_3prime']:
        return be1.target_seqsize - min(be1.target_start+be1.target_alnsize, be1.target_start)

    if be2 is not None and not ins['SR']['be2_is_3prime']:
        return be2.target_seqsize - min(be2.target_start+be2.target_alnsize, be2.target_start)

    # no idea
    return None


def best_ref(ins):
    be1 = None
    be2 = None

    if 'be1_bestmatch' in ins['SR']: be1 = ins['SR']['be1_bestmatch']
    if 'be2_bestmatch' in ins['SR']: be2 = ins['SR']['be2_bestmatch']

    if None not in (be1, be2):
        if be1.target_id == be2.target_id:
            return be1.target_id
        elif be1.score > be2.score:
            return be1.target_id
        elif be2.score < be1.score:
            return be2.target_id

    else:
        if be1 is not None: return be1.target_id
        if be2 is not None: return be2.target_id

    return None


def make_tmp_ref(ins, ref_fa, tmpdir='/tmp'):
    ref_id = best_ref(ins)
    if ref_id is None:
        return None

    inslib = load_inslib(ref_fa)
    assert ref_id in inslib, 'Reference is missing: %s' % ref_id

    tmp_ref = tmpdir + '/tebreak.ref.%s.%s.fa' % (ref_id, str(uuid4()))

    with open(tmp_ref, 'w') as fa:
        fa.write('>%s\n%s\n' % (ref_id, inslib[ref_id]))

    return prepare_ref(tmp_ref, refoutdir=tmpdir, makeLAST=False)


def remap_discordant(ins, ref_fa, tmpdir='/tmp'):
    if len(ins['READSTORE']) == 0:
        return None

    # make references for best target
    tmp_ref = make_tmp_ref(ins, ref_fa, tmpdir)

    if tmp_ref is None:
        return None

    tmp_fq  = '%s/tebreak.%s.discoremap.fq' % (tmpdir, str(uuid4()))
    tmp_sam = '.'.join(tmp_fq.split('.')[:-1]) + '.sam'
    tmp_bam = '.'.join(tmp_fq.split('.')[:-1]) + '.bam'
    tmp_srt = '.'.join(tmp_fq.split('.')[:-1]) + '.srt'

    with open(tmp_fq, 'w') as fq:
        for dr in ins['READSTORE']:
            fq.write(dr)

    sam_cmd = ['bwa', 'mem', '-k', '10', '-M', '-S', '-P', tmp_ref, tmp_fq]
    bam_cmd = ['samtools', 'view', '-bt', ref_fa + '.fai', '-o', tmp_bam, tmp_sam]
    srt_cmd = ['samtools', 'sort', tmp_bam, tmp_srt]
    idx_cmd = ['samtools', 'index', tmp_bam]

    with open(tmp_sam, 'w') as sam:
        p = subprocess.Popen(sam_cmd, stdout=subprocess.PIPE)
        for line in p.stdout:
            if not line.startswith('@'):
                if bin(int(line.split('\t')[1]) & 4) != bin(4): # not unmapped
                    sam.write(line)
            else:
                sam.write(line)

    subprocess.call(bam_cmd)
    
    subprocess.call(srt_cmd)
    shutil.move(tmp_srt+'.bam', tmp_bam)

    subprocess.call(idx_cmd)

    for ext in ('','.amb','.ann','.bck','.bwt','.des','.fai','.pac','.prj','.sa','.sds','.ssp','.suf','.tis'):
        if os.path.exists(tmp_ref+ext): os.remove(tmp_ref+ext)

    os.remove(tmp_sam)
    os.remove(tmp_fq)

    return tmp_bam


def resolve_insertion(args, ins, ref_fa):
    ''' add data based on alignments of library to consensus '''
    last_res = lastal_cons(ins, ref_fa)
    ins = add_insdata(ins, last_res)

    if 'best_ins_matchpct' in ins['SR']:
        if ins['DR']['dr_count'] > 0 and ins['SR']['best_ins_matchpct'] > 0.90: # change to parameter
            ins['DR']['tmpBAM'] = remap_discordant(ins, ref_fa, tmpdir=args.tmpdir)

    return ins


def resolve_unknown(insertions):
    ''' don't bother with alignment to consensus, attempt to assemble unknown insertions '''
    pass


def main(args):
    if args.verbose: logger.setLevel(logging.DEBUG)
    insertions = []

    with open(args.pickle, 'r') as pickin:
        insertions = pickle.load(pickin)

    ref_fa = prepare_ref(args.fasta, refoutdir=args.refoutdir, makeFAI=False, makeBWA=False)

    # parallelise insertions[]

    results = []

    for ins in insertions:
        results.append(resolve_insertion(args, ins, ref_fa))

    tebreak.text_summary(results) # debug


if __name__ == '__main__':
    # set up logger
    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(format=FORMAT)
 
    parser = argparse.ArgumentParser(description='Resolve TE insertions from TEbreak data')
    parser.add_argument('-p', '--pickle', required=True)
    parser.add_argument('-f', '--fasta', required=True, help='reference for insertions (not genome)')
    parser.add_argument('--refoutdir', default='tebreak_refs')

    parser.add_argument('--tmpdir', default='/tmp')
    parser.add_argument('-v', '--verbose', action='store_true', default=False)
    args = parser.parse_args()
    main(args)