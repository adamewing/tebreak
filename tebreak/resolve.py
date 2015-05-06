#!/usr/bin/env python

import os
import sys
import shutil
import cPickle as pickle
import argparse
import logging
import tebreak
import pysam
import subprocess
import traceback

import multiprocessing as mp

from uuid import uuid4
from collections import OrderedDict as od
from collections import defaultdict as dd
from bx.intervals.intersection import Intersecter, Interval # pip install bx-python

logger = logging.getLogger(__name__)


#######################################
## Classes                           ##
#######################################


class TEIns:
    def __init__(self, ins):
        self.ins = ins['INFO']
        self.out = od()
        self.end3 = 'NA'
        self.end5 = 'NA'

        self._fill_out()

    def _fill_out(self):
        self.out['Chromosome']    = self.ins['chrom']
        self.out['Left_Extreme']  = self.ins['min_supporting_base']
        self.out['Right_Extreme'] = self.ins['max_supporting_base']
        self.assign35ends()
        self.te_family()
        self.elt_coords()
        self.support()
        self.improve_cons()
        self.tsd()
        self.sample_list()
        self.consensus()

    def assign35ends(self):
        self.out['5_Prime_End'] = 'NA'
        self.out['3_Prime_End'] = 'NA'

        if 'be1_is_3prime' in self.ins:
            if self.ins['be1_is_3prime']:
                self.out['3_Prime_End'] = self.ins['be1_breakpos']
                self.out['5_Prime_End'] = self.ins['be2_breakpos']
                self.end3 = 'be1'
                self.end5 = 'be2'
            else:
                self.out['3_Prime_End'] = self.ins['be2_breakpos']
                self.out['5_Prime_End'] = self.ins['be1_breakpos']
                self.end3 = 'be2'
                self.end5 = 'be1'

    def te_family(self):
        ''' inslib input should have headers superfamily:subfamily e.g. L1:L1Ta or Alu:AluYa5 '''
        self.out['Superfamily'] = []
        self.out['Subfamily']   = []

        for be in ('be1', 'be2'):
            if be + '_bestmatch' in self.ins:
                superfam, subfam = self.ins[be + '_bestmatch'].target_id.split(':')
                self.out['Superfamily'].append(superfam)
                self.out['Subfamily'].append(subfam)

        self.out['Superfamily'] = ','.join(list(set(self.out['Superfamily'])))
        self.out['Subfamily']   = ','.join(list(set(self.out['Subfamily'])))

    def elt_coords(self):
        self.out['TE_Align_Start'] = 'NA'
        self.out['TE_Align_End']   = 'NA'

        self.out['Orient_5p'] = 'NA'
        self.out['Orient_3p'] = 'NA'

        if self.end5 + '_orient' in self.ins: self.out['Orient_5p'] = self.ins[self.end5 + '_orient']
        if self.end3 + '_orient' in self.ins: self.out['Orient_3p'] = self.ins[self.end3 + '_orient']

        if self.end5 + '_bestmatch' in self.ins: self.out['TE_Align_Start'] = self.ins[self.end5 + '_bestmatch'].target_start
        if self.end3 + '_bestmatch' in self.ins: self.out['TE_Align_End']   = self.ins[self.end3 + '_bestmatch'].target_end

        if 'remap_min_pos' in self.ins:
            if self.out['TE_Align_Start'] == 'NA' or self.out['TE_Align_Start'] > self.ins['remap_min_pos']:
                self.out['TE_Align_Start'] = self.ins['remap_min_pos']

        if 'remap_min_pos' in self.ins:
            
            if self.out['TE_Align_End'] == 'NA' or self.out['TE_Align_End'] < self.ins['remap_max_pos']:
                self.out['TE_Align_End'] = self.ins['remap_max_pos']

        self.out['Inversion'] = 'NA'
        if 'inversion' in self.ins:
            if self.ins['inversion']: self.out['Inversion'] = 'Y'
            else: self.out['Inversion'] = 'N'

    def junctions(self):
        j = []
        if 'be1_breakpos' in self.ins: j.append(self.ins['be1_breakpos'])
        if 'be2_breakpos' in self.ins: j.append(self.ins['be2_breakpos'])
        return j

    def support(self):
        self.out['5p_Elt_Match'] = 0.0
        self.out['3p_Elt_Match'] = 0.0

        if self.end5 + '_bestmatch' in self.ins: self.out['5p_Elt_Match'] = self.ins[self.end5 + '_bestmatch'].pct_match()
        if self.end3 + '_bestmatch' in self.ins: self.out['3p_Elt_Match'] = self.ins[self.end3 + '_bestmatch'].pct_match()

        self.out['5p_Genome_Match'] = 0.0
        self.out['3p_Genome_Match'] = 0.0

        if self.end5 + '_avgmatch' in self.ins: self.out['5p_Genome_Match'] = self.ins[self.end5 + '_avgmatch']
        if self.end3 + '_avgmatch' in self.ins: self.out['3p_Genome_Match'] = self.ins[self.end3 + '_avgmatch']

        self.out['Split_reads_3prime'] = 0
        self.out['Split_reads_5prime'] = 0

        if self.end3 + '_sr_count' in self.ins: self.out['Split_reads_3prime'] = self.ins[self.end3 + '_sr_count']
        if self.end5 + '_sr_count' in self.ins: self.out['Split_reads_5prime'] = self.ins[self.end5 + '_sr_count']

        self.out['Remapped_Discordant'] = 0
        self.out['Remapped_Splitreads'] = 0

        if 'remap_dr_count' in self.ins: self.out['Remapped_Discordant'] = self.ins['remap_dr_count']
        if 'remap_dr_count' in self.ins: self.out['Remapped_Splitreads'] = self.ins['remap_sr_count']

        self.out['5p_Cons_Len'] = 0
        self.out['3p_Cons_Len'] = 0

        if self.end5 + '_cons_seq' in self.ins: self.out['5p_Cons_Len'] = len(self.ins[self.end5 + '_cons_seq'])
        if self.end3 + '_cons_seq' in self.ins: self.out['3p_Cons_Len'] = len(self.ins[self.end3 + '_cons_seq'])

    def tsd(self):
        self.out['TSD_3prime'] = 'NA'
        self.out['TSD_5prime'] = 'NA'

        if self.end3 + '_end_over' in self.ins: self.out['TSD_3prime'] = self.ins[self.end3 + '_end_over']
        if self.end5 + '_end_over' in self.ins: self.out['TSD_5prime'] = self.ins[self.end5 + '_end_over']

    def improve_cons(self):
        self.out['5p_Improved'] = 'N'
        self.out['3p_Improved'] = 'N'

        if self.end5 + '_improved' in self.ins and self.ins[self.end5 + '_improved']:
            self.out['5p_Improved'] = 'Y'

        if self.end3 + '_improved' in self.ins and self.ins[self.end3 + '_improved']:
            self.out['3p_Improved'] = 'Y'

    def sample_list(self):
        samples = dd(int)

        for be in ('be1', 'be2'):
            if be + '_bf_count' in self.ins:
                for sample_support in self.ins[be + '_bf_count']:
                    sample, count = sample_support.split('|')
                    samples[sample] += int(count)

            self.out['Sample_count']   = len(samples) 
            self.out['Sample_support'] = ','.join(['%s|%d' % (sample, count) for sample, count in samples.iteritems()])

    def consensus(self):
        self.out['Consensus_5p'] = 'NA'
        self.out['Consensus_3p'] = 'NA'

        if self.end5 + '_cons_seq' in self.ins: self.out['Consensus_5p'] = self.ins[self.end5 + '_cons_seq']
        if self.end3 + '_cons_seq' in self.ins: self.out['Consensus_3p'] = self.ins[self.end3 + '_cons_seq']

    def polyA_filter(self):
        ''' return true if only supported by poly-A matches '''
        found_non_pA = False

        if 'be1_bestmatch' in self.ins:
            if poly_A_frac(self.ins['be1_bestmatch'].target_align) < 0.85: found_non_pA = True

        if 'be2_bestmatch' in self.ins:
            if poly_A_frac(self.ins['be2_bestmatch'].target_align) < 0.85: found_non_pA = True

        return not found_non_pA

    def genome_location_filter(self, forest):
        if self.ins['chrom'] in forest:
            hits = forest[self.ins['chrom']].find(self.out['Left_Junction'], self.out['Right_Junction']+1)
            for hit in hits:
                    return True

        return False

    def pass_filter(self, forest):
        passed = True
        if 'best_ins_matchpct' in self.ins:
            if self.ins['best_ins_matchpct'] < 0.9: passed = False
            if self.polyA_filter(): passed = False
        else: passed = False

        if len(self.ins['be1_prox_seq']) < 20: passed = False

        if forest is not None:
            if self.genome_location_filter(forest): passed = False

        if 'transducer' in self.ins: passed = False

        return passed

    def header(self):
        return '\t'.join(self.out.keys())

    def __lt__(self, other):
        if self.out['Chromosome'] == other.out['Chromosome']:
            return self.out['Left_Extreme'] < other.out['Left_Extreme']

        return self.out['Chromosome'] < other.out['Chromosome']

    def __str__(self):
        return '\t'.join(map(str, self.out.values()))


#######################################
## Functions                         ##
#######################################


def best_match(last_results, query_name, req_target=None, min_match=0.9):
    qres = []

    for res in sorted(last_results):
        if res.query_id == query_name:
            if req_target is None or req_target == res.target_id:
                if res.pct_match() > min_match:
                    return res

    return None


def prepare_ref(fasta, refoutdir='tebreak_refs', makeFAI=True, makeBWA=True, makeLAST=True):
    if not os.path.exists(refoutdir):
        os.mkdir(refoutdir)
        assert os.path.exists(refoutdir), 'could not create ref output directory: %s' % refoutdir

    ref_fa = refoutdir + '/' + os.path.basename(fasta)

    if not os.path.exists(ref_fa):
        logger.debug('Copying %s to %s ...' % (fasta, ref_fa))
        shutil.copy(fasta, ref_fa)

    if makeFAI:
        logger.debug('Samtools indexing %s ...' % ref_fa)
        subprocess.call(['samtools', 'faidx', ref_fa])
        assert os.path.exists(ref_fa + '.fai'), 'could not samtools faidx %s' % ref_fa

    if makeBWA:
        logger.debug('Create BWA db for %s ...' % ref_fa)
        p = subprocess.Popen(['bwa', 'index', ref_fa], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        for line in p.stdout: pass
        assert os.path.exists(ref_fa + '.bwt'), 'could not bwa index %s' % ref_fa

    if makeLAST:
        logger.debug('Create LAST db for %s ...' % ref_fa)
        subprocess.call(['lastdb', '-s', '4G', ref_fa, ref_fa])
        assert os.path.exists(ref_fa + '.tis'), 'could not lastdb -4G %s %s' % (ref_fa, ref_fa)

    return ref_fa


def last_alignment(ins, ref_fa, tmpdir='/tmp'):
    tmpfa = tmpdir + '/' + 'tebreak.resolve.%s.fa' % str(uuid4())
    with open(tmpfa, 'w') as fa:
        # realign all distal sequences
        if 'be1_dist_seq' in ins['INFO'] and ins['INFO']['be1_dist_seq'] is not None:
            for i, dist_seq in enumerate(ins['INFO']['be1_dist_seq'].split(',')):
                fa.write('>%s|%s\n%s\n' % (ins['INFO']['be1_obj_uuid'],str(i), dist_seq))

        if 'be2_dist_seq' in ins['INFO'] and ins['INFO']['be2_dist_seq'] is not None:
            for i, dist_seq in enumerate(ins['INFO']['be2_dist_seq'].split(',')):
                fa.write('>%s|%s\n%s\n' % (ins['INFO']['be2_obj_uuid'], str(i), dist_seq))

        # realign all unmapped sequences, use negative index (-1-i) as flag for unmapped
        if 'be1_umap_seq' in ins['INFO'] and ins['INFO']['be1_umap_seq'] is not None:
            for i, umap_seq in enumerate(ins['INFO']['be1_umap_seq'].split(',')):
                fa.write('>%s|%s\n%s\n' % (ins['INFO']['be1_obj_uuid'],str(-1-i), umap_seq))

        if 'be2_umap_seq' in ins['INFO'] and ins['INFO']['be2_umap_seq'] is not None:
            for i, umap_seq in enumerate(ins['INFO']['be2_umap_seq'].split(',')):
                fa.write('>%s|%s\n%s\n' % (ins['INFO']['be2_obj_uuid'], str(-1-i), umap_seq))

    # increasing -m accounts for poly-A tails
    cmd = ['lastal', '-e', '20', '-m', '100', ref_fa, tmpfa]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)

    maf_lines   = []
    maf_results = []

    for line in p.stdout:
        if not line.startswith('#'):
            if line.strip() != '':
                maf_lines.append(line.strip())

            else:
                maf_results.append(tebreak.LASTResult(maf_lines))
                maf_lines = []

    os.remove(tmpfa)

    return maf_results


def interval_forest(bed_file):
    ''' build dictionary of interval trees; value=LASTResult '''
    forest = dd(Intersecter)

    with open(bed_file, 'r') as bed:
        for line in bed:
            chrom, start, end = line.strip().split()[:3]
            info = line.strip().split()[3:]
            forest[chrom].add_interval(Interval(int(start), int(end), value=info))

    return forest


def poly_A_frac(seq):
    ''' home much of the sequence is A '''
    if seq is None: return 0.0

    f = []
    for b in ('A', 'T'):
        f.append(len([base for base in list(seq.upper()) if base == b]))

    return float(max(f)) / float(len(seq))


def bamcount(ins):
    bams = []
    for be in ('be1', 'be2'):
        if be + '_bf_count' in ins['INFO']:
            bams += [bam.split('|')[0] for bam in ins['INFO'][be + '_bf_count']]

    return len(list(set(bams)))


def add_insdata(ins, last_res, max_bam_count=0):
    be1_bestmatch = best_match(last_res, ins['INFO']['be1_obj_uuid'])
    be2_bestmatch = None
    if 'be2_obj_uuid' in ins['INFO'] and ins['INFO']['be2_obj_uuid'] != ins['INFO']['be1_obj_uuid']:
        be2_bestmatch = best_match(last_res, ins['INFO']['be2_obj_uuid'])

    better_match_maxdiff = 5

    if None not in (be1_bestmatch, be2_bestmatch):
        if be1_bestmatch.target_id != be2_bestmatch.target_id:
            if be1_bestmatch.score > be2_bestmatch.score:
                # try to get be2 target to match be1 target
                if max(poly_A_frac(ins['INFO']['be2_dist_seq']), poly_A_frac(ins['INFO']['be2_umap_seq'])) > 0.75:
                    better_match_maxdiff = 20 # if seq is poly_A make it easy to match

                newmatch = best_match(last_res, ins['INFO']['be2_obj_uuid'], req_target=be1_bestmatch.target_id)
                if newmatch is not None and be2_bestmatch.score - newmatch.score < better_match_maxdiff:
                    be2_bestmatch = newmatch

            else:
                # try to get be1 target to match be2 target
                if max(poly_A_frac(ins['INFO']['be1_dist_seq']), poly_A_frac(ins['INFO']['be1_umap_seq'])) > 0.75:
                    better_match_maxdiff = 20

                newmatch = best_match(last_res, ins['INFO']['be1_obj_uuid'], req_target=be2_bestmatch.target_id)
                if newmatch is not None and be1_bestmatch.score - newmatch.score < better_match_maxdiff:
                    be1_bestmatch = newmatch

    if be1_bestmatch is not None: ins['INFO']['be1_bestmatch'] = be1_bestmatch
    if be2_bestmatch is not None: ins['INFO']['be2_bestmatch'] = be2_bestmatch

    ins = assign_insertion_ends(ins)

    ins['INFO']['ins_length'] = infer_length(ins)

    ins = infer_orientation(ins)

    if be1_bestmatch is not None: ins['INFO']['best_ins_matchpct'] = be1_bestmatch.pct_match()

    if be2_bestmatch is not None:
        if 'best_ins_matchpct' in ins['INFO'] and be2_bestmatch.pct_match() > be1_bestmatch.pct_match():
            ins['INFO']['best_ins_matchpct'] = be2_bestmatch.pct_match()
        elif 'best_ins_matchpct' not in ins['INFO']:
            ins['INFO']['best_ins_matchpct'] = be2_bestmatch.pct_match()

    return ins


def assign_insertion_ends(ins):
    be1 = None
    be2 = None

    if 'be1_bestmatch' in ins['INFO']: be1 = ins['INFO']['be1_bestmatch']
    if 'be2_bestmatch' in ins['INFO']: be2 = ins['INFO']['be2_bestmatch']

    if None not in (be1, be2): # compare to decide 3' end
        ins['INFO']['be1_is_3prime'] = be1.target_start > be2.target_start
        ins['INFO']['be2_is_3prime'] = not ins['INFO']['be1_is_3prime']

    else: # for single end evidence, use coordinates to decide 3' end
        if be1 is not None:
            # 10bp to end of elt reference or has >= 10bp polyA we'll call it 3 prime
            ins['INFO']['be1_is_3prime'] = be1.target_seqsize - (be1.target_start+be1.target_alnsize) < 10
            if not ins['INFO']['be1_is_3prime']: ins['INFO']['be1_is_3prime'] = be1.query_align.endswith('A'*10)
            if not ins['INFO']['be1_is_3prime']: ins['INFO']['be1_is_3prime'] = be1.target_align.endswith('A'*10)
            ins['INFO']['be2_is_3prime'] = not ins['INFO']['be1_is_3prime']

        elif be2 is not None:
            ins['INFO']['be2_is_3prime'] = be2.target_seqsize - (be2.target_start+be2.target_alnsize) < 10
            if not ins['INFO']['be2_is_3prime']: ins['INFO']['be2_is_3prime'] = be2.query_align.endswith('A'*10)
            if not ins['INFO']['be2_is_3prime']: ins['INFO']['be2_is_3prime'] = be2.target_align.endswith('A'*10)
            ins['INFO']['be1_is_3prime'] = not ins['INFO']['be2_is_3prime']

    return ins


def infer_orientation(ins):

    # defaults
    ins['INFO']['be1_orient'] = None
    ins['INFO']['be2_orient'] = None
    
    # in case there is more than one distal mapping (e.g. transduced seq.)
    be1_distal_num = 0
    be2_distal_num = 0

    if 'be1_bestmatch' in ins['INFO'] and ins['INFO']['be1_bestmatch'] is not None:
        be1_distal_num = ins['INFO']['be1_bestmatch'].query_distnum

    if 'be2_bestmatch' in ins['INFO'] and ins['INFO']['be2_bestmatch'] is not None:
        be1_distal_num = ins['INFO']['be2_bestmatch'].query_distnum

    for be in ('be1', 'be2'):
        if be+'_prox_loc' in ins['INFO'] and len(ins['INFO'][be+'_prox_loc']) > 0:
            # work out which end of the consensus belongs to the distal sequence: _dist to proximal seq.
            # more than one proximal mapping may indicate the insertion is completely assembled (but not always)
            which_prox = 0
            if be+'_use_prox' in ins['INFO']: which_prox = ins['INFO'][be+'_use_prox']

            right_dist = len(ins['INFO'][be+'_prox_seq']) - max(ins['INFO'][be+'_prox_loc'][which_prox])
            left_dist  = min(ins['INFO'][be+'_prox_loc'][which_prox])

            distal_is_left = False

            if left_dist > right_dist:
                distal_is_left = True

            if be+'_is_3prime' in ins['INFO']:
                if ins['INFO'][be+'_is_3prime']:
                    if distal_is_left:
                        ins['INFO'][be+'_orient'] = '+'
                    else:
                        ins['INFO'][be+'_orient'] = '-'

                else:
                    if distal_is_left:
                        ins['INFO'][be+'_orient'] = '-'
                    else:
                        ins['INFO'][be+'_orient'] = '+'

                if ins['INFO'][be+'_prox_str'] == '-': ins['INFO'][be+'_orient'] = swapstrand(ins['INFO'][be+'_orient'])

    if None not in (ins['INFO']['be1_orient'], ins['INFO']['be2_orient']):
        ins['INFO']['inversion'] = ins['INFO']['be1_orient'] != ins['INFO']['be2_orient']

    return ins


def swapstrand(strand):
    if strand == '+': return '-'
    if strand == '-': return '+'


def infer_length(ins):
    be1 = None
    be2 = None

    if 'be1_bestmatch' in ins['INFO']: be1 = ins['INFO']['be1_bestmatch']
    if 'be2_bestmatch' in ins['INFO']: be2 = ins['INFO']['be2_bestmatch']

    # length is known if both ends are known
    if None not in (be1, be2):
        coords = (be1.target_start+be1.target_alnsize, be1.target_start, be2.target_start+be2.target_alnsize, be2.target_start)
        return max(coords) - min(coords)

    # can reasonably approximate length if only 5' end is known
    if be1 is not None and not ins['INFO']['be1_is_3prime']:
        return be1.target_seqsize - min(be1.target_start+be1.target_alnsize, be1.target_start)

    if be2 is not None and not ins['INFO']['be2_is_3prime']:
        return be2.target_seqsize - min(be2.target_start+be2.target_alnsize, be2.target_start)

    # no idea
    return None


def best_ref(ins):
    be1 = None
    be2 = None

    if 'be1_bestmatch' in ins['INFO']: be1 = ins['INFO']['be1_bestmatch']
    if 'be2_bestmatch' in ins['INFO']: be2 = ins['INFO']['be2_bestmatch']

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


def score_insertion(ins):
    score = 0
    bothends = 'be1_cons_seq' in ins['INFO'] and 'be2_cons_seq' in ins['INFO']

    if 'be1_sr_count'  in ins['INFO']: score += ins['INFO']['be1_sr_count']
    if 'be2_sr_count'  in ins['INFO']: score += ins['INFO']['be2_sr_count']
    if 'mapped_target' in ins['INFO']: score += ins['INFO']['mapped_target']

    if bothends: score *= 2

    return score


def make_tmp_ref(ins, ref_fa, tmpdir='/tmp'):
    ref_id = best_ref(ins)
    if ref_id is None:
        return None

    inslib = tebreak.load_falib(ref_fa)
    assert ref_id in inslib, 'Reference is missing: %s' % ref_id

    tmp_ref = tmpdir + '/tebreak.ref.%s.%s.fa' % (ref_id, str(uuid4()))

    with open(tmp_ref, 'w') as fa:
        fa.write('>%s\n%s\n' % (ref_id, inslib[ref_id]))

    return prepare_ref(tmp_ref, refoutdir=tmpdir, makeLAST=False)


def remap_discordant(ins, inslib_fa=None, useref=None, tmpdir='/tmp'):
    ''' will build temporary ref from inslib_fasta unless useref is specified '''
    if len(ins['READSTORE']) == 0:
        return None

    if inslib_fa is not None:
    # make references for best target
        tmp_ref = make_tmp_ref(ins, inslib_fa, tmpdir)
        if tmp_ref is None: return None

    if useref is not None: tmp_ref = useref

    if tmp_ref is None: return None

    tmp_fq  = '%s/tebreak.%s.discoremap.fq' % (tmpdir, str(uuid4()))
    tmp_sam = '.'.join(tmp_fq.split('.')[:-1]) + '.sam'
    tmp_bam = '.'.join(tmp_fq.split('.')[:-1]) + '.bam'
    tmp_srt = '.'.join(tmp_fq.split('.')[:-1]) + '.srt'

    with open(tmp_fq, 'w') as fq:
        for dr in ins['READSTORE']:
            fq.write(dr)

    sam_cmd = ['bwa', 'mem', '-v', '1', '-k', '10', '-M', '-S', '-P', tmp_ref, tmp_fq]
    bam_cmd = ['samtools', 'view', '-bt', tmp_ref + '.fai', '-o', tmp_bam, tmp_sam]
    srt_cmd = ['samtools', 'sort', tmp_bam, tmp_srt]
    idx_cmd = ['samtools', 'index', tmp_bam]

    FNULL = open(os.devnull, 'w')

    with open(tmp_sam, 'w') as sam:
        p = subprocess.Popen(sam_cmd, stdout=subprocess.PIPE, stderr=FNULL)
        for line in p.stdout:
            if not line.startswith('@'):
                if bin(int(line.split('\t')[1]) & 4) != bin(4): # not unmapped
                    sam.write(line)
            else:
                sam.write(line)

    p = subprocess.Popen(bam_cmd, stdout=subprocess.PIPE, stderr=FNULL)
    for line in p.stdout: pass
    
    subprocess.call(srt_cmd)
    shutil.move(tmp_srt+'.bam', tmp_bam)

    subprocess.call(idx_cmd)

    for ext in ('','.amb','.ann','.bck','.bwt','.des','.fai','.pac','.prj','.sa','.sds','.ssp','.suf','.tis'):
        if os.path.exists(tmp_ref+ext): os.remove(tmp_ref+ext)

    os.remove(tmp_sam)
    os.remove(tmp_fq)

    return tmp_bam


def identify_transductions(ins):
    bedict = {'be1': None, 'be2': None}

    if 'be1_bestmatch' in ins['INFO']: bedict['be1'] = ins['INFO']['be1_bestmatch']
    if 'be2_bestmatch' in ins['INFO']: bedict['be2'] = ins['INFO']['be2_bestmatch']

    for be in ('be1','be2'):
        if bedict[be] is not None:
            tr_seqs = []
            tr_locs = [] # chrom, start, end, 3p or 5p / unmap

            segtype = 'dist'
            if bedict[be].query_distnum < 0: segtype = 'umap'

            num_segs = 0
            if ins['INFO'][be+'_'+segtype+'_seq'] is not None: num_segs = len(ins['INFO'][be+'_'+segtype+'_seq'].split(','))

            # more than one sequence
            if num_segs > 1:
                for segnum in range(num_segs):
                    if segtype == 'umap': segnum = -1-segnum

                    if segnum != bedict[be].query_distnum:
                        if segtype == 'dist':
                            tr_chrom = ins['INFO'][be+'_'+segtype+'_chr'].split(',')[segnum]
                            tr_start = ins['INFO'][be+'_'+segtype+'_pos'].split(',')[segnum]
                            tr_end   = ins['INFO'][be+'_'+segtype+'_end'].split(',')[segnum]

                        tr_side  = '5p'
                        if ins['INFO'][be+'_is_3prime']: tr_side = '3p'

                        tr_seqs.append(ins['INFO'][be+'_'+segtype+'_seq'].split(',')[segnum])
                        if segtype == 'dist': tr_locs.append((tr_chrom, tr_start, tr_end, tr_side))
                        if segtype == 'umap': tr_locs.append(('unmap', 0, 0, tr_side))

            # mapped distal sequence, unmapped transduced seq
            if segtype == 'dist':
                if num_segs > 0 and ins['INFO'][be+'_umap_seq'] is not None: # unmapped transduced seqs
                    for unmap_seq in ins['INFO'][be+'_umap_seq'].split(','):
                        tr_side  = '5p'
                        if ins['INFO'][be+'_is_3prime']: tr_side = '3p'

                        tr_seqs.append(unmap_seq)
                        tr_locs.append(('unmap', 0, 0, tr_side))
            else:
                if num_segs > 0 and ins['INFO'][be+'_dist_seq'] is not None:
                    for i, mapped_seq in enumerate(ins['INFO'][be+'_dist_seq'].split(',')):
                        tr_side  = '5p'
                        if ins['INFO'][be+'_is_3prime']: tr_side = '3p'

                        tr_seqs.append(mapped_seq)
                        tr_locs.append((ins['INFO'][be+'_dist_chr'][i], 
                                        ins['INFO'][be+'_dist_pos'].split(',')[i], 
                                        ins['INFO'][be+'_dist_end'].split(',')[i], 
                                        tr_side
                                       )
                        )

            if len(tr_seqs) > 0:
                ins['INFO'][be+'_trans_seq'] = tr_seqs
                ins['INFO'][be+'_trans_loc'] = tr_locs

    return ins


def get_bam_info(bam, ins):
    max_positions = []
    min_positions = []
    sr_count = 0
    dr_count = 0

    for read in bam.fetch():
        max_positions.append(read.get_reference_positions()[-1])
        min_positions.append(read.get_reference_positions()[0])
        if read.qname.split('.')[-1] == 'DR': dr_count += 1
        if read.qname.split('.')[-1] == 'SR': sr_count += 1

    if len(min_positions) > 0: ins['INFO']['remap_min_pos'] = min(min_positions)
    if len(max_positions) > 0: ins['INFO']['remap_max_pos'] = max(max_positions)

    ins['INFO']['remap_sr_count'] = sr_count
    ins['INFO']['remap_dr_count'] = dr_count

    return ins


def resolve_insertion(args, ins, inslib_fa):
    ''' add data based on alignments of library to consensus '''
    try:
        last_res = last_alignment(ins, inslib_fa)
        ins = add_insdata(ins, last_res)

        if not args.skip_align and ins['INFO']['dr_count'] > 0:
            tmp_bam = remap_discordant(ins, inslib_fa=inslib_fa, tmpdir=args.tmpdir)

            if tmp_bam is not None:
                bam = pysam.AlignmentFile(tmp_bam, 'rb')
                ins['INFO']['support_bam_file'] = tmp_bam
                ins['INFO']['mapped_target'] = bam.mapped
                ins = get_bam_info(bam, ins)

                if not args.keep_tmp_bams:
                    if os.path.exists(tmp_bam): os.remove(tmp_bam)
                    if os.path.exists(tmp_bam + '.bai'): os.remove(tmp_bam + '.bai')

        if 'best_ins_matchpct' in ins['INFO']: ins = identify_transductions(ins)

        return ins

    except Exception, e:
        sys.stderr.write('*'*60 + '\tencountered error:\n')
        traceback.print_exc(file=sys.stderr)

        if ins and 'INFO' in ins and 'chrom' in ins['INFO'] and 'be1_breakpos' in ins['INFO']:
            sys.stderr.write("Insertion location: %s:%d\n" % (ins['INFO']['chrom'], ins['INFO']['be1_breakpos']))

        sys.stderr.write("*"*60 + "\n")

        return None


def resolve_transductions(insertions):
    ''' some insertions are actually transductions; try to work this out '''
    tr_coord_dict = dd(dict) # chrom --> pos --> [ins]

    for ins in insertions:
            for be in ('be1','be2'):
                if be+'_trans_loc' in ins['INFO']:
                    for (chrom, start, end, side) in ins['INFO'][be+'_trans_loc']:
                        for pos in map(int, (start, end)):
                            if pos not in tr_coord_dict[chrom]: tr_coord_dict[chrom][pos] = []
                            tr_coord_dict[chrom][pos].append(ins)

    for parent_ins in insertions:
        for be in ('be1','be2'):
            if parent_ins['INFO']['chrom'] in tr_coord_dict and int(parent_ins['INFO'][be+'_breakpos']) in tr_coord_dict[parent_ins['INFO']['chrom']]:
                for trduct_ins in tr_coord_dict[parent_ins['INFO']['chrom']][int(parent_ins['INFO'][be+'_breakpos'])]:

                    trduct_loc_string = '%s:%d-%d' % (trduct_ins['INFO']['chrom'], trduct_ins['INFO']['be1_breakpos'], trduct_ins['INFO']['be2_breakpos'])
                    parent_loc_string = '%s:%d-%d' % (parent_ins['INFO']['chrom'], parent_ins['INFO']['be1_breakpos'], parent_ins['INFO']['be2_breakpos'])

                    if score_insertion(parent_ins) < score_insertion(trduct_ins):
                        if 'transducer' not in parent_ins['INFO']: parent_ins['INFO']['transducer'] = []
                        if 'transduction' not in trduct_ins['INFO']: trduct_ins['INFO']['transduction'] = []

                        parent_ins['INFO']['transducer'].append(trduct_loc_string)
                        trduct_ins['INFO']['transduction'].append(parent_loc_string)

                    else:
                        if 'transducer' not in trduct_ins['INFO']: trduct_ins['INFO']['transducer'] = []
                        if 'transduction' not in parent_ins['INFO']: parent_ins['INFO']['transduction'] = []

                        parent_ins['INFO']['transduction'].append(trduct_loc_string)
                        trduct_ins['INFO']['transducer'].append(parent_loc_string)

    return insertions


def main(args):
    if args.verbose: logger.setLevel(logging.DEBUG)
    insertions = []

    with open(args.pickle, 'r') as pickin:
        insertions = pickle.load(pickin)

    inslib_fa = prepare_ref(args.inslib_fasta, refoutdir=args.refoutdir, makeFAI=False, makeBWA=False)

    forest = None
    # set up reference mask BED if available
    if args.filter_bed is not None:
        logger.debug('using BED as genome reference filter: %s' % args.filter_bed)
        forest = interval_forest(args.filter_bed)

    results = []

    pool = mp.Pool(processes=int(args.processes))

    for ins in insertions:
        if int(args.max_bam_count) == 0 or bamcount(ins) <= int(args.max_bam_count):
            res = pool.apply_async(resolve_insertion, [args, ins, inslib_fa])
            results.append(res)

    insertions = [res.get() for res in results if res is not None]

    insertions = resolve_transductions(insertions)

    tebreak.text_summary(insertions, outfile=args.detail_out) # debug

    te_insertions = [TEIns(ins) for ins in insertions]

    if len(te_insertions) > 0: print te_insertions[0].header()

    for te_ins in sorted(te_insertions):
        if te_ins.pass_filter(forest): print te_ins


if __name__ == '__main__':
    # set up logger
    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(format=FORMAT)
 
    parser = argparse.ArgumentParser(description='Resolve TE insertions from TEbreak data')
    parser.add_argument('-p', '--pickle', required=True)
    parser.add_argument('-t', '--processes', default=1, help='split work across multiple processes')
    parser.add_argument('-i', '--inslib_fasta', required=True, help='reference for insertions (not genome)')

    parser.add_argument('-m', '--filter_bed', default=None, help='MAF of -i/--inslib_fasta lastal vs. -d/--filter_ref_last_db')
    parser.add_argument('--max_bam_count', default=0)

    parser.add_argument('--refoutdir', default='tebreak_refs')

    parser.add_argument('--tmpdir', default='/tmp')
    parser.add_argument('--skip_align', action='store_true', default=False)
    parser.add_argument('-v', '--verbose', action='store_true', default=False)
    parser.add_argument('--keep_tmp_bams', action='store_true', default=False)
    parser.add_argument('--detail_out', default='resolve.out', help='file to write detailed output')
    args = parser.parse_args()
    main(args)
