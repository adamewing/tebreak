#!/usr/bin/env python

import os
import shutil
import cPickle as pickle
import argparse
import logging
import tebreak
import pysam
import subprocess

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
        self.end3 = None
        self.end5 = None

        self._fill_out()

    def _fill_out(self):
        self.out['Chromosome']     = self.ins['chrom']
        self.out['Left_Junction']  = min(self.junctions())
        self.out['Right_Junction'] = max(self.junctions())
        self.assign35ends()
        self.te_family()
        self.elt_coords()
        self.out['TE_Match_Pct'] = self.be_avg_pctmatch()

    def assign35ends(self):
        self.out['3_Prime_End'] = 'NA'
        self.out['5_Prime_End'] = 'NA'

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

        if None not in (self.end3, self.end5):
            if self.end3 + '_bestmatch' in self.ins: self.out['TE_Align_End']   = self.ins[self.end3 + '_bestmatch'].target_end
            if self.end5 + '_bestmatch' in self.ins: self.out['TE_Align_Start'] = self.ins[self.end5 + '_bestmatch'].target_start

    def junctions(self):
        j = []
        if 'be1_breakpos' in self.ins: j.append(self.ins['be1_breakpos'])
        if 'be2_breakpos' in self.ins: j.append(self.ins['be2_breakpos'])
        return j

    def be_avg_pctmatch(self):
        m = []
        if 'be1_bestmatch' in self.ins: m.append(self.ins['be1_bestmatch'].pct_match())
        if 'be2_bestmatch' in self.ins: m.append(self.ins['be2_bestmatch'].pct_match())
        if len(m) == 0: return 0.0
        return float(sum(m)) / float(len(m))


    def polyA_filter(self):
        ''' return true if only supported by poly-A matches '''
        found_non_pA = False

        if 'be1_bestmatch' in self.ins:
            if poly_A_frac(self.ins['be1_bestmatch'].target_align) < 0.95: found_non_pA = True

        if 'be2_bestmatch' in self.ins:
            if poly_A_frac(self.ins['be2_bestmatch'].target_align) < 0.95: found_non_pA = True

        return not found_non_pA

    def genome_location_filter(self, forest):
        if self.ins['chrom'] in forest:
            hits = forest[self.ins['chrom']].find(self.out['Left_Junction'], self.out['Right_Junction']+1)
            for hit in hits:
                superfamily = hit.query_id.split(':')[0]
                if superfamily in self.out['Superfamily'].split(','):
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

        return passed

    def header(self):
        return '\t'.join(self.out.keys())

    def __str__(self):
        return '\t'.join(map(str, self.out.values()))


#######################################
## Functions                         ##
#######################################


def best_match(last_results, query_name, req_target=None):
    qres = [res for res in sorted(last_results) if res.query_id == query_name and (req_target is None or req_target == res.target_id)]
    if len(qres) > 0:
        return qres[0]

    return None


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
        p = subprocess.Popen(['bwa', 'index', ref_fa], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        for line in p.stdout: pass
        assert os.path.exists(ref_fa + '.bwt'), 'could not bwa index %s' % ref_fa

    if not os.path.exists(ref_fa + '.tis') and makeLAST:
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

        # realign all unmapped sequences
        if 'be1_umap_seq' in ins['INFO'] and ins['INFO']['be1_umap_seq'] is not None:
            for i, umap_seq in enumerate(ins['INFO']['be1_umap_seq'].split(',')):
                fa.write('>%s|%s\n%s\n' % (ins['INFO']['be1_obj_uuid'],str(i), umap_seq))

        if 'be2_umap_seq' in ins['INFO'] and ins['INFO']['be2_umap_seq'] is not None:
            for i, umap_seq in enumerate(ins['INFO']['be2_umap_seq'].split(',')):
                fa.write('>%s|%s\n%s\n' % (ins['INFO']['be2_obj_uuid'], str(i), umap_seq))

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


def last_interval_forest(maf_file, padding=100):
    ''' build dictionary of interval trees; value=LASTResult '''
    forest = dd(Intersecter)

    maf_lines   = []
    maf_results = []

    with open(maf_file, 'r') as maf:
        for line in maf:
            if not line.startswith('#'):
                if line.strip() != '':
                    maf_lines.append(line.strip())

                else:
                    maf_results.append(tebreak.LASTResult(maf_lines))
                    maf_lines = []

    for res in maf_results:
        forest[res.target_id].add_interval(Interval(res.target_start, res.target_end, value=res))

    return forest


def poly_A_frac(seq):
    ''' home much of the sequence is A '''
    if seq is None: return 0.0

    f = []
    for b in ('A', 'T'):
        f.append(len([base for base in list(seq.upper()) if base == b]))

    return float(max(f)) / float(len(seq))


def add_insdata(ins, last_res):
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


def make_inslib_mask(ref_fa, lastdb, refoutdir):
    outmaf = refoutdir + '/' + '.'.join(os.path.basename(ref_fa).split('.')[:-1]) + '.last.maf'

    if os.path.exists(outmaf):
        sys.stderr.write('inslib mask MAF already exists, using %s\n' % outmaf)
        return outmaf

    cmd = ['lastal', '-e', '100', lastdb, ref_fa]

    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    with open(outmaf, 'w') as maf:
        for line in p:
            maf.write(line)

    return outfn


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


def identify_transductions(ins, minmapq=10):
    bedict = {'be1': None, 'be2': None}

    if 'be1_bestmatch' in ins['INFO']: bedict['be1'] = ins['INFO']['be1_bestmatch']
    if 'be2_bestmatch' in ins['INFO']: bedict['be2'] = ins['INFO']['be2_bestmatch']

    for be in ('be1','be2'):
        if bedict[be] is not None:
            tr_seqs = []
            tr_locs = [] # chrom, start, end, 3p or 5p / unmap

            num_segs = 0
            if ins['INFO'][be+'_dist_seq'] is not None: num_segs = len(ins['INFO'][be+'_dist_seq'].split(','))
            if num_segs > 1:
                for distnum in range(num_segs):
                    seg_mapq = int(ins['INFO'][be+'_dist_mpq'].split(',')[distnum])
                    if distnum != bedict[be].query_distnum and seg_mapq >= minmapq:
                        tr_chrom = ins['INFO'][be+'_dist_chr'].split(',')[distnum]
                        tr_start = ins['INFO'][be+'_dist_pos'].split(',')[distnum]
                        tr_end   = ins['INFO'][be+'_dist_end'].split(',')[distnum]
                        tr_side  = '5p'

                        if ins['INFO'][be+'_is_3prime']: tr_side = '3p'

                        tr_seqs.append(ins['INFO'][be+'_dist_seq'].split(',')[distnum])
                        tr_locs.append((tr_chrom, tr_start, tr_end, tr_side))

            if num_segs > 0 and ins['INFO'][be+'_umap_seq'] is not None: # unmapped (maybe) transduced seqs
                for unmap_seq in ins['INFO'][be+'_umap_seq'].split(','):
                    tr_side  = '5p'
                    if ins['INFO'][be+'_is_3prime']: tr_side = '3p'

                    tr_seqs.append(unmap_seq)
                    tr_locs.append(('unmap',0,0,tr_side))

            if len(tr_seqs) > 0:
                ins['INFO'][be+'_trans_seq'] = tr_seqs
                ins['INFO'][be+'_trans_loc'] = tr_locs

    return ins


def resolve_insertion(args, ins, inslib_fa):
    ''' add data based on alignments of library to consensus '''
    last_res = last_alignment(ins, inslib_fa)
    ins = add_insdata(ins, last_res)

    if 'best_ins_matchpct' in ins['INFO'] and not args.skip_align:
        if ins['INFO']['dr_count'] > 0 and ins['INFO']['best_ins_matchpct'] > 0.90: # change to parameter
            tmp_bam = remap_discordant(ins, inslib_fa=inslib_fa, tmpdir=args.tmpdir)

            if tmp_bam is not None:
                bam = pysam.AlignmentFile(tmp_bam, 'rb')
                ins['INFO']['support_bam_file'] = tmp_bam
                ins['INFO']['mapped_target'] = bam.mapped

    if 'best_ins_matchpct' in ins['INFO']: ins = identify_transductions(ins)

    return ins


def resolve_unknown(insertions):
    ''' don't bother with alignment to consensus, attempt to assemble unknown insertions '''
    pass


def main(args):
    if args.verbose: logger.setLevel(logging.DEBUG)
    insertions = []

    with open(args.pickle, 'r') as pickin:
        insertions = pickle.load(pickin)

    inslib_fa = prepare_ref(args.inslib_fasta, refoutdir=args.refoutdir, makeFAI=False, makeBWA=False)

    forest = None
    # set up reference mask MAF if available
    if args.filter_ref_maf is not None:
        logger.debug('using MAF as genome reference filter: %s' % args.filter_ref_maf)
        forest = last_interval_forest(args.filter_ref_maf)

    elif args.filter_ref_last_db is not None:
        logger.debug('building genome reference filter MAF from: %s vs %s' %(args.inslib_fasta, args.filter_ref_last_db))
        assert os.path.exists(args.filter_ref_last_db + '.suf'), 'reference %s not indexed with lastdb' % args.filter_ref_last_db
        filter_ref_maf = make_inslib_mask(args.inslib_fasta, args.filter_ref_last_db, args.refoutdir)

    results = []

    pool = mp.Pool(processes=int(args.processes))

    for ins in insertions:
        res = pool.apply_async(resolve_insertion, [args, ins, inslib_fa])
        results.append(res)

    insertions = [res.get() for res in results]

    tebreak.text_summary(insertions, outfile=args.detail_out) # debug

    for ins in insertions:
        te_ins = TEIns(ins)
        if te_ins.pass_filter(forest): print te_ins


if __name__ == '__main__':
    # set up logger
    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(format=FORMAT)
 
    parser = argparse.ArgumentParser(description='Resolve TE insertions from TEbreak data')
    parser.add_argument('-p', '--pickle', required=True)
    parser.add_argument('-t', '--processes', default=1, help='split work across multiple processes')
    parser.add_argument('-i', '--inslib_fasta', required=True, help='reference for insertions (not genome)')

    parser.add_argument('-d', '--filter_ref_last_db', default=None, help='reference LAST db')
    parser.add_argument('-m', '--filter_ref_maf', default=None, help='MAF of -i/--inslib_fasta lastal vs. -d/--filter_ref_last_db')

    parser.add_argument('--refoutdir', default='tebreak_refs')

    parser.add_argument('--tmpdir', default='/tmp')
    parser.add_argument('--skip_align', action='store_true', default=False)
    parser.add_argument('-v', '--verbose', action='store_true', default=False)
    parser.add_argument('--detail_out', default='tebreak.out', help='file to write detailed output')
    args = parser.parse_args()
    main(args)