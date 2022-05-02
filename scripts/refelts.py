#!/usr/bin/env python

import os
import sys
import pysam
import re
import subprocess
import argparse

import multiprocessing as mp

import scipy.stats as ss
import numpy as np

import skbio.alignment as skalign
import skbio.sequence as skseq

from collections import Counter
from uuid import uuid4
from time import sleep

import logging
FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


header = [
'Chromosome',
'RMSK_Start',
'RMSK_End',
'Orientation',
'Name',
'Junc_5p',
'Junc_3p',
'TSD_seq',
'Empty_Site_Consensus',
'AltCount_5p',
'AltCount_3p',
'RefCount_5p',
'RefCount_3p',
'VAF'
]


class Block:
    def __init__(self, tstart, tend):
        self.tstart = min(int(tstart), int(tend))
        self.tend   = max(int(tstart), int(tend))

    def length(self):
        return self.tend - self.tstart

    def __str__(self):
        return str(self.tstart) + ' ' +  str(self.tend)

    def __lt__(self, other):
        return self.tstart < other.tstart


class PSL:
    def __init__(self, rec, consensus):
        # psl spec: http://www.ensembl.org/info/website/upload/psl.html
        (self.matches, self.misMatches, self.repMatches, self.nCount, self.qNumInsert, self.qBaseInsert, 
         self.tNumInsert, self.tBaseInsert, self.strand, self.qName, self.qSize, self.qStart, self.qEnd,
         self.tName, self.tSize, self.tStart, self.tEnd, self.blockCount, self.blockSizes, self.qStarts,
         self.tStarts) = rec.strip().split()

        self.cons = consensus
        self.rec = rec.strip()

        self.tBlocks = []
        for bsize, tstart in zip(self.blockSizes.split(',')[:-1], self.tStarts.split(',')[:-1]): # [:-1] due to trailing comma
            self.tBlocks.append(Block(int(tstart), int(tstart) + int(bsize)))

        self.tBlocks.sort()

        self.tName = self.tName.replace('chr', '')

        self.tStart, self.tEnd, self.qStart, self.qEnd = map(int, (self.tStart, self.tEnd, self.qStart, self.qEnd))
        self.qSize, self.tSize = map(int, (self.qSize, self.tSize))
        
        if self.qStart > self.qEnd:
            self.qStart, self.qEnd = self.qEnd, self.qStart

        if self.tStart > self.tEnd:
            self.tStart, self.tEnd = self.tEnd, self.tStart

    def match(self, chrom, pos, window=0):
        ''' return True if chrom:pos intersects BLAT hit +/- window '''
        chrom = chrom.replace('chr', '')
        if chrom != self.tName:
            return False

        if int(pos) >= int(self.tStart)-window and int(pos) <= int(self.tEnd)+window:
            return True

        return False

    def refspan(self):
        ''' return footprint of match relative to referece genome '''
        return self.tEnd - self.tStart

    def score(self):
        ''' adapted from https://genome.ucsc.edu/FAQ/FAQblat.html#blat4 '''
        return (int(self.matches) + (int(self.repMatches)>>1)) - int(self.misMatches) - int(self.qNumInsert) - int(self.tNumInsert)

    def pctmatch(self):
        ''' adapted from https://genome.ucsc.edu/FAQ/FAQblat.html#blat4 '''
        qAliSize = int(self.qEnd) - int(self.qStart)
        tAliSize = int(self.tEnd) - int(self.tStart)
        if min(qAliSize, tAliSize) <= 0:
            return 0.0

        sizeDif = abs(qAliSize - tAliSize)
        total = int(self.matches) + int(self.repMatches) + int(self.misMatches)

        if total > 0:
            return 1.0-float((int(self.misMatches) + int(self.qNumInsert) + int(self.tNumInsert) + round(3*log(1+sizeDif)))) / float(total)

        return 0.0

    def __lt__(self, other):
        ''' used for ranking BLAT hits '''
        return self.score() > other.score()

    def __str__(self):
        return self.rec


class SortableRead:
    def __init__(self, read):
        self.read = read
        self.seq  = read.seq
        self.seqstart = read.reference_start-read.query_alignment_start

    def __gt__(self, other):
        if self.read.tid == other.read.tid:
            return self.seqstart > other.seqstart
        else:
            return self.read.tid > other.read.ti


class SplitRead:
    ''' store information about split read alignment '''
    def __init__(self, chrom, read, bamfn, minqual):
        self.uuid  = str(uuid4())
        self.chrom = chrom
        self.read  = read
        self.bamfn = os.path.basename(bamfn)
        self.minqual = minqual

        self.cliplen = len(read.seq) - len(read.query_alignment_sequence)

        self.breakleft  = False
        self.breakright = False
        self.breakpos   = None
 
        if read.qstart < read.rlen - read.qend:
            self.breakpos   = read.get_reference_positions()[-1] # breakpoint on right
            self.breakright = True
 
        else:
            self.breakpos  = read.get_reference_positions()[0] # breakpoint on left
            self.breakleft = True
 
        assert self.breakpos is not None
        assert self.breakleft != self.breakright
 
    def getRG(self):
        ''' return read group from RG aux tag '''
        for tag, val in self.read.tags:
            if tag == 'RG': return val
        return None

    def __gt__(self, other):
        ''' enables sorting of SplitRead objects '''
        if self.chrom == other.chrom:
            return self.breakpos > other.breakpos
        else:
            return self.chrom > other.chrom
 
    def __str__(self):
        dir = 'left'
        if self.breakright: dir='right'
        return ' '.join(map(str, ('SplitRead:', self.chrom, self.breakpos, self.cliplen, self.read.qname, dir)))


class ReadCluster:
    ''' parent class for read clusters '''
    def __init__(self, firstread=None):
        self.uuid   = str(uuid4())
        self.reads  = []
        self.start  = 0
        self.end    = 0
        self.median = 0
        self.chrom  = None

        if firstread is not None:
            self.add_read(firstread)

    def add_read(self, r):
        ''' add a read and update '''
        self.reads.append(r)
 
        if self.chrom is None: self.chrom = r.chrom
 
        assert self.chrom == r.chrom # clusters can't include > 1 chromosome
 
        ''' update statistics '''
        positions = []
        positions += [pos for r in self.reads for pos in r.read.positions]

        self.reads.sort()
        self.start  = max(positions)
        self.end    = min(positions)
        self.median = int(np.median(positions))

    def readgroups(self):
        c = Counter([r.getRG() for r in self.reads])
        return [str(k[0]) + '|' + str(k[1]) for k in zip(c.keys(), c.values())]

    def bamfiles(self):
        c = Counter([r.bamfn for r in self.reads])
        return [str(k[0]) + '|' + str(k[1]) for k in zip(c.keys(), c.values())]

    def find_extrema(self):
        ''' return leftmost and rightmost aligned positions in cluster vs. reference '''
        positions = []
        positions += [pos for r in self.reads for pos in r.read.positions]
        return min(positions), max(positions)
 
    def avg_matchpct(self):
        return np.mean([read_matchpct(r.read) for r in self.reads])
 
    def __len__(self):
        return len(self.reads)


class SplitCluster(ReadCluster):
    ''' store and manipulate groups of SplitRead objects '''

    def add_splitread(self, sr):
        ''' add a SplitRead and update '''
        self.reads.append(sr)
 
        if self.chrom is None: self.chrom = sr.chrom
 
        assert self.chrom == sr.chrom # clusters can't include > 1 chromosome
 
        ''' update statistics '''
        self.reads.sort()
        self.start  = self.reads[0].breakpos
        self.end    = self.reads[-1].breakpos
        self.median = self.reads[int(len(self)/2)].breakpos
 
    def subcluster_by_breakend(self, breakends, direction='both'):
        ''' return a new cluster containing only reads with breakpoints in passed list '''
        new = SplitCluster()
        assert direction in ('both', 'left', 'right')
        
        if direction == 'both':
            [new.add_splitread(sr) for sr in self.reads if sr.breakpos in breakends]
 
        if direction == 'left':
            [new.add_splitread(sr) for sr in self.reads if sr.breakpos in breakends and sr.breakleft]
 
        if direction == 'right':
            [new.add_splitread(sr) for sr in self.reads if sr.breakpos in breakends and sr.breakright]
 
        return new

    def consensus(self, minscore = 0.9):
        ''' build consensus from sorted aligned reads iteratively '''

        S = -np.ones((256, 256)) + 2 * np.identity(256)
        S = S.astype(np.int16)

        minqual = self.reads[0].minqual

        sortable_reads = [SortableRead(sr.read) for sr in self.reads]
        seqs = [qualtrim(sorted_read.read, minqual=minqual) for sorted_read in sorted(sortable_reads)]
        seqs = [s for s in seqs if len(s) > 20]

        if len(seqs) == 0:
            return '', 0.0

        if len(seqs) == 1: # no consensus necessary
            return seqs[0], 1.0

        uniq_seqs = [seqs[0]]
        for i, seq in enumerate(seqs[1:], start=1):
            if seq != seqs[i-1]:
                uniq_seqs.append(seq)

        if len(uniq_seqs) == 1: # all seqs were the same!
            return uniq_seqs[0], 1.0

        cons = uniq_seqs[0]
        scores = []

        if len(uniq_seqs) > 1000:
            uniq_seqs = [uniq_seqs[u] for u in sorted(np.random.choice(range(len(uniq_seqs)), size=1000))]

        for seq in uniq_seqs[1:]:

            cons = cons.replace('N', 'A')
            seq = seq.replace('N', 'A')

            s1 = skseq.DNA(cons)
            s2 = skseq.DNA(seq)

            try:
                aln_res = skalign.local_pairwise_align_ssw(s1, s2)
            except:
                return cons, 0.0

            aln_tab = aln_res[0]

            s1_aln, s2_aln = aln_res[2]

            a1 = cons[s1_aln[0]:s1_aln[1]+1]

            score = 0.0
            if aln_tab.shape.position > 10: # param?
                score = sum(aln_tab.conservation(gap_mode='include')==1.)/aln_tab.shape.position

            if re.search(a1, cons):
                cons_start, cons_end = s1_aln[0], s1_aln[1]+1

                if score >= minscore and cons_end > len(cons)-5:
                    scores.append(score)
                    align_end = s2_aln[1]+1 
                    cons += seq[align_end:]

        if scores:
            return cons, np.mean(scores)

        else:
            return cons, 0.0

    def all_breakpoints(self):
        ''' returns uniquified list of breakpoints '''
        return list(set([read.breakpos for read in self.reads]))
 
    def median_D(self):
        return np.median([splitqual(sr.read) for sr in self.reads])

    def min_cliplen(self):
        return min([sr.cliplen for sr in self.reads])

    def max_cliplen(self):
        return max([sr.cliplen for sr in self.reads])
 
    def __str__(self):
        break_count = Counter([read.breakpos for read in self.reads])
        return '\t'.join(map(str, ('SplitCluster:', self.chrom, self.start, self.end, len(self.reads), break_count)))


class BreakEnd:
    ''' coallate information about a breakend '''
    def __init__(self, chrom, breakpos, cluster, consensus, score, direction):
        self.uuid      = str(uuid4())
        self.cluster   = cluster
        self.chrom     = chrom
        self.breakpos  = breakpos
        self.consensus = consensus
        self.consscore = score
        self.direction = direction
 
    def __len__(self):
        return len(self.cluster)

    def __str__(self):
        return '%s:%d:%s:%s' % (self.chrom, self.breakpos, self.consensus, self.direction)


def splitqual(read):
    ''' return Mann-Whitney P for clipped vs unclipped quals '''
    
    breakpos = None
 
    breakpos = read.get_aligned_pairs()[-1][0] # breakpoint on right
 
    q1 = list(map(ord, list(read.qual[:breakpos])))
    q2 = list(map(ord, list(read.qual[breakpos:])))

    if min(q1) == max(q1) == min(q2) == max(q2):
        return 1.0
 
    return ss.mannwhitneyu(q1, q2)[1]


def guess_minqual(bam):
    minscore = None
    n = 0

    for read in bam.fetch():
        n += 1
        m = min([ord(q) for q in list(read.qual)])
        if minscore is None or minscore > m:
            minscore = m

        if n > 10000: break

    return minscore


def fetch_clipped_reads(bams, chrom, start, end):
    ''' Return list of SplitRead objects '''
 
    splitreads = []
 
    start = int(start)
    end   = int(end)
 
    assert start < end

    if start < 0: start = 0

    minqual = guess_minqual(bams[0]) # used for quality trimming when building consensus

    for bam in bams:
        if chrom not in bam.references:
            return splitreads


    for bam in bams:
        for read in bam.fetch(chrom, start, end):

            if not read.is_unmapped and not read.is_duplicate and read.mapq > 0:
                if read.rlen - read.alen >= 3: # 'soft' clipped?
     
                    # length of 'minor' clip
                    altclip = min(read.qstart, read.rlen-read.qend)

                    # junk bases
                    N_count = 0
                    if 'N' in read.seq: N_count = Counter(read.seq)['N']
     
                    if altclip <= 2: # could add as a filter
                        if N_count <= 2 and splitqual(read) >= 0.01:
                            chrom = str(bam.getrname(read.tid))

                            if len(read.get_reference_positions()) > 0:
                                splitreads.append(SplitRead(chrom, read, bam.filename, minqual))

    return splitreads


def build_sr_clusters(splitreads, searchdist=100): # TODO PARAM, 
    ''' cluster SplitRead objects into Cluster objects and return a list of them '''
    clusters  = []
 
    for sr in splitreads:
        if len(clusters) == 0:
            clusters.append(SplitCluster(sr))
 
        elif clusters[-1].chrom != sr.chrom:
            clusters.append(SplitCluster(sr))
 
        else:
            if abs(clusters[-1].median - sr.breakpos) > searchdist:
                clusters.append(SplitCluster(sr))
 
            else:
                clusters[-1].add_splitread(sr)

    return clusters


def build_breakends(cluster, tmpdir='/tmp'):
    ''' returns list of breakends from cluster '''
    breakends = []

    for breakpos in cluster.all_breakpoints():
        for dir in ('left', 'right'):
            subcluster = cluster.subcluster_by_breakend([breakpos], direction=dir)
            if len(subcluster) >= 4 and subcluster.max_cliplen() >= 10:
                seq     = subcluster.reads[0].read.seq
                score   = 1.0

                if len(subcluster) > 1: seq, score = subcluster.consensus()
 
                N_count = 0
                if 'N' in seq: N_count = Counter(seq)['N']

                if seq != '' and score >= 0.9 and N_count <= 3:
                    breakends.append(BreakEnd(cluster.chrom, breakpos, subcluster, seq, score, dir))
 
    return breakends


def qualtrim(read, minqual=35):
    ''' return quality-trimmed sequence given a pysam.AlignedSegment '''
    q = [ord(b)-minqual for b in list(read.qual)]

    for i in range(0,len(q)-4): # sliding window, 4bp
        if np.mean(q[i:i+4]) < 5:
            return read.seq[:i]

    return read.seq


def locate_subseq(longseq, shortseq):
    ''' return (start, end) of shortseq in longseq '''
    assert len(longseq) >= len(shortseq), 'orient_subseq: %s < %s' % (longseq, shortseq)
 
    match = re.search(shortseq, longseq)
    if match is not None:
        return sorted((match.start(0), match.end(0)))
 
    return None


def start_blat_server(blatref, port=9999):
    # parameters from https://genome.ucsc.edu/FAQ/FAQblat.html#blat5

    server_up = True

    dummy_fa  = '/tmp/' + str(uuid4()) + '.fa'
    dummy_psl = dummy_fa.replace('.fa', '.psl')

    with open(dummy_fa, 'w') as dout:
        dout.write('>\n' + 'A'*100)

    poll_cmd = ['gfClient', 'localhost', str(port), os.path.dirname(blatref), dummy_fa, dummy_psl]
    poll_time = 10

    t = subprocess.Popen(poll_cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    for line in t.stderr:
        line = line.decode()
        if line.startswith('Sorry'):
            server_up = False
            logger.info("No BLAT server found, starting one up...")
            cmd = ['gfServer', 'start', 'localhost', str(port), '-stepSize=5', blatref]
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        else:
            logger.info("Found BLAT server!")



    while not server_up:
        started = True
        t = subprocess.Popen(poll_cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE)

        for line in t.stderr:
            line = line.decode()
            if line.startswith('Sorry'):
                started = False

        if not started:
            logger.info("waiting for BLAT server to start...")
            sleep(poll_time)
        else:
            server_up=True
            logger.info("BLAT for %s server up" % blatref)

    return port


def blat(fasta, outpsl, blatrefdir, port=9999, minScore=0, maxIntron=None):
    ''' BLAT using gfClient utility '''
    cmd  = ['gfClient', 'localhost', str(port), '-nohead']

    if maxIntron is not None:
        cmd.append('-maxIntron=' + str(maxIntron))

    if minScore is not None:
        cmd.append('-minScore=' + str(minScore))

    cmd += [blatrefdir, fasta, outpsl]

    p = subprocess.Popen(cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    for line in p.stdout:
        pass


def eval_break(blatrefdir, breakend, direct, elt_chrom, elt_start, elt_end):
    if breakend.direction != direct:
        return False

    out_psl = None

    fa_tmp = '/tmp/' + str(uuid4()) + '.fa'
    psl_tmp = '/tmp/' + str(uuid4()) + '.psl'

    with open(fa_tmp, 'w') as fa_out:
        fa_out.write('>%s_%d\n%s\n' % (breakend.chrom, breakend.cluster.start, breakend.consensus))

    blat(fa_tmp, psl_tmp, blatrefdir)

    with open(psl_tmp, 'r') as psl:
        for line in psl:
            rec = PSL(line.strip(), breakend.consensus)
            if float(rec.matches) / len(breakend.consensus) > 0.9:
                # BLAT strips 'chr' prefix from PSL output
                if elt_chrom.lstrip('chr') == rec.tName and int(rec.tStart) < elt_start + 50 and int(rec.tEnd) > elt_end-50:
                    out_psl = rec

        os.remove(fa_tmp)
        os.remove(psl_tmp)

        return out_psl


def tsd(psl, ref, b_left_init=0, b_right_init=0, max_iter=100):
    b_left_pos  = b_left_init
    b_right_pos = b_right_init

    if int(psl.blockCount) < 2:
        return 0,0,0,0,'NA'

    # pick largest gap between blocks
    gap = 0
    for i, block in enumerate(psl.tBlocks[1:], 1):
        prev_block = psl.tBlocks[i-1]

        if block.tstart - prev_block.tend > gap:
            gap = block.tstart - prev_block.tend
            b_left_pos  = b_left_init + prev_block.tend
            b_right_pos = b_right_init + block.tstart


    chrom = 'chr' + psl.tName

    b_left  = b_left_pos 
    b_right = b_right_pos

    if b_left < 0:
        b_left = 0

    if b_right < 0:
        b_right = 0

    nt_l = ref.fetch(chrom, b_left, b_left+1)
    nt_r = ref.fetch(chrom, b_right, b_right+1)

    if nt_l != nt_r:
        return b_left, b_left, b_right, b_right, 'NA'

    else:
        i = 0
        while nt_l == nt_r and i < max_iter:
            b_left  -= 1
            b_right -= 1

            if b_left < 0:
                b_left = 0

            if b_right < 0:
                b_right = 0

            i += 1

            nt_l = ref.fetch(chrom, b_left, b_left+1)
            nt_r = ref.fetch(chrom, b_right, b_right+1)

    if i >= max_iter:
        return b_left, b_left, b_right, b_right, 'NA'

    l_b_start = b_left+1
    r_b_start = b_right+1

    b_left  = b_left_pos
    b_right = b_right_pos

    nt_l = ref.fetch(chrom, b_left, b_left+1)
    nt_r = ref.fetch(chrom, b_right, b_right+1)

    i = 0

    while nt_l == nt_r and i < max_iter:
        b_left += 1
        b_right += 1

        if b_left < 0:
            b_left = 0

        if b_right < 0:
            b_right = 0

        i += 1

        nt_l = ref.fetch(chrom, b_left, b_left+1)
        nt_r = ref.fetch(chrom, b_right, b_right+1)

    if i >= max_iter:
        return b_left, b_left, b_right, b_right, 'NA'

    l_b_end = b_left
    r_b_end = b_right

    if l_b_start < 0:
        l_b_start = 0

    tsd_seq = ref.fetch(chrom, l_b_start, l_b_end)

    if tsd_seq:
        return l_b_start, l_b_end, r_b_start, r_b_end, tsd_seq

    else:
        return l_b_left, l_b_left, r_b_start, r_b_end, 'NA'


def getVAF(bam, chrom, poslist):
    ''' return number of reads supporting alt (insertion), ref (reference) and vaf (variant allele fraction) '''
    poslist = map(int, poslist)
    alt, ref = break_count(bam, chrom, poslist)
    vaf = 0.0 

    if float(ref+alt) > 0:
        vaf = float(alt)/float(alt+ref)

    return alt, ref, vaf


def break_count(bam, chrom, poslist, minpad=5, flex=1, minmapq=10):
    ''' ref = number of reads spanning TSD, alt = number of reads clipped at breakpoint in poslist '''
    altcount = 0
    refcount = 0
    discards = 0

    poslist = list(poslist)

    tsd_start = min(poslist)
    tsd_end   = max(poslist)

    tsd_len = tsd_end - tsd_start

    if tsd_start < minpad: tsd_start = minpad

    for read in bam.fetch(chrom, tsd_start-minpad, tsd_end+minpad):
        if read.is_unmapped or read.is_duplicate:
            continue

        if read.mapq < minmapq:
            continue

        rclip = len(read.seq) - read.query_alignment_end 
        lclip = read.query_alignment_start

        rbreak = 0
        lbreak = 0

        if rclip > max(tsd_len, minpad):
            rbreak = read.reference_end

        if lclip > max(tsd_len, minpad):
            lbreak = read.reference_start

        support_alt = False

        for pos in poslist: # does this read support a breakpoint in the list?
            if (rbreak >= pos-flex and rbreak <= pos+flex) or (lbreak >= pos-flex and lbreak <= pos+flex):
                support_alt = True

        if support_alt:
            altcount += 1

        else:
            #for pos in poslist: # does this read span a breakpoint in the list?
            if read.alen == len(read.seq):
                if read.reference_start < tsd_start and read.reference_end > tsd_end: # span TSD
                    refcount += 1


    return altcount, refcount


def junc_search(bams, chrom, pos, end, window=200, minmapq=10):

    junc_pos = []

    for bam in bams:
        for read in bam.fetch(chrom, pos-window, pos+window):
            if read.mapq < minmapq:
                continue

            if end == 'L':
                if read.cigartuples[-1][0] == 4:
                    junc_pos.append(read.reference_end)

            if end == 'R':
                if read.cigartuples[0][0] == 4:
                    junc_pos.append(read.reference_start)

    if len(junc_pos) < 4:
        logger.warning('no good junction candidate: %s %d %s' % (chrom, pos, end))
        return pos

    c = Counter(junc_pos)

    best_junc = c.most_common()[0][0]
    best_count = c.most_common()[0][1]

    if best_count/len(junc_pos) < .5:
        logger.warning('no good junction candidate: %s %d %s' % (chrom, pos, end))
        return pos

    return best_junc


def support(bams, chrom, pos, end, tsd_len, margin=5, minmapq=10):
    alt_count = 0
    ref_count = 0

    for bam in bams:
        for read in bam.fetch(chrom, pos-margin-tsd_len, pos+margin+tsd_len):
            if read.mapq < minmapq:
                continue

            if not read.is_unmapped and not read.is_duplicate:
                if read.alen == read.rlen: # not clipped:
                    if read.reference_start < pos-margin and read.reference_end > pos+margin:
                        ref_count += 1

                else:
                    if end == 'L':
                        if read.cigartuples[-1][0] == 4:
                            #print(end, pos, pos-margin-tsd_len, pos+margin+tsd_len, read.reference_end)
                            if pos-margin-tsd_len < read.reference_end and pos+margin+tsd_len > read.reference_end:
                                alt_count += 1

                    if end == 'R':
                        if read.cigartuples[0][0] == 4:
                            #print(end, pos, pos-margin-tsd_len, pos+margin+tsd_len, read.reference_start)
                            if pos-margin-tsd_len < read.reference_start and pos+margin+tsd_len > read.reference_start:
                                alt_count += 1

    return alt_count, ref_count


def ref_no_junc(args, chrom, start, end, orient, name):
    bams = [pysam.AlignmentFile(bam) for bam in args.bam.split(',')]
    ref = pysam.Fastafile(args.fastaref)


    l_altcount, l_refcount = support(bams, chrom, int(start), 'L', 0)
    r_altcount, r_refcount = support(bams, chrom, int(end), 'R', 0)

    vaf = 0.0
    if l_altcount + r_altcount + l_refcount + r_refcount > 0:
        vaf = float(l_altcount + r_altcount) / float(l_altcount + r_altcount + l_refcount + r_refcount)

    altcount_5p = l_altcount
    altcount_3p = r_altcount

    refcount_5p = l_refcount
    refcount_3p = r_refcount

    junc_5p = start
    junc_3p = end

    if orient == '-':
        junc_5p, junc_3p = junc_3p, junc_5p
        altcount_5p, altcount_3p = altcount_3p, altcount_5p
        refcount_5p, refcount_3p = refcount_3p, refcount_5p

    if args.persample is not None:
        vaf = []

        altcount_5p = []
        altcount_3p = []
        refcount_5p = []
        refcount_3p = []

        with open(args.persample) as samples:
            for line in samples:
                sbamfn, sname = line.strip().split()
                sbam = [pysam.AlignmentFile(sbamfn)]

                l_altcount, l_refcount = support(sbam, chrom, start, 'L', 0) 
                r_altcount, r_refcount = support(sbam, chrom, end, 'R', 0)

                svaf = 0.0
                if l_altcount + r_altcount + l_refcount + r_refcount > 0:
                    svaf = float(l_altcount + r_altcount) / float(l_altcount + r_altcount + l_refcount + r_refcount)

                vaf.append('%s|%.3f' % (sname, svaf))

                s_altcount_5p = l_altcount
                s_altcount_3p = r_altcount
                s_refcount_5p = l_refcount
                s_refcount_3p = r_refcount

                if orient == '-':
                    s_altcount_5p, s_altcount_3p = s_altcount_3p, s_altcount_5p
                    s_refcount_5p, s_refcount_3p = s_refcount_3p, s_refcount_5p

                #print(sbamfn, l_altcount, l_refcount, r_altcount, r_refcount)

                altcount_5p.append('%s|%d' % (sname, s_altcount_5p))
                altcount_3p.append('%s|%d' % (sname, s_altcount_3p))
                refcount_5p.append('%s|%d' % (sname, s_refcount_5p))
                refcount_3p.append('%s|%d' % (sname, s_refcount_3p))

        vaf = ','.join(vaf)
        altcount_5p = ','.join(altcount_5p)
        altcount_3p = ','.join(altcount_3p)
        refcount_5p = ','.join(refcount_5p)
        refcount_3p = ','.join(refcount_3p)

    out = (chrom, start, end, orient, name, junc_5p, junc_3p, 'NA', 'NA', altcount_5p, altcount_3p, refcount_5p, refcount_3p, vaf)
    return out


def ref_ins(args, chrom, start, end, orient, name):
    bams = [pysam.AlignmentFile(bam) for bam in args.bam.split(',')]
    ref = pysam.Fastafile(args.fastaref)

    # Find the junction

    start_splits = fetch_clipped_reads(bams, chrom, start-25, start+25)
    start_splits.sort()
    start_clusters = build_sr_clusters(start_splits)
    start_breaks = [build_breakends(c) for c in start_clusters]

    psl_rec = None

    for be in start_breaks:
        for b in be:
            psl_rec = eval_break(os.path.dirname(args.blatref), b, 'right', chrom, start, end)

    if psl_rec is None:
        end_splits = fetch_clipped_reads(bams, chrom, end-25, end+25)
        end_splits.sort()
        end_clusters = build_sr_clusters(end_splits)
        end_breaks = [build_breakends(c) for c in end_clusters]

        for be in end_breaks:
            for b in be:
                psl_rec = eval_break(os.path.dirname(args.blatref), b, 'left', chrom, start, end)

    # Locate TSD if possible:
    # may need to jiggle the start location for TSD search
    tries = [[0,0], [1,1], [-1,-1], [0,-1], [-1,0], [1,0], [0,1]]

    if psl_rec:
        max_tsd = -1
        best_tsd = []
        while len(tries) > 0:
            tsd_result = tsd(psl_rec, ref, b_left_init=tries[0][0], b_right_init=tries[0][1])
            if tsd_result[1] - tsd_result[0] > max_tsd:
                best_tsd = tsd_result
                max_tsd = tsd_result[1] - tsd_result[0]
            
            if len(tries) > 1:
                tries = tries[1:]

            else:
                tries = []


        if 0 in best_tsd:
            return ref_no_junc(args, chrom, start, end, orient, name)

        l_tsd_start, l_tsd_end, r_tsd_start, r_tsd_end, tsd_seq = best_tsd

        l_junc = junc_search(bams, chrom, start, 'L')
        r_junc = junc_search(bams, chrom, end, 'R')

        l_altcount, l_refcount = support(bams, chrom, l_junc, 'L', len(tsd_seq)) 
        r_altcount, r_refcount = support(bams, chrom, r_junc, 'R', len(tsd_seq))

        logger.info('%s:%d-%d  %s:%d-%d  %s  %d' % (chrom, start, end, chrom, l_tsd_end, r_tsd_start, psl_rec.cons, len(psl_rec.cons)))

        vaf = 0.0
        if l_altcount + r_altcount + l_refcount + r_refcount > 0:
            vaf = float(l_altcount + r_altcount) / float(l_altcount + r_altcount + l_refcount + r_refcount)

        if vaf == 0.0:
            return ref_no_junc(args, chrom, start, end, orient, name)

        vaf = str(vaf)

        junc_5p = l_junc
        junc_3p = r_junc

        tsd_start_5p = l_tsd_start
        tsd_end_5p   = l_tsd_end

        tsd_start_3p = r_tsd_start
        tsd_end_3p   = r_tsd_end

        altcount_5p = l_altcount
        altcount_3p = r_altcount

        refcount_5p = l_refcount
        refcount_3p = r_refcount

        if orient == '-':
            junc_5p, junc_3p = junc_3p, junc_5p
            tsd_start_5p, tsd_start_3p = tsd_start_3p, tsd_start_5p
            tsd_end_5p, tsd_end_3p = tsd_end_3p, tsd_end_5p
            altcount_5p, altcount_3p = altcount_3p, altcount_5p
            refcount_5p, refcount_3p = refcount_3p, refcount_5p
            

        if args.persample is not None:
            vaf = []

            altcount_5p = []
            altcount_3p = []
            refcount_5p = []
            refcount_3p = []

            with open(args.persample) as samples:
                for line in samples:
                    sbamfn, sname = line.strip().split()
                    sbam = [pysam.AlignmentFile(sbamfn)]

                    l_altcount, l_refcount = support(sbam, chrom, l_junc, 'L', len(tsd_seq)) 
                    r_altcount, r_refcount = support(sbam, chrom, r_junc, 'R', len(tsd_seq))

                    svaf = 0.0
                    if l_altcount + r_altcount + l_refcount + r_refcount > 0:
                        svaf = float(l_altcount + r_altcount) / float(l_altcount + r_altcount + l_refcount + r_refcount)

                    vaf.append('%s|%.3f' % (sname, svaf))

                    s_altcount_5p = l_altcount
                    s_altcount_3p = r_altcount
                    s_refcount_5p = l_refcount
                    s_refcount_3p = r_refcount

                    if orient == '-':
                        s_altcount_5p, s_altcount_3p = s_altcount_3p, s_altcount_5p
                        s_refcount_5p, s_refcount_3p = s_refcount_3p, s_refcount_5p

                    altcount_5p.append('%s|%d' % (sname, s_altcount_5p))
                    altcount_3p.append('%s|%d' % (sname, s_altcount_3p))
                    refcount_5p.append('%s|%d' % (sname, s_refcount_5p))
                    refcount_3p.append('%s|%d' % (sname, s_refcount_3p))

            vaf = ','.join(vaf)
            altcount_5p = ','.join(altcount_5p)
            altcount_3p = ','.join(altcount_3p)
            refcount_5p = ','.join(refcount_5p)
            refcount_3p = ','.join(refcount_3p)

        out = (chrom, start, end, orient, name, junc_5p, junc_3p, tsd_seq.upper(), psl_rec.cons, altcount_5p, altcount_3p, refcount_5p, refcount_3p, vaf)

        return out

    return ref_no_junc(args, chrom, start, end, orient, name)


def main(args):

    logger.info('%s started with cmd: %s' % (sys.argv[0], ' '.join(sys.argv)))

    p = start_blat_server(args.blatref)

    print('\t'.join(header))

    pool = mp.Pool(processes=int(args.procs))
    results = []

    with open(args.ins) as bed:
        for line in bed:
            chrom, start, end, name, orient = line.strip().split()[:5]
            start = int(start)
            end   = int(end)

            assert orient in ('+','-'), orient

            res = pool.apply_async(ref_ins, [args, chrom, start, end, orient, name])
            results.append(res)

    output = []

    for res in results:
        out = res.get()
        output.append(out)

    for out in output:
        print('\t'.join(map(str, out)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Find breakpoints, TSD, VAF, and read counts for reference insertions')
    parser.add_argument('-b', '--bam', required=True, help='initial BAM(s) for deletion discovery')
    parser.add_argument('-i', '--ins', required=True, help='insertion locations (five columns required: chrom, start, end, strand, annotation')
    parser.add_argument('-r', '--blatref', required=True, help='BLAT reference')
    parser.add_argument('-f', '--fastaref', required=True, help='samtools faidx indexed genome fasta')
    parser.add_argument('-p', '--procs', default=1, help='split work across multiple processes')
    parser.add_argument('--port', default=9999)
    parser.add_argument('--persample', default=None, help='List of files (2 column: BAM, Name) for per-sample information')
    args = parser.parse_args()
    main(args)

