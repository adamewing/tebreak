#!/usr/bin/env python
 
import os
import re
import sys
import time
import pysam
import align
import random
import argparse
import logging
import subprocess
import itertools
import traceback
import cPickle as pickle
import multiprocessing as mp

import numpy as np
import scipy.stats as ss
 
from uuid import uuid4
from string import maketrans
from operator import itemgetter
from collections import Counter
from collections import OrderedDict as od
from collections import defaultdict as dd
from bx.intervals.intersection import Intersecter, Interval # pip install bx-python

import profile
 
#######################################
## Classes                           ##
#######################################


class Genome:
    def __init__(self, gfn):
        ''' gfn = genome file name (.fai or chrom, length tsv) '''
        self.chrlen = {} # length of each chromosome
        self.chrmap = [] # used for picking chromosomes
 
        self.bp = 0
        bins = 100000
 
        with open(gfn, 'r') as g:
            for line in g:
                if not line.startswith('#'):
                    chrom, length = line.strip().split()[:2]
                    self.chrlen[chrom] = int(length)
                    self.bp += int(length)
     
        for chrom, length in self.chrlen.iteritems():
            self.chrmap += [chrom] * int(float(length) / float(self.bp) * bins)
 
    def pick(self):
        ''' return a random chromosome and position '''
        rchrom = random.choice(self.chrmap)
        rpos   = int(random.uniform(1, self.chrlen[rchrom]))
 
        return rchrom, rpos
 
    def addpad(self, interval, pad):
        ''' pad interval such that it doesn't go out of bounds '''
        chrom, start, end = interval
        start = int(start) - int(pad)
        end   = int(end) + int(pad)

        assert chrom in self.chrlen, "error padding interval %s, %s not a known chromosome" % (str(interval), chrom)

        if start < 0: start = 0
        if end > self.chrlen[chrom]: end = self.chrlen[chrom]

        return (chrom, start, end)


    def chunk(self, n, seed=None, sorted=False, pad=0):
        ''' break genome into n evenly-sized chunks, return n lists of (chrom, start, end) '''
        chunklen = int(self.bp/n)
        
        chunks = []
        intervals = []
 
        chunkleft = chunklen # track how much genome needs to go into each chunk
 
        chromlist = self.chrlen.keys()
 
        if sorted:
            chromlist.sort()
        else:
            if seed is not None: random.seed(seed)
            random.shuffle(chromlist)
 
        for chrom in chromlist:
            length = self.chrlen[chrom]
 
            lenleft = length
            if length <= chunkleft:
                chunkleft -= length
                lenleft -= length
                intervals.append( self.addpad((chrom, 0, length), pad) )
                assert lenleft == 0
 
                if chunkleft == 0:
                    chunkleft = chunklen
                    chunks.append(intervals)
                    intervals = []
            else:
                while lenleft > 0:
                    if lenleft >= chunkleft:
                        intervals.append( self.addpad((chrom, length-lenleft, length-lenleft+chunkleft), pad) )
                        lenleft -= chunkleft
 
                        chunkleft = chunklen
                        chunks.append(intervals)
                        intervals = []
 
                    else: # lenleft < chunkleft
                        intervals.append( self.addpad((chrom, length-lenleft, length), pad) )
                        chunkleft -= lenleft
                        lenleft -= lenleft
 
        return list(itertools.chain.from_iterable(chunks)) # flatten list


class LASTResult:
    def __init__(self, res):
        self.raw = res
        self.score = int(res[0].split()[1].replace('score=', ''))

        self.target_id      = res[1].split()[1]
        self.target_start   = int(res[1].split()[2])
        self.target_alnsize = int(res[1].split()[3])
        self.target_end     = self.target_start + self.target_alnsize
        self.target_strand  = res[1].split()[4]
        self.target_seqsize = int(res[1].split()[5])
        self.target_align   = res[1].split()[6]

        self.query_id = res[2].split()[1]
        self.query_distnum = 0
        if '|' in self.query_id:
            self.query_id = res[2].split()[1].split('|')[0]
            self.query_distnum = int(res[2].split()[1].split('|')[-1]) # track which distal read this came from

        self.query_start   = int(res[2].split()[2])
        self.query_alnsize = int(res[2].split()[3])
        self.query_end     = self.query_start + self.query_alnsize
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


class SortableRead:
    def __init__(self, read):
        self.read = read
        self.seq  = read.seq
        self.seqstart = read.reference_start-read.query_alignment_start

    def __gt__(self, other):
        if self.read.tid == other.read.tid:
            return self.seqstart > other.seqstart
        else:
            return self.read.tid > other.read.tid

 
class SplitRead:
    ''' store information about split read alignment '''
    def __init__(self, chrom, read, bamfn):
        self.uuid  = str(uuid4())
        self.chrom = chrom
        self.read  = read
        self.bamfn = os.path.basename(bamfn)

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
 
 
class DiscoRead:
    ''' store information about discordant pair alignment '''
    def __init__(self, chrom, read, bamfn, mate_chrom=None):
        self.chrom = chrom
        self.read  = read
        self.bamfn = os.path.basename(bamfn)

        self.mate_chrom = mate_chrom # can be None
        self.mate_read  = None  # set later
 
    def mate_mapped(self):
        return self.mate_read is not None and not self.mate_read.is_unmapped

    def getRG(self):
        ''' return read group from RG aux tag '''
        for tag, val in self.read.tags:
            if tag == 'RG': return val
        return None

    def __gt__(self, other):
        if self.mate_mapped() and other.mate_mapped():
            return self.mate_read.get_reference_positions()[0] > other.mate_read.get_reference_positions()[0]
        else:
            return self.read.get_reference_positions()[0] > other.read.get_reference_positions()[0]

    def __str__(self):
        return  ' '.join(map(str, (self.chrom, self.read, self.mate_chrom, self.mate_start)))


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
        self.median = self.reads[len(self)/2].breakpos
 
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

        sortable_reads = [SortableRead(sr.read) for sr in self.reads]
        seqs = [sorted_read.seq for sorted_read in sorted(sortable_reads)]

        cons = seqs[0]
        scores = []

        for seq in seqs[1:]:

            s1 = align.string_to_alignment(cons)
            s2 = align.string_to_alignment(seq)

            (s, a1, a2) = align.align(s1, s2, -2, -2, S, local=True)
            a1 = align.alignment_to_string(a1)
            a2 = ''.join([b for b in list(align.alignment_to_string(a2)) if b != '-'])

            score = float(len(a1) - (len(a1)-s)) / float(len(a1))
            scores.append(score)

            # print('%s\n%s\nScore: %f' % (a1, a2, score))

            if score < minscore and len(cons) == len(seq):
                cons = seq

            if score >= minscore:
                align_end = locate_subseq(seq, a2)[1]
                cons += seq[align_end:]

        return cons, np.mean(score)

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


class DiscoCluster(ReadCluster):
    ''' store and manipulate groups of DiscoRead objects '''

    def overlap_insertion(self, insertion):
        iv_ins = [insertion.min_supporting_base(), insertion.max_supporting_base()]
        iv_dsc = self.find_extrema()

        return min(iv_ins[1], iv_dsc[1]) - max(iv_ins[0], iv_dsc[0]) > 0

    def find_mate_extrema(self):
        ''' return leftmost and rightmost aligned positions in cluster vs. reference '''
        positions = []
        for r in self.reads:
            if r.mate_mapped():
                positions += [pos for pos in r.mate_read.positions]

        if len(positions) == 0: return -1, -1
        return min(positions), max(positions)

    def summary_tuple(self):
        return (self.reads[0].mate_chrom, self.find_mate_extrema()[0], self.find_mate_extrema()[1], self.readgroups(), self.bamfiles())

 
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
 
        self.mappings = []
 
    def microhomology(self):
        pass # placeholder
 
    def proximal_subread(self):
        ''' return mapping(s) containing breakpoint '''
        return [read for read in self.mappings if self.breakpos in read.get_reference_positions()]
 
    def distal_subread(self):
        ''' return mapping(s) not containing breakpoint '''
        return [read for read in self.mappings if not read.is_unmapped and self.breakpos not in read.get_reference_positions()]
 
    def unmapped_subread(self):
        ''' returns list of intervals and corresponding subseqs '''
        covered_bases = []
        for subseq in map(lambda x : x.query_alignment_sequence, self.mappings):
            subseq = orient_subseq(self.consensus, subseq)
            covered_bases += range(*locate_subseq(self.consensus, subseq))
 
        subseqs   = []
        intervals = []
 
        interval = []
        subseq   = []
 
        for p, b in enumerate(list(self.consensus)):
            if p not in covered_bases:
                if len(interval) > 0 and interval[-1]+1 < p:
                    intervals.append( (min(interval), max(interval)) )
                    subseqs.append(''.join(subseq))
                    interval = []
                    subseq   = []
 
                else:
                    interval.append(p)
                    subseq.append(b)
 
        if len(interval) > 0:
            intervals.append( (min(interval), max(interval)) )
            subseqs.append(''.join(subseq))
 
        return intervals, subseqs
 
    def __len__(self):
        return len(self.cluster)

    def __str__(self):
        return '%s:%d:%s:%s' % (self.chrom, self.breakpos, self.consensus, self.direction)

class Insertion:
    ''' store and compile information about an insertion with 1 or 2 breakpoints '''
    def __init__(self, be1=None, be2=None):
        self.uuid = str(uuid4())

        self.be1 = None
        self.be2 = None
 
        if be1 is not None:
            self.be1 = be1
        if be2 is not None:
            self.be2 = be2
 
        self.be1_alt = None
        self.be2_alt = None

        self.be1_improved_cons = False
        self.be2_improved_cons = False

        if self.paired():
            if self.be1.breakpos > self.be2.breakpos:
                self.be1, self.be2 = self.be2, self.be1 # keep breakends in position order

        self.info = od() # set with self.compile_info()
        self.discoreads = []
        #self.dr_clusters = []
        self.fastqrecs = []
        self.mappability = None

    def __len__(self):
        ''' return total number of reads (sr+dr) associated with insertion '''
        return self.num_sr() + len(self.discoreads)

    def num_sr(self):
        ''' return number of split reads supporting insertion '''
        l = 0
        if self.be1 is not None: l += len(self.be1.cluster)
        if self.be2 is not None: l += len(self.be2.cluster)
        return l

    def paired(self):
        return None not in (self.be1, self.be2)
 
    def breakend_overlap(self):
        if not self.paired(): return None
        if len(self.be1.proximal_subread()) == 0 or len(self.be2.proximal_subread()) == 0: return None
        return ref_dist(self.be1.proximal_subread()[0], self.be2.proximal_subread()[0])
 
    def min_supporting_base(self):
        ''' return leftmost supporting reference position covered '''
        sites = []
        for be in (self.be1, self.be2):
            if be is not None:
                if len(be.proximal_subread()) > 0:
                    for proxread in be.proximal_subread(): 
                        sites += proxread.get_reference_positions()

        if len(sites) == 0:
            return None

        return min(sites)

    def max_supporting_base(self):
        ''' return rightmost supporting reference position covered '''
        sites = []
        for be in (self.be1, self.be2):
            if be is not None:
                if len(be.proximal_subread()) > 0:
                    for proxread in be.proximal_subread(): 
                        sites += proxread.get_reference_positions()
                    
        if len(sites) == 0:
            return None

        return max(sites)

    def tsd(self, be1_use_prox=0, be2_use_prox=0):
        ''' target site duplication '''
        if not self.paired():
            return None
 
        if self.breakend_overlap() > 0:
            return None

        else:
            if len(self.be1.proximal_subread()) == 0 or len(self.be2.proximal_subread()) == 0: return None

            junc1 = self.be1.proximal_subread()[be1_use_prox]
            junc2 = self.be2.proximal_subread()[be2_use_prox]
 
            tsd_ref_interval = ref_overlap(junc1, junc2)

            if tsd_ref_interval is None: return None

            tsd_ref_interval[1] += 1
 
            tsdseq1 = ''
            tsdseq2 = ''
 
            #for (qrypos, refpos) in junc1.get_aligned_pairs(): # broken by pysam 8.3
            for qrypos, refpos in enumerate(junc1.get_reference_positions()):
                if refpos in range(*tsd_ref_interval):
                    if qrypos is not None:
                        tsdseq1 += junc1.seq[qrypos+junc1.qstart]
 
            #for (qrypos, refpos) in junc2.get_aligned_pairs(): # broken by pysam 8.3
            for qrypos, refpos in enumerate(junc2.get_reference_positions()):
                if refpos in range(*tsd_ref_interval):
                    if qrypos is not None:
                        tsdseq2 += junc2.seq[qrypos+junc2.qstart]

            return tsdseq1, tsdseq2

    def fetch_discordant_reads(self, bams, isize=10000):
        ''' Return list of DiscoRead objects '''
        chrom = self.be1.chrom
        start = self.min_supporting_base()
        end   = self.max_supporting_base()

        if None in (start, end): return []
     
        assert start < end, 'fetch_discordant_reads: start > end'

        mapped   = {}
        unmapped = {}

        for bam in bams:
            for read in bam.fetch(chrom, start, end):
                if read.is_paired and not read.is_unmapped and not read.is_duplicate:
                    chrom = str(bam.getrname(read.tid))
         
                    if read.mate_is_unmapped:
                        unmapped[read.qname] = DiscoRead(chrom, read, bam.filename)
         
                    else:
                        pair_dist = abs(read.reference_start - read.next_reference_start)
                        if read.tid != read.next_reference_id or pair_dist > isize:
                            mate_chrom = str(bam.getrname(read.next_reference_id))
                            mapped[read.qname] = DiscoRead(chrom, read, bam.filename, mate_chrom)

            # get mate info

            # mate mapped
            for qname, dr in mapped.iteritems():
                for read in bam.fetch(dr.mate_chrom, dr.read.next_reference_start, dr.read.next_reference_start+1):
                    if read.qname == qname and not read.is_secondary and not is_supplementary(read):
                        if read.seq != mapped[qname].read.seq:
                            mapped[qname].mate_read = read

            # mate unmapped
            for read in bam.fetch(chrom, start, end):
                if read.is_unmapped and read.qname in unmapped:
                    if not read.is_secondary and not is_supplementary(read):
                        unmapped[read.qname].mate_read = read

        self.discoreads = mapped.values() + unmapped.values()


    def align_filter(self, insref_fa, tmpdir='/tmp'):
        ''' check whether insertion consensus aligns to '''
        if not os.path.exists(insref_fa + '.tis'): build_last_db(insref_fa)

        out_fa = '%s/tebreak.align_filter.%s.fa' % (tmpdir, str(uuid4()))

        d1, d2, u1, u2 = (None, None, None, None)

        if self.be1 is not None:
            d1 = self.be1.distal_subread()       # pysam.AlignedSegment or None
            u1 = self.be1.unmapped_subread()[1]  # seq strings or None
            if d1 is not None: d1 = map(lambda x: x.seq, d1)

        if self.be2 is not None:
            d2 = self.be2.distal_subread()
            u2 = self.be2.unmapped_subread()[1]
            if d2 is not None: d2 = map(lambda x: x.seq, d2)

        seqsizes = [20]

        with open(out_fa, 'w') as fa:
            for seqlist in (d1, d2, u1, u2):
                if seqlist is not None:
                    for seq in seqlist:
                        if len(seq) >= 10:
                            seqsizes.append(len(seq))
                            fa.write('>%s\n%s\n' % (str(uuid4()), seq))

        la_results = align_last(out_fa, insref_fa, e=min(seqsizes)-2)

        os.remove(out_fa)

        if len(la_results) == 0: return True

        return False


    def improve_consensus(self, ctg_fa, bwaref, tmpdir='/tmp'):
        ''' attempt to use assembly of all supporting reads to build a better consensus '''
        cons_fasta = self.consensus_fasta(tmpdir)

        ctg_falib = load_falib(ctg_fa)

        build_last_db(ctg_fa)

        la_results = align_last(cons_fasta, ctg_fa, e=20)

        if self.be1 is not None:
            # find corresponding contig, if possible

            for res in la_results:
                if res.query_id == self.be1.uuid:
                    # criteria for deciding the new consensus is better
                    if res.target_seqsize > len(self.be1.consensus) and len(self.be1.consensus) - res.query_alnsize < 10 and res.pct_match() > 0.95:
                        self.be1_alt = self.be1
                        self.be1.consensus = ctg_falib[res.target_id]
                        self.be1_improved_cons = True


        if self.be2 is not None:
            # find corresponding contig, if possible

            for res in la_results:
                if res.query_id == self.be2.uuid:
                    # criteria for deciding the new consensus is better
                    if res.target_seqsize > len(self.be2.consensus) and len(self.be2.consensus) - res.query_alnsize < 10 and res.pct_match() > 0.95:
                        self.be2_alt = self.be2
                        self.be2.consensus = ctg_falib[res.target_id]
                        self.be2_improved_cons = True

        for ext in ('','.amb','.ann','.bck','.bwt','.des','.fai','.pac','.prj','.sa','.sds','.ssp','.suf','.tis'):
            if os.path.exists(ctg_fa+ext): os.remove(ctg_fa+ext)

        if os.path.exists(cons_fasta): os.remove(cons_fasta)

        return self.be1_improved_cons, self.be2_improved_cons


    def supportreads_fastq(self, outdir, min_readlen=50, limit=1000):
        ''' discordant support reads marked DR, split support reads marked SR '''
        assert os.path.exists(outdir)

        outreads  = od()
        usedreads = {}

        out_fastq = outdir + '/' + '.'.join(('supportreads', self.be1.chrom, str(self.be1.breakpos), 'fq'))
        with open(out_fastq, 'w') as out:
            for readstore in (self.be1, self.be2, self.discoreads):
                if readstore:
                    try:
                        rtype = 'SR'
                        readlist = readstore.cluster.reads
                    except:
                        rtype = 'DR'
                        readlist = readstore

                    for r in readlist:
                        read = r.read
                        name = read.qname
                        unseen = True
                        if read.is_read1:
                            if name + '/1' in usedreads: unseen = False
                            usedreads[name + '/1'] = True
                            name += '.%s/1' % rtype

                        if read.is_read2:
                            if name + '/2' in usedreads: unseen = False
                            usedreads[name + '/2'] = True
                            name += '.%s/2' % rtype

                        if len(read.seq) > min_readlen and unseen: outreads[name] = read.seq + '\n+\n' + read.qual

                        if rtype == 'DR' and r.mate_read is not None: # get discordant mates
                            read = r.mate_read
                            unseen = True

                            if read.is_read1:
                                if name + '/1' in usedreads: unseen = False
                                usedreads[name + '/1'] = True
                                name += '.%s/1' % rtype

                            if read.is_read2:
                                if name + '/2' in usedreads: unseen = False
                                usedreads[name + '/2'] = True
                                name += '.%s/2' % rtype

                            if len(read.seq) > min_readlen and unseen: outreads[name] = read.seq + '\n+\n' + read.qual

            if len(outreads) >= limit or len(outreads) == 0: return None

            for name, data in outreads.iteritems():
                out.write('@%s\n%s\n' % (name, data))
                self.fastqrecs.append('@%s\n%s\n' % (name, data))

        return out_fastq

    def consensus_fasta(self, tmpdir='/tmp'):
        assert os.path.exists(tmpdir)

        out_fasta = tmpdir + '/' + '.'.join(('consensus', self.be1.chrom, str(self.be1.breakpos), str(uuid4()), 'fa'))
        with open(out_fasta, 'w') as out:
            out.write('>%s\n%s\n' % (self.be1.uuid, self.be1.consensus))
            if self.be2 is not None:
                out.write('>%s\n%s\n' % (self.be2.uuid, self.be2.consensus))

        return out_fasta

    def compile_info(self, bams):
        ''' fill self.info with summary info, needs original bam for chromosome lookup '''
        if self.be1 == None and self.be2 == None:
            return None

        self.info['ins_uuid'] = self.uuid
        self.info['chrom'] = self.be1.chrom
        self.info['min_supporting_base'] = self.min_supporting_base()
        self.info['max_supporting_base'] = self.max_supporting_base()
        self.info['mappability'] = self.mappability

        self.info['be1_breakpos'] = self.be1.breakpos
        self.info['be1_obj_uuid'] = self.be1.uuid
        
        #seqs
        self.info['be1_cons_seq'] = self.be1.consensus
        self.info['be1_prox_seq'] = ','.join(map(lambda x : x.query_alignment_sequence, self.be1.proximal_subread()))
        self.info['be1_dist_seq'] = ','.join(map(lambda x : x.query_alignment_sequence, self.be1.distal_subread()))
        self.info['be1_umap_seq'] = ','.join(self.be1.unmapped_subread()[1])

        if self.be1.proximal_subread() and self.be1.proximal_subread()[0].is_reverse:
            self.info['be1_prox_str'] = '-'
        else:
            self.info['be1_prox_str'] = '+'

        self.info['be1_prox_loc'] = []
        for subseq in self.info['be1_prox_seq'].split(','):
            self.info['be1_prox_loc'].append(locate_subseq(self.be1.consensus, orient_subseq(self.be1.consensus, subseq)))

        # stats
        self.info['be1_sr_count'] = len(self.be1)
        self.info['be1_num_maps'] = len(self.be1.mappings)
        self.info['be1_cons_scr'] = self.be1.consscore
        self.info['be1_median_D'] = self.be1.cluster.median_D()
        self.info['be1_avgmatch'] = self.be1.cluster.avg_matchpct()
        self.info['be1_rg_count'] = self.be1.cluster.readgroups()
        self.info['be1_bf_count'] = self.be1.cluster.bamfiles()
        self.info['be1_prox_mpq'] = ','.join(map(lambda x : str(x.mapq), self.be1.proximal_subread()))
        self.info['be1_improved'] = self.be1_improved_cons

        if self.info['be1_dist_seq'] == '':
            self.info['be1_dist_seq'] = None
        else:
            self.info['be1_dist_chr'] = ','.join(map(lambda x : bams[0].getrname(x.tid), self.be1.distal_subread()))
            self.info['be1_dist_pos'] = ','.join(map(lambda x : str(x.get_reference_positions()[0]), self.be1.distal_subread()))
            self.info['be1_dist_end'] = ','.join(map(lambda x : str(x.get_reference_positions()[-1]), self.be1.distal_subread()))
            self.info['be1_dist_mpq'] = ','.join(map(lambda x : str(x.mapq), self.be1.distal_subread()))

        if self.info['be1_umap_seq'] == '':
            self.info['be1_umap_seq'] = None

        if self.be2 is not None:
            self.info['be2_breakpos'] = self.be2.breakpos
            self.info['be2_obj_uuid'] = self.be2.uuid
            self.info['be2_cons_seq'] = self.be2.consensus
            self.info['be2_prox_seq'] = ','.join(map(lambda x : x.query_alignment_sequence, self.be2.proximal_subread()))
            self.info['be2_dist_seq'] = ','.join(map(lambda x : x.query_alignment_sequence, self.be2.distal_subread()))
            self.info['be2_umap_seq'] = ','.join(self.be2.unmapped_subread()[1])

            if self.be2.proximal_subread() and self.be2.proximal_subread()[0].is_reverse:
                self.info['be2_prox_str'] = '-'
            else:
                self.info['be2_prox_str'] = '+'


            self.info['be2_prox_loc'] = []
            for subseq in self.info['be2_prox_seq'].split(','):
                self.info['be2_prox_loc'].append(locate_subseq(self.be2.consensus, orient_subseq(self.be2.consensus, subseq)))

            # stats
            self.info['be2_sr_count'] = len(self.be2)
            self.info['be2_num_maps'] = len(self.be2.mappings)
            self.info['be2_cons_scr'] = self.be2.consscore
            self.info['be2_median_D'] = self.be2.cluster.median_D()
            self.info['be2_avgmatch'] = self.be2.cluster.avg_matchpct()
            self.info['be2_rg_count'] = self.be2.cluster.readgroups()
            self.info['be2_bf_count'] = self.be2.cluster.bamfiles()
            self.info['be2_prox_mpq'] = ','.join(map(lambda x : str(x.mapq), self.be2.proximal_subread()))
            self.info['be2_improved'] = self.be2_improved_cons

            if self.info['be2_dist_seq'] == '':
                self.info['be2_dist_seq'] = None
            else:
                self.info['be2_dist_chr'] = ','.join(map(lambda x: bams[0].getrname(x.tid), self.be2.distal_subread()))
                self.info['be2_dist_pos'] = ','.join(map(lambda x: str(x.get_reference_positions()[0]), self.be2.distal_subread()))
                self.info['be2_dist_end'] = ','.join(map(lambda x: str(x.get_reference_positions()[-1]), self.be2.distal_subread()))
                self.info['be2_dist_mpq'] = ','.join(map(lambda x : str(x.mapq), self.be2.distal_subread()))                

            if self.info['be2_umap_seq'] == '':
                self.info['be2_umap_seq'] = None

            self.info['be1_use_prox'] = 0
            self.info['be2_use_prox'] = 0

            if self.info['be1_prox_loc'] == self.info['be2_prox_loc']: # insertion may be completely assembled
                if len(self.info['be1_prox_loc']) > 1:
                    self.info['be1_use_prox'] = 1

                elif len(self.info['be2_prox_loc']) > 1:
                    self.info['be2_use_prox'] = 1

            tsdpair = self.tsd(be1_use_prox=self.info['be1_use_prox'], be2_use_prox=self.info['be2_use_prox'])
            if tsdpair is not None:
                self.info['be1_end_over'], self.info['be2_end_over'] = tsdpair

        if 'be2_breakpos' not in self.info:
            self.info['be2_breakpos'] = self.info['be1_breakpos']

        if 'be2_sr_count' not in self.info:
            self.info['be2_sr_count'] = 0

        self.info['dr_count'] = len(self.discoreads)
        self.info['dr_unmapped_mates'] = len([dr for dr in self.discoreads if dr.mate_read is not None and dr.mate_read.is_unmapped])

 
#######################################
## Functions                         ##
#######################################
 
 
def rc(dna):
    ''' reverse complement '''
    complements = maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
    return dna.translate(complements)[::-1]
 

def read_matchpct(read):
    ''' return number of mismatches / aligned length of read '''
    nm = [value for (tag, value) in read.tags if tag == 'NM'][0]
    return 1.0 - (float(nm)/float(read.alen))

 
def ref_overlap(read1, read2):
    ''' return overlapping interval in ref. coords (not chrom), none otherwise '''
 
    if read1 is None or read2 is None:
        return None
 
    if read1.is_unmapped or read2.is_unmapped:
        return None

    iv1 = sorted((read1.get_reference_positions()[0], read1.get_reference_positions()[-1]))
    iv2 = sorted((read2.get_reference_positions()[0], read2.get_reference_positions()[-1]))
 
    if min(iv1[1], iv2[1]) - max(iv1[0], iv2[0]) > 0: # is there overlap?
        return [max(iv1[0], iv2[0]), min(iv1[1], iv2[1])]
 
    return None
 
 
def ref_dist(read1, read2):
    ''' return distance between intervals in ref., overlapping = negative values '''
 
    if read1 is None or read2 is None:
        return None
 
    iv1 = sorted((read1.get_reference_positions()[0], read1.get_reference_positions()[-1]))
    iv2 = sorted((read2.get_reference_positions()[0], read2.get_reference_positions()[-1]))
 
    return max(iv1[0], iv2[0]) - min(iv1[1], iv2[1])
 
 
def orient_subseq(longseq, shortseq):
    ''' return shortseq in same orientation as longseq '''
    assert len(longseq) >= len(shortseq), 'orient_subseq: %s < %s' % (longseq, shortseq)
 
    if re.search(shortseq, longseq):
        return shortseq
    else:
        assert re.search(rc(shortseq), longseq), "orient_subseq: %s not a subseq of %s" %(shortseq, longseq)
        return rc(shortseq)
 
 
def locate_subseq(longseq, shortseq):
    ''' return (start, end) of shortseq in longseq '''
    assert len(longseq) >= len(shortseq), 'orient_subseq: %s < %s' % (longseq, shortseq)
 
    match = re.search(shortseq, longseq)
    if match is not None:
        return sorted((match.start(0), match.end(0)))
 
    return None
 
 
def is_primary(read):
    ''' not included in pysam '''
    return not (read.is_secondary or is_supplementary(read))
 
 
def is_supplementary(read):
    ''' pysam does not currently include a check for this flag '''
    return bin(read.flag & 2048) == bin(2048)

 
def fetch_clipped_reads(bams, chrom, start, end, filters):
    ''' Return list of SplitRead objects '''
    assert filters['min_minclip'] > 2
 
    splitreads = []
 
    start = int(start)
    end   = int(end)
 
    assert start < end
 
    for bam in bams:
        for read in bam.fetch(chrom, start, end):
            masked = False
            if filters['genome_mask'] is not None and chrom in filters['genome_mask']:
                if filters['genome_mask'][chrom].find(read.pos, read.pos+1): masked = True

            if not masked and not read.is_unmapped and not read.is_duplicate: #and read.mapq > 0:
                if read.rlen - read.alen >= int(filters['min_minclip']): # 'soft' clipped?
     
                    # length of 'minor' clip
                    altclip = min(read.qstart, read.rlen-read.qend)

                    # junk bases
                    N_count = 0
                    if 'N' in read.seq: N_count = Counter(read.seq)['N']
     
                    if altclip <= 2: # could add as a filter
                        if N_count <= filters['max_N_consensus'] and splitqual(read) <= filters['max_D_score']:
                            chrom = str(bam.getrname(read.tid))
                            if len(read.get_reference_positions()) > 0:
                                splitreads.append(SplitRead(chrom, read, bam.filename))
 
    return splitreads
 
 
def splitqual(read):
    ''' return KS-test D value for clipped vs. unclipped bases'''
    
    breakpos = None
 
    breakpos = read.get_aligned_pairs()[-1][0] # breakpoint on right
 
    q1 = map(ord, list(read.qual[:breakpos]))
    q2 = map(ord, list(read.qual[breakpos:]))
 
    return ss.ks_2samp(q1, q2)[0]
 

def load_falib(infa):
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


def build_breakends(cluster, filters, tmpdir='/tmp'):
    ''' returns list of breakends from cluster '''
    breakends = []

    for breakpos in cluster.all_breakpoints():
        for dir in ('left', 'right'):
            subcluster = cluster.subcluster_by_breakend([breakpos], direction=dir)
            if len(subcluster) >= filters['min_sr_per_break'] and subcluster.max_cliplen() >= filters['min_maxclip']:
                seq     = subcluster.reads[0].read.seq
                score   = 1.0

                if len(subcluster) > 1: seq, score = subcluster.consensus()
 
                N_count = 0
                if 'N' in seq: N_count = Counter(seq)['N']

                if seq != '' and score >= filters['min_consensus_score'] and N_count <= filters['max_N_consensus']:
                    breakends.append(BreakEnd(cluster.chrom, breakpos, subcluster, seq, score, dir))
 
    return breakends


def map_breakends(breakends, db, tmpdir='/tmp'):
    ''' remap consensus sequences stored in BreakEnd objects '''
    tmp_fa = tmpdir + '/' + '.'.join(('tebreak', str(uuid4()), 'be.fa'))
    breakdict = {} # for faster lookup
 
    with open(tmp_fa, 'w') as out:
        for be in breakends:
            be.mappings = []
            breakdict[be.uuid] = be
            qual = 'I' * len(be.consensus)
            out.write('>%s\n%s\n+\n%s\n' % (be.uuid, be.consensus, qual))
 
    tmp_sam = '.'.join(tmp_fa.split('.')[:-1]) + '.sam'

    FNULL = open(os.devnull, 'w')

    with open(tmp_sam, 'w') as out:
        sam_cmd  = ['bwa', 'mem', '-k', '10', '-w', '500', '-M', '-v', '0', db, tmp_fa]
        p = subprocess.Popen(sam_cmd, stdout=subprocess.PIPE, stderr=FNULL)

        for line in p.stdout:
            out.write(line)
 
    sam = pysam.AlignmentFile(tmp_sam)
 
    passed_parse = False

    # rarely, pysam will fail to parse the bwa-mem generated SAM and I haven't worked out why... workaround for now
    while not passed_parse:
        try:
            for i, read in enumerate(sam.fetch(until_eof=True)):
                breakdict[read.qname].mappings.append(read)
            passed_parse = True

        except IOError as e:
            sys.stderr.write("warning: pysam failed parse at read %d, modifying file and re-try...\n" % i)
            sam.close()
            lines = []
            with open(tmp_sam, 'r') as sam:
                lines = [line for line in sam]

            with open(tmp_sam, 'w') as sam:
                for n, line in enumerate(lines):
                    if line.startswith('@'): n -= 1 # header
                    if n < i: sam.write(line)

            sam = pysam.AlignmentFile(tmp_sam)

    sam.close()

    os.remove(tmp_fa)
    os.remove(tmp_sam)

    return breakdict.values()


def build_last_db(fa):
    ''' make db for LAST alignments '''
    subprocess.call(['lastdb', '-s', '4G', fa, fa])
    assert os.path.exists(fa + '.tis'), 'could not lastdb -4G %s %s' % (fa, fa)


def align_last(fa, db, e=20):
    ''' returns list of LASTResult objects '''
    assert os.path.exists(db + '.tis'), 'not a lastdb: %s' % db

    last_cmd = ['lastal', '-e', str(e), db, fa]

    la_lines   = []
    la_results = []

    p = subprocess.Popen(last_cmd, stdout=subprocess.PIPE)

    for line in p.stdout:
        if not line.startswith('#'):
            if line.strip() != '':
                la_lines.append(line.strip())

            else:
                la_results.append(LASTResult(la_lines))
                la_lines = []

    return la_results


def score_breakend_pair(be1, be2, k=2.5, s=3.0):
    ''' assign a score to a breakend, higher is "better" '''
    prox1 = be1.proximal_subread()
    prox2 = be2.proximal_subread()

    if prox1 and prox2:
        prox1 = prox1[0]
        prox2 = prox2[0]
        overlap = abs(min(0, ref_dist(prox1, prox2))) # overlap = negative distance between proximal read mappings i.e. potential TSD
        weighted_overlap = ss.gamma(k, scale=s).pdf(overlap) * float(overlap)**2 # TSD length distribution taken into account
        distance_penalty = 0

        if overlap > 0:  distance_penalty = abs(abs(be1.breakpos-be2.breakpos) - overlap) # disagreement in TSD length
        if overlap == 0: distance_penalty = abs(be1.breakpos-be2.breakpos) # no TSD

        score = weighted_overlap - distance_penalty + len(be1) + len(be2)
        return score
 
    return None


def checkref(ref_fasta):
    assert os.path.exists(ref_fasta), 'reference not found: %s' % ref_fasta
    assert os.path.exists(ref_fasta + '.fai'), 'please run samtools faidx on %s' % ref_fasta
    assert os.path.exists(ref_fasta + '.bwt'), 'please run bwa index on %s' % ref_fasta


def build_insertions(breakends, maxdist=100):
    ''' return list of Insertion objects '''
    insertions = []
    be_dict = dict([(be.uuid, be) for be in breakends])

    be_itree = Intersecter() # interval tree

    for be in breakends:
        be_itree.add_interval(Interval(be.breakpos-maxdist, be.breakpos+maxdist, value=be.uuid))

    pair_scores = []

    checked_pairs = {}

    for be1 in breakends:
        for be2_coords in be_itree.find(be1.breakpos, be1.breakpos+1):
            be2 = be_dict[be2_coords.value]

            pair_name = '-'.join(sorted((be1.uuid, be2.uuid)))

            if be1.uuid != be2.uuid and be1.direction != be2.direction and pair_name not in checked_pairs:
                if abs(be1.breakpos-be2.breakpos) <= maxdist:
                    pair_scores.append((be1.uuid, be2.uuid, score_breakend_pair(be1, be2)))

            checked_pairs[pair_name] = True

    # sort breakends by score, descending
    pair_scores = [score for score in pair_scores if score[2] is not None]
    pair_scores.sort(key=itemgetter(2),reverse=True)

    used = {} # each breakend can only be used once
    for be1_uuid, be2_uuid, score in pair_scores:
        if be1_uuid not in used and be2_uuid not in used:
            insertions.append(Insertion(be_dict[be1_uuid], be_dict[be2_uuid]))
            #print "debug:, insertion: %s and %s" % (be_dict[be1_uuid], be_dict[be2_uuid])
        used[be1_uuid] = True
        used[be2_uuid] = True

    # single-end detections
    for be in breakends:
        if be.uuid not in used:
            insertions.append(Insertion(be))
            used[be.uuid] = True
     
    return insertions


def minia(fq, tmpdir='/tmp'):
    ''' sequence assembly '''
    # minia temp files don't seem to be compatabile with concurrency, workaround w/ temp cwd
    oldcwd = os.getcwd()
    tmpcwd = '%s/%s' % (tmpdir, 'tebreak.'+str(uuid4()))
    os.mkdir(tmpcwd)
    assert os.path.exists(tmpcwd), 'cannot create temp dir: %s' % tmpcwd

    os.chdir(tmpcwd)

    if not os.path.exists(fq):
        fq = oldcwd + '/' + fq

    ctgbase = tmpdir + '/tebreak.minia.%s' % str(uuid4())
    
    cmd = ['minia', '-in', fq, '-abundance-min', '1', '-no-length-cutoff', '-verbose', '0', '-nb-cores', '1', '-out', ctgbase]

    FNULL = open(os.devnull, 'w')
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=FNULL)
    for line in p.stdout: pass

    if os.path.exists(ctgbase + '.h5'): os.remove(ctgbase + '.h5')

    os.chdir(oldcwd)
    os.rmdir(tmpcwd)

    return ctgbase + '.contigs.fa'


def build_mask(bedfile):
    ''' return a dictionary of interval trees '''
    forest = dd(Intersecter)

    with open(bedfile, 'r') as bed:
        for line in bed:
            chrom, start, end = line.strip().split()[:3]
            start = int(start)
            end   = int(end)

            forest[chrom].add_interval(Interval(start, end))

    return forest


def avgmap(maptabix, chrom, start, end):
    ''' return average mappability across chrom:start-end region; maptabix = pysam.Tabixfile'''
    scores = []

    if None in (start, end): return None

    if chrom in maptabix.contigs:
        for rec in maptabix.fetch(chrom, int(start), int(end)):
            mchrom, mstart, mend, mscore = rec.strip().split()
            mstart, mend = int(mstart), int(mend)
            mscore = float(mscore)

            while mstart < mend and mstart:
                mstart += 1
                if mstart >= int(start) and mstart <= int(end):
                    scores.append(mscore)

        if len(scores) > 0:
            return sum(scores) / float(len(scores))
        else:
            return 0.0
    else:
        return 0.0


def summarise_insertion(ins):
    ''' returns a pickleable version of the insertion information '''
    pi = dd(dict)

    pi['INFO'] = ins.info
    pi['READSTORE'] = ins.fastqrecs

    return pi


def filter_insertions(insertions, filters, tmpdir='/tmp'):
    filtered = []
    for ins in insertions:
        exclude = False

        mapq = map(lambda x : str(x.mapq), ins.be1.proximal_subread())
        rgs  = [rg.split('|')[0] for rg in ins.be1.cluster.readgroups()]
        bams = [bam.split('|')[0] for bam in ins.be1.cluster.bamfiles()]

        if ins.be2 is not None and len(ins.be2.proximal_subread()) > 0:
            mapq += map(lambda x : str(x.mapq), ins.be2.proximal_subread())
            rgs  += [rg.split('|')[0] for rg in ins.be2.cluster.readgroups()]
            bams += [bam.split('|')[0] for bam in ins.be2.cluster.bamfiles()]

        bams = list(set(bams)) # uniqify

        if len(ins) >= filters['max_ins_reads']: exclude = True
        if max(mapq) < filters['min_prox_mapq']: exclude = True

        if len(ins.discoreads) < filters['min_discordant_reads']: exclude = True
        if ins.num_sr() < filters['min_split_reads']: exclude = True

        if filters['exclude_bam']:
            for bam in bams:
                if bam in filters['exclude_bam']: exclude = True

        if filters['exclude_readgroup']:
            for rg in rgs:
                if rg in filters['exclude_rg']: exclude = True

        if filters['max_bam_count'] > 0:
            if len(bams) > filters['max_bam_count']: exclude = True

        if filters['insertion_library'] is not None and not exclude:
            if ins.align_filter(filters['insertion_library'], tmpdir=tmpdir): exclude = True

        if filters['map_tabix'] is not None and not exclude:
            if ins.be1.chrom in filters['map_tabix'].contigs:
                ins.mappability = avgmap(filters['map_tabix'], ins.be1.chrom, ins.min_supporting_base(), ins.max_supporting_base())
            else:
                ins.mappability = 0.0

            if ins.mappability < filters['min_mappability']: exclude = True

        if not exclude: filtered.append(ins)

    return filtered


def postprocess_insertions(insertions, filters, bwaref, bams, tmpdir='/tmp'):
    for ins in insertions:
        support_fq  = ins.supportreads_fastq(tmpdir, limit=filters['max_ins_reads'])
        if support_fq is None: return insertions

        cwd = os.getcwd()
        support_asm = minia(support_fq, tmpdir=tmpdir)

        retry_counter = 0 # minia might not be the most reliable option...
        while not os.path.exists(support_asm) and retry_counter < 10:
            retry_counter += 1
            sys.stderr.write('***Assembly retry: %s:%d\n' % (ins.be1.chrom, ins.be1.breakpos))
            support_asm = minia(support_fq, tmpdir=tmpdir)

        if not os.path.exists(support_asm):
            sys.stderr.write('***Assembly failed!: %s:%d\n' % (ins.be1.chrom, ins.be1.breakpos))

        else:
            #sys.stderr.write('Assembled: %s:%d, filename: %s\n' % (ins.be1.chrom, ins.be1.breakpos, support_asm))
            ins.improve_consensus(support_asm, bwaref, tmpdir=tmpdir)

        if os.path.exists(support_fq): os.remove(support_fq)
        if os.path.exists(support_asm): os.remove(support_asm)

        if os.getcwd() != cwd: os.chdir(cwd)

    # collect altered breakends
    alt_be_list = []
    for ins in insertions:
        if ins.be1_improved_cons: alt_be_list.append(ins.be1)
        if ins.be2_improved_cons: alt_be_list.append(ins.be2)

    remap_be_dict = {}
    for be in map_breakends(alt_be_list, bwaref, tmpdir=tmpdir):
        remap_be_dict[be.uuid] = be

    for ins in insertions:
        # use previous mapping if new consensus did not map
        if ins.be1_improved_cons:
            if ins.be1.uuid in remap_be_dict:
                if len(remap_be_dict[ins.be1.uuid].proximal_subread()) > 0:
                    ins.be1 = remap_be_dict[ins.be1.uuid]
                else:
                    ins.be1 = ins.be1_alt
                    ins.be1_improved_cons = False

        if ins.be2_improved_cons:
            if ins.be2.uuid in remap_be_dict:
                if len(remap_be_dict[ins.be2.uuid].proximal_subread()) > 0:
                    ins.be2 = remap_be_dict[ins.be2.uuid]
                else:
                    ins.be2 = ins.be2_alt
                    ins.be2_improved_cons = False

        if ins.be1_improved_cons or ins.be2_improved_cons:
            ins.compile_info(bams)

    return insertions



def run_chunk(args, exp_rpkm, chrom, start, end):
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    chunkname = '%s:%d-%d' % (chrom, start, end)

    try:
        bams = [pysam.AlignmentFile(bam, 'rb') for bam in args.bam.split(',')]

        # would do this outside but can't pass a non-pickleable object
        if args.mask is not None: args.mask = build_mask(args.mask)

        if args.map_tabix is not None: args.map_tabix = pysam.Tabixfile(args.map_tabix)

        start = int(start)
        end   = int(end)

        filters = {
            'min_maxclip':           int(args.min_maxclip),
            'min_minclip':           int(args.min_minclip),
            'min_sr_per_break':      int(args.min_sr_per_break),
            'min_consensus_score':   int(args.min_consensus_score),
            'max_D_score':           float(args.maxD),
            'max_ins_reads':         int(args.max_ins_reads),
            'min_split_reads':       int(args.min_split_reads),
            'min_discordant_reads':  int(args.min_discordant_reads),
            'min_prox_mapq':         int(args.min_prox_mapq),
            'max_N_consensus':       int(args.max_N_consensus),
            'max_rpkm':              int(args.max_fold_rpkm)*exp_rpkm,
            'exclude_bam':           [],
            'exclude_readgroup':     [],
            'max_bam_count':         int(args.max_bam_count),
            'insertion_library':     args.insertion_library,
            'genome_mask':           args.mask,
            'map_tabix':             args.map_tabix,
            'min_mappability':       float(args.min_mappability)
        }

        if args.exclude_bam is not None: filters['exclude_bam'] = map(os.path.basename, args.exclude_bam.split(','))
        if args.exclude_readgroup is not None: filters['exclude_readgroup'] = args.exclude_readgroup.split(',')

        insertions = []
     
        logger.debug('Processing chunk: %s ...' % chunkname)
        logger.debug('Chunk %s: Parsing split reads from bam(s): %s ...' % (chunkname, args.bam))
        sr = fetch_clipped_reads(bams, chrom, start, end, filters)
        sr.sort()


        logger.debug('Chunk %s: Building clusters from %d split reads ...' % (chunkname, len(sr)))
        clusters = build_sr_clusters(sr)
     
        logger.debug('Chunk %s: Building breakends...' % chunkname)

        breakends = []

        for cluster in clusters:
            cl_readcount = 0
            cl_min, cl_max = cluster.find_extrema()

            for bam in bams:
                cl_readcount += sum([not read.is_unmapped for read in bam.fetch(cluster.chrom, cl_min, cl_max)])

            rpkm = cl_readcount/((cl_max-cl_min)/1000.)

            if filters['max_rpkm'] == 0 or rpkm < filters['max_rpkm']:
                breakends += build_breakends(cluster, filters, tmpdir=args.tmpdir)

            else:
                logger.debug('Chunk %s, cluster %d-%d over max RPKM with %f' % (chunkname, cl_min, cl_max, rpkm))
     
        logger.debug('Chunk %s: Mapping %d breakends ...' % (chunkname, len(breakends)))
        if len(breakends) > 0:
            breakends = map_breakends(breakends, args.bwaref, tmpdir=args.tmpdir)
         
            logger.debug('Chunk %s: Building insertions...' % chunkname)

            insertions = build_insertions(breakends)
            logger.debug('Chunk %s: Processing and filtering %d potential insertions ...' % (chunkname, len(insertions)))

            insertions = [ins for ins in insertions if len(ins.be1.proximal_subread()) > 0] # remove bogus insertions

            for ins in insertions:
                ins.fetch_discordant_reads(bams)
                ins.compile_info(bams)

            insertions = filter_insertions(insertions, filters, tmpdir=args.tmpdir)

            logger.debug('Chunk %s: Postprocessing %d filtered insertions, trying to improve consensus breakend sequences ...' % (chunkname, len(insertions)))
            processed_insertions  = postprocess_insertions(insertions, filters, args.bwaref, bams, tmpdir=args.tmpdir)

            logger.debug('Chunk %s: Summarising insertions ...' % chunkname)
            summarised_insertions = [summarise_insertion(ins) for ins in processed_insertions]

            logger.debug('Finished chunk: %s' % chunkname)

            for bam in bams: bam.close()
            return summarised_insertions

        else:
            for bam in bams: bam.close()
            return []

    except Exception, e:
        sys.stderr.write('*'*60 + '\tencountered error in chunk: %s\n' % chunkname)
        traceback.print_exc(file=sys.stderr)
        sys.stderr.write("*"*60 + "\n")

        return []


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


def text_summary(insertions, outfile='tebreak.out'):
    with open(outfile, 'w') as out:
        for ins in insertions:
            out.write('#BEGIN\n')
            for label, value in ins['INFO'].iteritems():
                out.write('%s: %s\n' % (label, str(value)))
            out.write('#END\n')

            out.write('\n')


def expected_rpkm(bam_files, genome, intervals=None):
    ''' expected reads per kilobase mapped '''
    bams = [pysam.AlignmentFile(bam_file, 'rb') for bam_file in bam_files]
    total_mapped_reads = sum([bam.mapped for bam in bams])
    km = genome.bp/1000.

    if intervals is not None:
        total_length = 0
        total_mapped_reads = 0

        with open(intervals, 'r') as bed:
            for line in bed:
                chrom, start, end = line.strip().split()
                start = int(start)
                end   = int(end)

                total_length += end - start

                for bam in bams:
                    total_mapped_reads += sum([not read.is_unmapped for read in bam.fetch(chrom, start, end)])

        km = total_length/1000.

    for bam in bams: bam.close()

    return total_mapped_reads/km


def main(args):
    ''' housekeeping '''
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)

    checkref(args.bwaref)

    if not args.no_shared_mem:
        sys.stderr.write("loading bwa index %s into shared memory ...\n" % args.bwaref)
        p = subprocess.Popen(['bwa', 'shm', args.bwaref], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        for line in p.stdout: pass # wait for bwa to load

    ''' Chunk genome or use input BED '''
    
    procs = int(args.processes)
    chunk_count = int(args.chunks)

    if chunk_count < procs: chunks = procs

    pool = mp.Pool(processes=procs)

    genome = Genome(args.bwaref + '.fai')

    chunks = []
    exp_rpkm = 0

    if args.interval_bed is None or args.wg_rpkm:
        chunks = genome.chunk(chunk_count, sorted=True, pad=5000)

        if args.rpkm_bam:
            exp_rpkm = expected_rpkm(args.rpkm_bam.split(','), genome)
        else:
            exp_rpkm = expected_rpkm(args.bam.split(','), genome)

    else:
        if args.rpkm_bam:
            exp_rpkm = expected_rpkm(args.rpkm_bam.split(','), genome, intervals=args.interval_bed)
        else:
            exp_rpkm = expected_rpkm(args.bam.split(','), genome, intervals=args.interval_bed)


    if args.interval_bed is not None:
        with open(args.interval_bed, 'r') as bed:
            chunks = [(line.strip().split()[0], int(line.strip().split()[1]), int(line.strip().split()[2])) for line in bed]


    if exp_rpkm < 10:
        sys.stderr.write("expected RPKM is less than 10, ignoring high RPKM cutoffs...\n")
        exp_rpkm = 0

    reslist = []

    for chunk in chunks:
        # run_chunk(args, exp_rpkm, chunk[0], chunk[1], chunk[2]) # uncomment for mp debug
        res = pool.apply_async(run_chunk, [args, exp_rpkm, chunk[0], chunk[1], chunk[2]])
        reslist.append(res)

    insertions = []
    for res in reslist:
        insertions += res.get()
    
    insertions = resolve_duplicates(insertions)
    text_summary(insertions, outfile=args.detail_out)

    pickoutfn = re.sub('.bam$', '.tebreak.pickle', os.path.basename(args.bam))
    if args.pickle is not None: pickoutfn = args.pickle

    with open(pickoutfn, 'w') as pickout:
        pickle.dump(insertions, pickout)

    #if not args.no_shared_mem:
    #    sys.stderr.write("unloading bwa index %s from shared memory ...\n" % args.bwaref)
    #    p = subprocess.Popen(['bwa', 'shm', '-d'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #    for line in p.stdout: pass # wait for bwa to unload

    logger.debug('Pickled to %s' % pickoutfn)

 
if __name__ == '__main__':
    # set up logger
    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(format=FORMAT)
 
    parser = argparse.ArgumentParser(description='Find inserted sequences vs. reference')
    parser.add_argument('-b', '--bam', required=True, help='target BAM(s): can be comma-delimited list')
    parser.add_argument('-r', '--bwaref', required=True, help='bwa/samtools indexed reference genome')
    parser.add_argument('-p', '--processes', default=1, help='split work across multiple processes')
    parser.add_argument('-c', '--chunks', default=1, help='split genome into chunks (default = # processes), helps control memory usage')
    parser.add_argument('-i', '--interval_bed', default=None, help='BED file with intervals to scan')
    parser.add_argument('-D', '--maxD', default=0.8, help='maximum value of KS D statistic for split qualities (default = 0.8)')
    parser.add_argument('--min_minclip', default=3, help='min. shortest clipped bases per cluster (default = 3)')
    parser.add_argument('--min_maxclip', default=10, help='min. longest clipped bases per cluster (default = 10)')
    parser.add_argument('--min_sr_per_break', default=1, help='minimum split reads per breakend (default = 1)')
    parser.add_argument('--min_consensus_score', default=0.95, help='quality of consensus alignment (default = 0.95)')
    parser.add_argument('-m', '--mask', default=None, help='BED file of masked regions')

    parser.add_argument('--rpkm_bam', default=None, help='use alternate BAM(s) for RPKM calculation: use original BAMs if using reduced BAM(s) for -b/--bam')
    parser.add_argument('--max_fold_rpkm', default=10, help='ignore insertions supported by rpkm*max_fold_rpkm reads (default = 10)')
    parser.add_argument('--max_ins_reads', default=1000, help='maximum number of reads per insertion call (default = 1000)')
    parser.add_argument('--min_split_reads', default=4, help='minimum total split reads per insertion call (default = 4)')
    parser.add_argument('--min_discordant_reads', default=4, help='minimum discordant read count (default = 4)')
    parser.add_argument('--min_prox_mapq', default=10, help='minimum map quality for proximal subread (default = 10)')
    parser.add_argument('--max_N_consensus', default=4, help='exclude breakend seqs with > this number of N bases (default = 4)')
    parser.add_argument('--exclude_bam', default=None, help='may be comma delimited')
    parser.add_argument('--exclude_readgroup', default=None, help='may be comma delimited')
    parser.add_argument('--max_bam_count', default=0, help='maximum number of bams supporting per insertion')
    parser.add_argument('--insertion_library', default=None, help='for pre-selecting insertion types')
    parser.add_argument('--map_tabix', default=None, help='tabix-indexed BED of mappability scores')
    parser.add_argument('--min_mappability', default=0.5, help='minimum mappability (default = 0.5; only matters with --map_tabix)')

    parser.add_argument('--tmpdir', default='/tmp', help='temporary directory (default = /tmp)')
    parser.add_argument('--pickle', default=None, help='pickle output name')
    parser.add_argument('--detail_out', default='tebreak.out', help='file to write detailed output')
 
    parser.add_argument('--wg_rpkm', default=False, action='store_true', help='force calculate rpkm over whole genome')
    parser.add_argument('--no_shared_mem', default=False, action='store_true')
 
    args = parser.parse_args()
    main(args)
