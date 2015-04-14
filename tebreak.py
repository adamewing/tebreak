#!/usr/bin/env python
 
import os
import re
import sys
import pysam
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
from collections import Counter
from collections import OrderedDict as od
from collections import defaultdict as dd


 
#######################################
## Classes                           ##
#######################################
 
 
class AlignedColumn:
    ''' used by MSA class to store aligned bases '''
    def __init__(self):
        self.bases = od() # species name --> base
        self.annotations = od() # metadata
 
    def gap(self):
        if '-' in self.bases.values(): return True
        return False
 
    def subst(self):
        if self.gap(): return False
        if len(set(self.bases.values())) > 1: return True
        return False
 
    def cons(self):
        ''' consensus prefers bases over gaps '''
        for base, count in Counter(map(str.upper, self.bases.values())).most_common():
            if base != '-': return base
 
    def score(self):
        topbase = self.cons()
        nongaps = [b for b in map(str.upper, self.bases.values()) if b != '-']
        matches = [b for b in map(str.upper, self.bases.values()) if b == topbase]

        if len(nongaps) == 0: return 0.0
 
        return float(len(matches)) / float(len(nongaps))
 
 
    def __str__(self):
        return str(self.bases)
 
 
class MSA:
    ''' multiple sequence alignment class '''
    def __init__(self, infile=None):
        self.columns = []
        self.ids     = []
        self.seqs    = od()
 
        if infile is not None: self.readFastaMSA(infile)
 
    def __len__(self):
        return len(self.columns)
 
    def readFastaMSA(self, infile):
        id   = None
        seq  = ''
 
        with open(infile, 'r') as fasta:
            for line in fasta:
                line = line.strip()
                if line.startswith('>'):
                    if id is not None:
                        self.seqs[id] = seq
                    seq = ''
                    id = line.lstrip('>')
                    self.ids.append(id)
                else:
                    seq += line
            self.seqs[id] = seq
 
        first = True
        colen = 0
        for ID, seq in self.seqs.iteritems():
            if first:
                colen = len(seq)
                for base in list(seq):
                    ac = AlignedColumn()
                    ac.bases[ID] = base
                    self.columns.append(ac)
                first = False
            else:
                assert len(seq) == colen
                pos = 0
                for base in list(seq):
                    ac = self.columns[pos]
                    ac.bases[ID] = base
                    pos += 1
 
    def consensus(self):
        ''' compute consensus '''
        bases  = [column.cons() for column in self.columns]
        scores = [column.score() for column in self.columns]
 
        if bases is not None and None not in bases:
            return ''.join(bases), np.mean(scores)
        else:
            return '', np.mean(scores)


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

    def make_fasta(self, outdir):
        ''' for downstream consensus building '''
        out_fasta = outdir + '/' + '.'.join(('tebreak', str(uuid4()), 'fasta'))
        with open(out_fasta, 'w') as out:
            for sr in self.reads:
                read = sr.read
                name = read.qname
                if read.is_read1:
                    name += '/1'
                if read.is_read2:
                    name += '/2'
 
                out.write('>%s\n%s\n' % (name, read.seq))
 
        return out_fasta
 
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
    def __init__(self, chrom, breakpos, cluster, consensus, score):
        self.uuid      = str(uuid4())
        self.cluster   = cluster
        self.chrom     = chrom
        self.breakpos  = breakpos
        self.consensus = consensus
        self.consscore = score
 
        self.mappings = []
 
    def microhomology(self):
        pass # placeholder
 
    def proximal_subread(self):
        ''' return mapping(s) containing breakpoint '''
        return [read for read in self.mappings if self.breakpos in read.get_reference_positions()]
 
    def distal_subread(self):
        ''' return mapping(s) not containing breakpoint '''
        return [read for read in self.mappings if self.breakpos not in read.get_reference_positions()]
 
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
        return '%s\t%d' % (self.chrom, self.breakpos)


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
        self.dr_clusters = []
        self.fastqrecs = []
        self.cons_fasta = None

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
        be2 = []
        if self.be2 is not None:
            be2 = self.be2.proximal_subread()[0].get_reference_positions()

        return min(self.be1.proximal_subread()[0].get_reference_positions() + be2)

    def max_supporting_base(self):
        ''' return rightmost supporting reference position covered '''
        be2 = []
        if self.be2 is not None:
            be2 = self.be2.proximal_subread()[0].get_reference_positions()

        return max(self.be1.proximal_subread()[0].get_reference_positions() + be2)

    def tsd(self, be1_use_prox=0, be2_use_prox=0):
        ''' target site duplication '''
        if not self.paired():
            return None
 
        if self.breakend_overlap() > 0:
            return None

        else:
            if len(self.be1.distal_subread()) == 0 or len(self.be2.distal_subread()) == 0: return None
            if len(self.be1.proximal_subread()) == 0 or len(self.be2.proximal_subread()) == 0: return None

            junc1 = self.be1.proximal_subread()[be1_use_prox]
            junc2 = self.be2.proximal_subread()[be2_use_prox]
 
            tsd_ref_interval = ref_overlap(junc1, junc2)

            if tsd_ref_interval is None: return None

            tsd_ref_interval[1] += 1
 
            tsdseq1 = ''
            tsdseq2 = ''
 
            for (qrypos, refpos) in junc1.get_aligned_pairs():
                if refpos in range(*tsd_ref_interval):
                    if qrypos is not None:
                        tsdseq1 += junc1.seq[qrypos+junc1.qstart]
 
            for (qrypos, refpos) in junc2.get_aligned_pairs():
                if refpos in range(*tsd_ref_interval):
                    if qrypos is not None:
                        tsdseq2 += junc2.seq[qrypos+junc2.qstart]
 
            return tsdseq1, tsdseq2

    def fetch_discordant_reads(self, bams, isize=10000):
        ''' Return list of DiscoRead objects '''
        chrom = self.be1.chrom
        start = self.min_supporting_base()
        end   = self.max_supporting_base()
     
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


    def improve_consensus(self, ctg_fa, bwaref):
        ''' attempt to use assembly of all supporting reads to build a better consensus '''
        assert self.cons_fasta is not None, "need to write consensus fasta before trying to improve consensus"

        ctg_falib = load_falib(ctg_fa)

        subprocess.call(['lastdb', '-s', '4G', ctg_fa, ctg_fa])
        assert os.path.exists(ctg_fa + '.tis'), 'could not lastdb -4G %s %s' % (ctg_fa, ctg_fa)

        last_cmd = ['lastal', '-e', '20', ctg_fa, self.cons_fasta]

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

        if self.be1 is not None:
            # find corresponding contig, if possible
            be_name = 'tebreak:%s:%d' % (self.be1.chrom, self.be1.breakpos)

            for res in la_results:
                if res.query_id == be_name:
                    print "LASTres:", res
                    # criteria for deciding the new consensus is better
                    if res.target_seqsize > len(self.be1.consensus) and len(self.be1.consensus) - res.query_alnsize < 10 and res.pct_match() > 0.95:
                        self.be1_alt = self.be1
                        self.be1.consensus = ctg_falib[res.target_id]
                        self.be1_improved_cons = True


        if self.be2 is not None:
            # find corresponding contig, if possible
            be_name = 'tebreak:%s:%d' % (self.be2.chrom, self.be2.breakpos)

            for res in la_results:
                if res.query_id == be_name:
                    # criteria for deciding the new consensus is better
                    if res.target_seqsize > len(self.be2.consensus) and len(self.be2.consensus) - res.query_alnsize < 10 and res.pct_match() > 0.95:
                        self.be2_alt = self.be2
                        self.be2.consensus = ctg_falib[res.target_id]
                        self.be2_improved_cons = True

        for ext in ('','.amb','.ann','.bck','.bwt','.des','.fai','.pac','.prj','.sa','.sds','.ssp','.suf','.tis'):
            if os.path.exists(ctg_fa+ext): os.remove(ctg_fa+ext)


        return self.be1_improved_cons, self.be2_improved_cons


    def unmapped_fastq(self, outdir):
        ''' for downstream analysis of unmapped paired reads '''
        assert os.path.exists(outdir), 'cannot find path: %s' % outpath

        out_fastq = outdir + '/' + '.'.join(('disc_unmap', self.be1.chrom, str(self.be1.breakpos), 'fq'))
        with open(out_fastq, 'w') as out:
            for dr in self.discoreads:
                read = dr.mate_read
                if read is not None and read.is_unmapped:
                    name = read.qname
                    if read.is_read1:
                        name += '/1'
                    if read.is_read2:
                        name += '/2'
 
                    out.write('@%s\n%s\n+\n%s\n' % (name, read.seq, read.qual))
 
        return out_fastq

    def supportreads_fastq(self, outdir):
        ''' discordant support reads marked DR, split support reads marked SR '''
        assert os.path.exists(outdir)

        outreads = od()

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
                        if read.is_read1:
                            name += '.%s/1' % rtype
                        if read.is_read2:
                            name += '.%s/2' % rtype

                        outreads[name] = read.seq + '\n+\n' + read.qual

                        if rtype == 'DR' and r.mate_read is not None: # get discordant mates
                            read = r.mate_read
                            if read.is_read1:
                                name += '.%s/1' % rtype
                            if read.is_read2:
                                name += '.%s/2' % rtype

                            outreads[name] = read.seq + '\n+\n' + read.qual

            for name, data in outreads.iteritems():
                out.write('@%s\n%s\n' % (name, data))
                self.fastqrecs.append('@%s\n%s\n' % (name, data))

        return out_fastq

    def consensus_fasta(self, outdir):
        assert os.path.exists(outdir)

        out_fasta = outdir + '/' + '.'.join(('consensus', self.be1.chrom, str(self.be1.breakpos), 'fa'))
        with open(out_fasta, 'w') as out:
            out.write('>tebreak:%s:%d\n%s\n' % (self.be1.chrom, self.be1.breakpos, self.be1.consensus))
            if self.be2 is not None:
                out.write('>tebreak:%s:%d\n%s\n' % (self.be2.chrom, self.be2.breakpos, self.be2.consensus))

        self.cons_fasta = out_fasta

        return out_fasta

    def compile_info(self, bams):
        ''' fill self.info with summary info, needs original bam for chromosome lookup '''
        if self.be1 == None and self.be2 == None:
            return None

        self.info['ins_uuid'] = self.uuid
        self.info['chrom'] = self.be1.chrom

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

        # discordant read info
        self.dr_clusters = build_dr_clusters(self)

        self.info['dr_count'] = len(self.discoreads)
        self.info['dr_clusters']  = map(lambda x : x.summary_tuple(), self.dr_clusters)
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
 
 
def fetch_clipped_reads(bams, chrom, start, end, minclip=3, maxaltclip=2, maxD=0.8):
    ''' Return list of SplitRead objects '''
    assert minclip > maxaltclip
 
    splitreads = []
 
    start = int(start)
    end   = int(end)
 
    assert start < end
 
    for bam in bams:
        for read in bam.fetch(chrom, start, end):
            if not read.is_unmapped and not read.is_duplicate:
     
                if read.rlen - read.alen >= int(minclip): # 'soft' clipped?
     
                    # length of 'minor' clip
                    altclip = min(read.qstart, read.rlen-read.qend)
     
                    if altclip <= maxaltclip:
                        if splitqual(read) <= maxD:
                            chrom = str(bam.getrname(read.tid))
                            splitreads.append(SplitRead(chrom, read, bam.filename))
 
    return splitreads
 
 
def splitqual(read):
    ''' return KS-test D value for clipped vs. unclipped bases'''
    
    breakpos = None
 
    breakpos = read.get_aligned_pairs()[-1][0] # breakpoint on right
 
    q1 = map(ord, list(read.qual[:breakpos]))
    q2 = map(ord, list(read.qual[breakpos:]))
 
    return ss.ks_2samp(q1, q2)[0]
 
 
def mafft(infafn, iterations=100, threads=1, tmpdir='/tmp'):
    ''' use MAFFT to create MSA '''
 
    outfafn = tmpdir + '/' + str(uuid4()) + '.aln.fa'
 
    args = ['mafft', '--localpair', '--maxiterate', str(iterations), '--thread', str(threads), infafn]
 
    FNULL = open(os.devnull, 'w')
 
    with open(outfafn, 'w') as outfa:
        p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=FNULL)
        for line in p.stdout:
            outfa.write(line)
 
    return outfafn
 

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
 

def build_sr_clusters(splitreads, searchdist=100): # TODO PARAM
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


def build_dr_clusters(insertion, searchdist=100): # TODO PARAM
    ''' cluster discordant read ends assocaited with insertion '''

    clusters  = []
 
    for dr in sorted(insertion.discoreads):
        if len(clusters) == 0:
            clusters.append(DiscoCluster(dr))
 
        elif clusters[-1].chrom != dr.mate_chrom:
            clusters.append(DiscoCluster(dr))
 
        else:
            if ref_overlap(clusters[-1].reads[-1].mate_read, dr.mate_read) is None:
                clusters.append(DiscoCluster(dr))
 
            else:
                clusters[-1].add_read(dr)

    return clusters
 
 
def build_breakends(cluster, minsr_break, mincs, min_maxclip=10, tmpdir='/tmp'):
    ''' returns list of breakends from cluster '''
    breakends = []
 
    for breakpos in cluster.all_breakpoints():
        subclusters = (cluster.subcluster_by_breakend([breakpos], direction='left'), cluster.subcluster_by_breakend([breakpos], direction='right'))
        for subcluster in subclusters:
            if len(subcluster) >= minsr_break:
                seq     = subcluster.reads[0].read.seq
                score   = 1.0
                maxclip = subcluster.max_cliplen()

                if len(subcluster) > 1:
                    out_fa     = subcluster.make_fasta(tmpdir)
                    align_fa   = mafft(out_fa, tmpdir=tmpdir)
                    msa        = MSA(align_fa)
                    seq, score = msa.consensus()

                    os.remove(out_fa)
                    os.remove(align_fa)

                #print '*'*60
                #print subcluster
                #print seq
                #print score, maxclip
 
                if seq != '' and score >= mincs and min_maxclip <= maxclip:
                    breakends.append(BreakEnd(cluster.chrom, breakpos, subcluster, seq, score))
 
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
 
    with open(tmp_sam, 'w') as out:
        sam_cmd  = ['bwa', 'mem', '-k', '10', '-w', '500', '-M', '-v', '0', db, tmp_fa]
        p = subprocess.Popen(sam_cmd, stdout=subprocess.PIPE)
        for line in p.stdout:
            out.write(line)
 
    sam = pysam.AlignmentFile(tmp_sam)
 
    for read in sam.fetch(until_eof=True):
        breakdict[read.qname].mappings.append(read)
 
    os.remove(tmp_fa)
    os.remove(tmp_sam)
 
    return breakdict.values()
 
 
def score_breakend_pair(be1, be2):
    ''' assign a score to a breakend, higher is "better" '''
    prox1 = be1.proximal_subread()[0]
    prox2 = be2.proximal_subread()[0]

    score =  (1.0 - ref_dist(prox1, prox2)) + np.sqrt(len(be1.cluster) + len(be2.cluster))
 
    if prox1 and prox2:
        # prefer overlapping breakend reads (TSDs), but small TSDs can be overridden with high enough read counts
        print "prox1, prox2: (%d-%d) vs (%d-%d) " % (prox1.get_reference_positions()[0], prox1.get_reference_positions()[-1], prox2.get_reference_positions()[0], prox2.get_reference_positions()[-1])
        print "score       :",score
        return score 
 
    return 0
 
 
def checkref(ref_fasta):
    assert os.path.exists(ref_fasta), 'reference not found: %s' % ref_fasta
    assert os.path.exists(ref_fasta + '.fai'), 'please run samtools faidx on %s' % ref_fasta
    assert os.path.exists(ref_fasta + '.bwt'), 'please run bwa index on %s' % ref_fasta
 
 
def build_insertions(breakends, maxdist=100):
    ''' return list of Insertion objects '''
    insertions = []
    paired = {}
    for be1 in breakends:
        best_match = None
        best_pair  = None
 
        for be2 in breakends:
            if be1.proximal_subread() and be2.proximal_subread():
                dist = ref_dist(be1.proximal_subread()[0], be2.proximal_subread()[0])
 
                if be1.uuid != be2.uuid and dist <= maxdist and be2.uuid not in paired and be1.uuid not in paired:
                    if best_match is None:
                        best_match = score_breakend_pair(be1, be2)
                        best_pair  = (be1, be2)
 
                    else:
                        if score_breakend_pair(be1, be2) > best_match:
                            best_match = score_breakend_pair(be1, be2)
                            best_pair  = (be1, be2)
 
        if best_pair is not None:
            paired[best_pair[0].uuid] = True
            paired[best_pair[1].uuid] = True
            insertions.append(Insertion(*best_pair))
 
        else:
            if be1.uuid not in paired:
                paired[be1.uuid] = True
                insertions.append(Insertion(be1))
     
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
    cmd = ['minia', '-in', fq, '-abundance-min', '1', '-no-length-cutoff', '-verbose', '0', '-out', ctgbase]
    subprocess.call(cmd)
    os.remove(ctgbase + '.h5')

    os.chdir(oldcwd)
    os.rmdir(tmpcwd)

    return ctgbase + '.contigs.fa'


def summarise_insertion(ins):
    ''' returns a pickleable version of the insertion information '''
    pi = dd(dict)

    pi['INFO'] = ins.info
    pi['READSTORE'] = ins.fastqrecs

    return pi


def filter_insertions(insertions, filters):
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

        if len(ins) >= filters['max_ins_reads']: exclude = True
        if max(mapq) < filters['min_prox_mapq']: exclude = True

        if len(ins.discoreads) < filters['min_discordant_reads']: exclude = True
        if ins.num_sr() < filters['min_split_reads']: exclude = True

        if filters['restrict_to_bam'] is not None:
            if len(bams) > 1 or filters['restrict_to_bam'] not in bams: exclude = True

        if filters['restrict_to_readgroup'] is not None:
            if len(rgs) > 1 or filters['restrict_to_readgroup'] not in rgs: exclude = True

        if not exclude: filtered.append(ins)

    return filtered


def postprocess_insertions(insertions, bwaref, outpath, bams, tmpdir='/tmp'):
    for ins in insertions:
        support_fq  = ins.supportreads_fastq(outpath)
        support_asm = minia(support_fq, tmpdir=tmpdir)

        retry_counter = 0 # minia might not be the most reliable option...
        while not os.path.exists(support_asm) and retry_counter < 10:
            retry_counter += 1
            sys.stderr.write('***Assembly retry: %s:%d\n' % (ins.be1.chrom, ins.be1.breakpos))
            support_asm = minia(support_fq, tmpdir=tmpdir)

        if not os.path.exists(support_asm):
            sys.stderr.write('***Assembly failed!: %s:%d\n' % (ins.be1.chrom, ins.be1.breakpos))

        else:
            sys.stderr.write('Assembled: %s:%d, filename: %s\n' % (ins.be1.chrom, ins.be1.breakpos, support_asm))
            ins.improve_consensus(support_asm, bwaref)

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
            ins.consensus_fasta(outpath)

    return insertions



def run_chunk(args, chrom, start, end):
    logger = logging.getLogger(__name__)
    if args.verbose: logger.setLevel(logging.DEBUG)

    try:
        bams = [pysam.AlignmentFile(bam, 'rb') for bam in args.bam.split(',')]
        minsr_break = int(args.min_sr_per_break)
        mincs = float(args.min_consensus_score)
        maxD  = float(args.maxD)

        start = int(start)
        end   = int(end)

        filters = {
            'max_ins_reads':         int(args.max_ins_reads),
            'min_split_reads':       int(args.min_split_reads),
            'min_discordant_reads':  int(args.min_discordant_reads),
            'min_prox_mapq':         int(args.min_prox_mapq),
            'restrict_to_bam':       args.restrict_to_bam,
            'restrict_to_readgroup': args.restrict_to_readgroup
        }

        insertions = []
     
        chunkname = '%s:%d-%d' % (chrom, start, end)

        logger.debug('Processing chunk: %s ...' % chunkname)
        logger.debug('Chunk %s: Parsing split reads from bam(s): %s ...' % (chunkname, args.bam))
        sr = fetch_clipped_reads(bams, chrom, start, end, minclip=int(args.min_minclip), maxD=maxD)
        sr.sort()

        logger.debug('Chunk %s: Building clusters from %d split reads ...' % (chunkname, len(sr)))
        clusters = build_sr_clusters(sr)
     
        breakends = []
        for cluster in clusters:
            print cluster
            breakends += build_breakends(cluster, minsr_break, mincs, min_maxclip=int(args.min_maxclip), tmpdir=args.tmpdir)
     
        logger.debug('Chunk %s: Mapping %d breakends ...' % (chunkname, len(breakends)))
        if len(breakends) > 0:
            breakends = map_breakends(breakends, args.bwaref, tmpdir=args.tmpdir)
         
            insertions = build_insertions(breakends)
            logger.debug('Chunk %s: Processing %d insertions ...' % (chunkname, len(insertions)))

            insertions = [ins for ins in insertions if len(ins.be1.proximal_subread()) > 0] # remove bogus insertions

            for ins in insertions:
                ins.fetch_discordant_reads(bams)
                ins.compile_info(bams)

            insertions = filter_insertions(insertions, filters)

            for ins in insertions:
                # various FASTA/FASTQ outputs
                ins.unmapped_fastq(args.fasta_out_path)
                ins.consensus_fasta(args.fasta_out_path)

            logger.debug('Chunk %s: Postprocessing, trying to improve consensus breakend sequences ...' % chunkname)
            processed_insertions  = postprocess_insertions(insertions, args.bwaref, args.fasta_out_path, bams, tmpdir=args.tmpdir)

            logger.debug('Chunk %s: Summarising insertions ...' % chunkname)
            summarised_insertions = [summarise_insertion(ins) for ins in processed_insertions]

            logger.debug('Finished chunk: %s' % chunkname)

            return summarised_insertions

        else:
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


def text_summary(insertions):
    for ins in insertions:
        print '#BEGIN'
        for label, value in ins['INFO'].iteritems():
            print '%s: %s' % (label, str(value))
        print '#END'

        # some test filters...
        bestmatch = ins['INFO']['be1_avgmatch']
        if 'be2_avgmatch' in ins['INFO']:
            bestmatch = max(ins['INFO']['be1_avgmatch'], ins['INFO']['be2_avgmatch'])

        print "\n"

 
def main(args):
    ''' housekeeping '''
    logger = logging.getLogger(__name__)
    if args.verbose: logger.setLevel(logging.DEBUG)

    checkref(args.bwaref)

    if not os.path.exists(args.fasta_out_path):
        os.mkdir(args.fasta_out_path)
        assert os.path.exists(args.fasta_out_path), "could not create directory: %s" % args.fasta_out_path

    ''' Chunk genome or use input BED '''
    
    procs = int(args.processes)
    pool = mp.Pool(processes=procs)
    chunks = []

    if args.interval_bed is None:
        genome = Genome(args.bwaref + '.fai')
        chunks = genome.chunk(procs, sorted=True, pad=5000)

    else:
        with open(args.interval_bed, 'r') as bed:
            chunks = [(line.strip().split()[0], int(line.strip().split()[1]), int(line.strip().split()[2])) for line in bed]

    reslist = []

    for chunk in chunks:
        # ins_info = run_chunk(args, *chunk) # uncomment for mp debug
        # text_summary(ins_info)             # uncomment for mp debug
        res = pool.apply_async(run_chunk, [args, chunk[0], chunk[1], chunk[2]])
        reslist.append(res)

    insertions = []
    for res in reslist:
        insertions += res.get()
    
    insertions = resolve_duplicates(insertions)
    text_summary(insertions)

    pickoutfn = re.sub('.bam$', '.tebreak.pickle', os.path.basename(args.bam))

    with open(pickoutfn, 'w') as pickout:
        pickle.dump(insertions, pickout)

    logger.debug('Pickled to %s' % pickoutfn)

 
if __name__ == '__main__':
    # set up logger
    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(format=FORMAT)
 
    parser = argparse.ArgumentParser(description='Find inserted sequences vs. reference')
    parser.add_argument('-b', '--bam', required=True, help='target BAM(s): can be comma-delimited list')
    parser.add_argument('-r', '--bwaref', required=True, help='bwa/samtools indexed reference genome')
    parser.add_argument('-p', '--processes', default=1, help='split work across multiple processes')
    parser.add_argument('-i', '--interval_bed', default=None, help='BED file with intervals to scan')
    parser.add_argument('-D', '--maxD', default=0.8, help='maximum value of KS D statistic for split qualities (default = 0.8)')
    parser.add_argument('--min_minclip', default=3, help='min. shortest clipped bases per cluster (default = 3)')
    parser.add_argument('--min_maxclip', default=10, help='min. longest clipped bases per cluster (default = 10)')
    parser.add_argument('--min_sr_per_break', default=1, help='minimum split reads per breakend (default = 1)')
    parser.add_argument('--min_consensus_score', default=0.95, help='quality of consensus alignment (default = 0.95)')

    parser.add_argument('--max_ins_reads', default=1000)
    parser.add_argument('--min_split_reads', default=4)
    parser.add_argument('--min_discordant_reads', default=4)
    parser.add_argument('--min_prox_mapq', default=10)
    parser.add_argument('--restrict_to_bam', default=None)
    parser.add_argument('--restrict_to_readgroup', default=None)

    parser.add_argument('--tmpdir', default='/tmp', help='temporary directory (default = /tmp)')
    parser.add_argument('-o', '--fasta_out_path', default='tebreak_seqdata', help='path for FASTA output')
 
    parser.add_argument('-v', '--verbose', action='store_true')
 
    args = parser.parse_args()
    main(args)
