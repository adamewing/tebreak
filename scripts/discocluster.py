#!/usr/bin/env python

import os
import pysam
import random
import logging
import argparse
import itertools
logger = logging.getLogger(__name__)

import multiprocessing as mp

from collections import defaultdict as dd
from bx.intervals.intersection import Intersecter, Interval


''' identify clusters of discordant read ends where one end is in BED file '''


class Genome:
    def __init__(self, bamfn):
        bam = pysam.AlignmentFile(bamfn, 'rb')
        self.chrlen = {r: l for r,l in zip(bam.references, bam.lengths)}
        self.bp = sum(bam.lengths)


    def addpad(self, interval, pad):
        ''' pad interval such that it doesn't go out of bounds '''
        chrom, start, end = interval
        start = int(start) - int(pad)
        end   = int(end) + int(pad)

        assert chrom in self.chrlen, "error padding interval %s, %s not a known chromosome" % (str(interval), chrom)

        if start < 0: start = 0
        if end > self.chrlen[chrom]: end = self.chrlen[chrom]

        return (chrom, start, end)


    def chunk(self, n, seed=None, sorted=True, pad=0):
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
 
        return chunks


class Coord:
    def __init__(self, chrom, start, end, strand, mchrom, mstart, mend, mstrand, label, bam_name):
        self.chrom   = chrom
        self.start   = int(start)
        self.end     = int(end)
        self.strand  = strand
        self.mchrom  = mchrom
        self.mstart  = int(mstart)
        self.mend    = int(mend)
        self.mstrand = mstrand
        self.label   = label
        self.bam     = bam_name

        # if strand of genome element is '-', flip apparent mate strand
        elt_str = self.label.split('|')[-1]
        assert elt_str in ('+', '-'), 'malformed input BED: last three cols need to be class, family, orientation (+/-)'

        if elt_str == '-': self.mstrand = flip(self.mstrand)


    def __gt__(self, other):
        if self.chrom == other.chrom:
            return self.start > other.start
        else:
            return self.chrom > other.chrom


    def __str__(self):
        return '\t'.join(map(str, (self.bam, self.label, self.chrom, self.start, self.end, self.strand, self.mchrom, self.mstart, self.mend, self.mstrand)))


class InsCall:
    def __init__(self, coord_list, chrom, start, end, strand, bamlist, mapscore, nonref):
        self.coord_list = coord_list
        self.chrom      = chrom
        self.start      = int(start)
        self.end        = int(end)
        self.strand     = strand
        self.bamlist    = bamlist
        self.length     = len(coord_list)
        self.mapscore   = mapscore
        self.nonref     = nonref


    def out(self, verbose=True):
        output = ['#BEGIN']
        output.append('%s\t%d\t%d\t%s\t%s\t%d\t%0.3f\t%s' % (self.chrom, self.start, self.end, self.strand, self.bamlist, self.length, self.mapscore, self.nonref))
        if verbose:
            for c in self.coord_list: output.append(str(c))
        output.append('#END')

        return '\n'.join(output)


    def overlaps(self, other):
        ''' return true if overlap > 0 '''
        return min(self.end, other.end) - max(self.start, other.start) > 0


    def __gt__(self, other):
        if self.chrom == other.chrom:
            return self.start > other.start
        else:
            return self.chrom > other.chrom


    def __str__(self):
        return self.out(verbose=False)


def flip(strand):
    if strand == '+': return '-'
    if strand == '-': return '+'


def avgmap(maptabix, chrom, start, end):
    ''' return average mappability across chrom:start-end region; maptabix = pysam.Tabixfile '''
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


def interval_forest(bed_file):
    ''' build dictionary of interval trees '''
    forest = dd(Intersecter)

    with open(bed_file, 'r') as bed:
        for line in bed:
            chrom, start, end = line.strip().split()[:3]
            label = '|'.join(line.strip().split())
            forest[chrom].add_interval(Interval(int(start), int(end), value=label))

    return forest


def read_gen(bam, chrom=None, start=None, end=None):
    if None in (chrom, start, end):
        for read in bam.fetch():
            yield read

    else:
        for read in bam.fetch(chrom, start, end):
            yield read



def get_coords(forest, bams, chrom=None, start=None, end=None, min_mapq=1, min_dist=10000):
    coords = []

    for bam in bams:
        tick = 10000000
        try:
            tick = int((bam.mapped + bam.unmapped) * 0.01)
            if tick == 0: tick = 1
            logger.debug('outputting status every %d reads (1 pct)' % tick)

        except ValueError as e:
            logger.debug('no index found, outputting status every %d reads' % tick)

        #for i, read in enumerate(bam.fetch()):
        for i, read in enumerate(read_gen(bam, chrom=chrom, start=start, end=end)):
            if not read.is_unmapped and not read.mate_is_unmapped and not read.is_duplicate:

                rchrom = bam.getrname(read.reference_id)
                rstart = read.reference_start
                rend   = read.reference_end

                rstr = '+'
                if read.is_reverse: rstr = '-'

                mdist = abs(read.reference_start-read.next_reference_start)
                if read.reference_id != read.next_reference_id: mdist=3e9

                if read.mapq >= min_mapq and mdist >= min_dist:
                    mchrom = bam.getrname(read.next_reference_id)
                    mstart = read.next_reference_start
                    mend   = mstart + len(read.seq)

                    mstr = '+'
                    if read.mate_is_reverse: mstr = '-'

                    if mchrom in forest:
                        for rec in forest[mchrom].find(mstart, mend):
                            coords.append(Coord(rchrom, rstart, rend, rstr, mchrom, mstart, mend, mstr, rec.value, os.path.basename(bam.filename)))
                            break

            if i % tick == 0:
                if read.is_unmapped:
                    logger.debug('parsed %d reads, last position unmapped' % i)
                else:
                    logger.debug('parsed %d reads, last position: %s:%d' % (i, bam.getrname(read.tid), read.pos))

    return coords


def subcluster_by_label(cluster):
    subclusters = dd(list)
    for c in cluster:
        subclusters[c.label.split('|')[3]].append(c)

    return subclusters.values()


def eval_strands(s):
    left = 0
    for i in range(1,len(s)):
        if s[i] != s[0]:
            left = i
            break

    right = 0
    for i in range(len(s)-1, 0, -1):
        if s[i] != s[-1]:
            right = i+1
            break

    if left == right: return left
    
    return 0


def filter_cluster(cluster):
    s1 = eval_strands([c.strand for c in cluster])
    s2 = eval_strands([c.mstrand for c in cluster])

    if s1 == s2 and s1 > 0: return False

    return True 


def infer_strand(cluster):
    c1 = [c.strand for c in cluster]
    c2 = [c.mstrand for c in cluster]

    if c1[0] == c2[0] and c1[-1] == c2[-1]: return '-'
    if c1[0] != c2[0] and c1[-1] != c2[-1]: return '+'

    return 'NA'


def output_cluster(cluster, forest, mapping, nonref, min_size=4):
    if len(cluster) >= min_size:
        cluster_chrom = cluster[0].chrom
        cluster_start = cluster[0].start
        if cluster_start < 0: cluster_start = 0

        cluster_end = cluster[-1].end

        bamlist = ','.join(list(set([c.bam for c in cluster])))

        if cluster_chrom not in forest or len(list(forest[cluster_chrom].find(cluster_start, cluster_end))) == 0:
            map_score = 0.0
            if mapping is not None:
                map_score = avgmap(mapping, cluster_chrom, cluster_start, cluster_end)

            if not filter_cluster(cluster) and (map_score >= 0.95 or mapping is None):
                nr = ['NA']

                if nonref is not None:
                    if cluster_chrom in nonref.contigs:
                        nr = ['|'.join(te.split()) for te in nonref.fetch(cluster_chrom, cluster_start, cluster_end)]
                    else:
                        nr = ['NA']

                    if not nr: nr = ['NA']

                nr = ','.join(nr)

                return InsCall(cluster, cluster_chrom, cluster_start, cluster_end, infer_strand(cluster), bamlist, map_score, nr)


def cluster(forest, coords, mapping, nonref, min_size=4, max_spacing=250):
    logger.debug('sorting coordinates')
    coords.sort()

    cluster = []
    insertion_list = []

    for c in coords:
        if len(cluster) == 0:
            cluster = [c]
        else:
            if c.chrom == cluster[-1].chrom and c.start - cluster[-1].end <= max_spacing:
                cluster.append(c)
            else:
                for cluster in subcluster_by_label(cluster):
                    i = output_cluster(cluster, forest, mapping, nonref, min_size=min_size)
                    if i is not None: insertion_list.append(i)
                cluster = [c]

    for cluster in subcluster_by_label(cluster):
        i = output_cluster(cluster, forest, mapping, nonref, min_size=min_size)
        if i is not None: insertion_list.append(i)

    return insertion_list


def run_chunk(args, chunk):
    ''' chunk is a list of (chrom, start, end) tuples '''

    bams = [pysam.AlignmentFile(bam, 'rb') for bam in args.bam.split(',')]

    mapping = None
    if args.mapping is not None:
        mapping = pysam.Tabixfile(args.mapping)

    nonref = None
    if args.nonref is not None:
        nonref = pysam.Tabixfile(args.nonref)

    logger.debug('building interval trees for %s' % args.bed)
    forest = interval_forest(args.bed)

    coords = []

    # operate over intervals in chunk
    for interval in chunk:
        chrom, start, end = interval

        logger.debug('%s:%d-%d: fetching coordinates from %s' % (chrom, start, end, args.bam))

        coords += get_coords(forest, bams, chrom=chrom, start=start, end=end)

        logger.debug('%s:%d-%d: found %d anchored reads' % (chrom, start, end, len(coords)))

    return cluster(forest, coords, mapping, nonref, min_size=int(args.minsize), max_spacing=int(args.maxspacing))


def resolve_dups(ins_list):
    ''' resolve cases where the same insertion has been called in multiple chunks '''
    ins_list.sort()
    new_list = []

    last = None

    for ins in ins_list:
        if last is None:
            last = ins

        elif last.overlaps(ins):
            if ins.length > last.length:
                last = ins

        else:
            new_list.append(last)
            last = ins

    if last is not None:
        new_list.append(last)

    return new_list


def main(args):
    logger.setLevel(logging.DEBUG)

    g = Genome(args.bam.split(',')[0])

    chunks = g.chunk(int(args.procs), pad=5000)

    pool = mp.Pool(processes=int(args.procs))

    reslist = []
    for chunk in chunks:
        res = res = pool.apply_async(run_chunk, [args, chunk])
        reslist.append(res)

    ins_list = []
    for res in reslist:
        ins_list += res.get()

    ins_list = resolve_dups(ins_list)
    for i in ins_list: print i.out()


if __name__ == '__main__':
    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(format=FORMAT)

    parser = argparse.ArgumentParser(description='identify clusters of discordant read ends where one end is in BED file')
    parser.add_argument('--bam', required=True, help='can be comma-delimited for multiple BAMs')
    parser.add_argument('--bed', required=True, help='locations of source locations (e.g. reference TEs) in genome')
    parser.add_argument('--mapping', default=None, help='mappability track tabix')
    parser.add_argument('--nonref', default=None, help='known nonreference element annotation')
    parser.add_argument('--minsize', default=4, help='minimum cluster size to output')
    parser.add_argument('--maxspacing', default=250, help='maximum spacing between support reads (default=250)')
    parser.add_argument('-p', '--procs', default=1, help='split work over multiple processes')
  
    args = parser.parse_args()
    main(args)
