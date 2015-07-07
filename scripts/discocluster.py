#!/usr/bin/env python

import os
import pysam
import argparse
import logging
logger = logging.getLogger(__name__)

from collections import defaultdict as dd
from bx.intervals.intersection import Intersecter, Interval




''' identify clusters of discordant read ends where one end is in BED file '''


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


def get_coords(forest, bams, min_mapq=1, min_dist=10000):
    coords = []

    for bam in bams:
        tick = 10000000
        try:
            tick = int((bam.mapped + bam.unmapped) * 0.01)
            if tick == 0: tick = 1
            logger.debug('outputting status every %d reads (1 pct)' % tick)

        except ValueError as e:
            logger.debug('no index found, outputting status every %d reads' % tick)

        for i, read in enumerate(bam.fetch()):
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


def eval_cluster(cluster):
    ''' check for strand switch, strand consistency '''
    pass


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

                print '#BEGIN'
                print '%s\t%d\t%d\t%s\t%s\t%d\t%0.3f\t%s' % (cluster_chrom, cluster_start, cluster_end, infer_strand(cluster), bamlist, len(cluster), map_score, nr)
                for c in cluster: print c
                print '#END\n'


def cluster(forest, coords, mapping, nonref, min_size=4, max_spacing=250):
    logger.debug('sorting coordinates')
    coords.sort()

    cluster = []

    for c in coords:
        #print c
        if len(cluster) == 0:
            cluster = [c]
        else:
            if c.chrom == cluster[-1].chrom and c.start - cluster[-1].end <= max_spacing:
                cluster.append(c)
            else:
                for cluster in subcluster_by_label(cluster):
                    output_cluster(cluster, forest, mapping, nonref, min_size=min_size)
                cluster = [c]

    for cluster in subcluster_by_label(cluster):
        output_cluster(cluster, forest, mapping, nonref, min_size=min_size)


def main(args):
    logger.setLevel(logging.DEBUG)
    assert args.bam.endswith('.bam'), "not a BAM file: %s" % args.bam
    bams = [pysam.AlignmentFile(bam, 'rb') for bam in args.bam.split(',')]

    mapping = None
    if args.mapping is not None:
        mapping = pysam.Tabixfile(args.mapping)

    nonref = None
    if args.nonref is not None:
        nonref = pysam.Tabixfile(args.nonref)

    logger.debug('building interval trees for %s' % args.bed)
    forest = interval_forest(args.bed)

    logger.debug('fetching coordinates from %s' % args.bam)
    coords = get_coords(forest, bams)
    logger.debug('found %d anchored reads' % len(coords))

    cluster(forest, coords, mapping, nonref, min_size=int(args.minsize), max_spacing=int(args.maxspacing))


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
  
    args = parser.parse_args()
    main(args)
