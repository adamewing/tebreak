#!/usr/bin/env python

from __future__ import print_function

import os
import pysam
import random
import logging
import argparse
import itertools
import subprocess

logger = logging.getLogger(__name__)

from uuid import uuid4

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

        if start < 0:
            start = 0

        if end > self.chrlen[chrom]:
            end = self.chrlen[chrom]

        return (chrom, start, end)


    def chunk(self, n, seed=None, sorted=True, pad=0):
        ''' break genome into n evenly-sized chunks, return n lists of (chrom, start, end) '''
        chunklen = int(self.bp/n)
        
        chunks = []
        intervals = []
 
        chunkleft = chunklen # track how much genome needs to go into each chunk
 
        chromlist = list(self.chrlen.keys())
 
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
                intervals.append(self.addpad((chrom, 0, length), pad))
                assert lenleft == 0
 
                if chunkleft == 0:
                    chunkleft = chunklen
                    chunks.append(intervals)
                    intervals = []
            else:
                while lenleft > 0:
                    if lenleft >= chunkleft:
                        intervals.append(self.addpad((chrom, length-lenleft, length-lenleft+chunkleft), pad))
                        lenleft -= chunkleft
 
                        chunkleft = chunklen
                        chunks.append(intervals)
                        intervals = []
 
                    else: # lenleft < chunkleft
                        intervals.append(self.addpad((chrom, length-lenleft, length), pad))
                        chunkleft -= lenleft
                        lenleft -= lenleft
 
        return chunks


class DiscoCoord:
    def __init__(self, chrom, start, end, strand, mchrom, mstart, mend, mstrand, label, rname):
        self.chrom   = chrom
        self.start   = int(start)
        self.end     = int(end)
        self.strand  = strand
        self.mchrom  = mchrom
        self.mstart  = int(mstart)
        self.mend    = int(mend)
        self.mstrand = mstrand
        self.label   = label
        self.rname   = rname

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
        return '\t'.join(map(str, (self.rname, self.label, self.chrom, self.start, self.end, self.strand, self.mchrom, self.mstart, self.mend, self.mstrand)))


class DiscoInsCall:
    def __init__(self, coord_list, chrom, start, end, strand, mapscore, nonref):
        self.uuid = str(uuid4())
        self.coord_list = coord_list
        self.chrom      = chrom
        self.start      = int(start)
        self.end        = int(end)
        self.strand     = strand
        self.length     = len(coord_list)
        self.mapscore   = mapscore
        self.nonref     = nonref
        self.clip_reads = []
        self.align_count = 0


    def getlabel(self):
        return ','.join(list(set([coord.label.split('|')[4] for coord in self.coord_list])))


    def out(self, verbose=False):
        output = []

        if verbose:
            output = ['#BEGIN']

        output.append('%s\t%s\t%d\t%d\t%s\t%d\t%0.3f\t%d\t%s' % (self.getlabel(), self.chrom, self.start, self.end, self.strand, self.length, self.mapscore, self.align_count, self.nonref))

        if verbose:
            for c in self.coord_list:
                output.append(str(c))

            output.append('#END')

        return '\n'.join(output)


    def overlaps(self, other):
        ''' return true if overlap > 0 '''
        return min(self.end, other.end) - max(self.start, other.start) > 0


    def align_split(self, bam, teindex, minclip=15, max_altclip=2):
        sr = []
        bp = []

        assert os.path.exists(teindex + '.sa'), 'not indexed: %' % teindex

        with open(self.uuid + '.fq', 'w') as out:
            for read in bam.fetch(self.chrom, self.start, self.end):

                if None in (read.rlen, read.alen):
                    continue

                if read.rlen - read.alen >= minclip: # soft-clipped
                    altclip = min(read.qstart, read.rlen-read.qend)

                    if altclip > max_altclip:
                        continue

                    sr.append(SplitRead(read))
                    bp.append(sr[-1].genome_breakpos)

                    out.write('@%s\n%s\n+\n%s\n' % (sr[-1].uuid,sr[-1].unmapped_seq(),sr[-1].unmapped_qual()))

        bwa_cmd = ['bwa', 'mem', '-k', str(minclip), teindex, self.uuid+'.fq']

        FNULL = open(os.devnull, 'w')
        aln = subprocess.Popen(bwa_cmd, stdout=subprocess.PIPE, stderr=FNULL)

        self.align_count = 0

        for line in aln.stdout:
            line=line.decode()
            if line.startswith('@'):
                continue
            c = line.strip().split('\t')
            if c[2] != '*':
                self.align_count += 1

            elif c[9].startswith('A'*10): # poly-A might not align well but don't want to penalise
                self.align_count += 1

        os.remove(self.uuid+'.fq')

        return self.align_count


    def __gt__(self, other):
        if self.chrom == other.chrom:
            return self.start > other.start
        else:
            return self.chrom > other.chrom


    def __str__(self):
        return self.out(verbose=False)



class SplitRead:
    def __init__(self, read):
        self.uuid  = str(uuid4())
        self.read  = read

        self.cliplen = len(read.seq) - len(read.query_alignment_sequence)

        self.breakleft  = False
        self.breakright = False

        self.genome_breakpos = None
        self.query_breakpos  = None
 
        if read.qstart < read.rlen - read.qend:
            self.genome_breakpos = read.get_reference_positions()[-1] # breakpoint on right
            self.query_breakpos = read.query_alignment_end
            self.breakright = True
 
        else:
            self.genome_breakpos = read.get_reference_positions()[0] # breakpoint on left
            self.query_breakpos = read.query_alignment_start
            self.breakleft = True
 
        assert None not in (self.genome_breakpos, self.query_breakpos)
        assert self.breakleft != self.breakright
 

    def __gt__(self, other):
        return self.breakpos > other.breakpos


    def unmapped_seq(self):
        if self.breakright:
            return self.read.seq[self.query_breakpos:]
        else:
            return self.read.seq[:self.query_breakpos]

    def unmapped_qual(self):
        if self.breakright:
            return self.read.qual[self.query_breakpos:]
        else:
            return self.read.qual[:self.query_breakpos]


def avgmap(maptabix, chrom, start, end):
    ''' return average mappability across chrom:start-end region; maptabix = pysam.Tabixfile '''
    scores = []

    if None in (start, end):
        return None

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


def flip(strand):
    if strand == '+':
        return '-'

    if strand == '-':
        return '+'


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


def disco_get_coords(forest, bam, chrom=None, start=None, end=None, min_mapq=1, min_dist=10000):
    coords = []

    logger.info('parsing %s:%d-%d' % (chrom, start, end))

    for i, read in enumerate(read_gen(bam, chrom=chrom, start=start, end=end)):
        if not read.is_unmapped and not read.mate_is_unmapped and not read.is_duplicate:

            rchrom = bam.getrname(read.reference_id)
            rstart = read.reference_start
            rend   = read.reference_end

            rstr = '+'
            if read.is_reverse:
                rstr = '-'

            mdist = abs(read.reference_start-read.next_reference_start)

            if read.reference_id != read.next_reference_id:
                mdist=3e9

            if read.mapq >= min_mapq and mdist >= min_dist:
                mchrom = bam.getrname(read.next_reference_id)
                mstart = read.next_reference_start
                mend   = mstart + len(read.seq)

                mstr = '+'
                if read.mate_is_reverse:
                    mstr = '-'

                if mchrom in forest:
                    for rec in forest[mchrom].find(mstart, mend):
                        coords.append(DiscoCoord(rchrom, rstart, rend, rstr, mchrom, mstart, mend, mstr, rec.value, read.qname))
                        break

    return coords


def disco_subcluster_by_label(cluster):
    subclusters = dd(list)
    for c in cluster:
        subclusters[c.label.split('|')[3]].append(c)

    return subclusters.values()


def disco_eval_strands(s):
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

    if left == right:
        return left
    
    return 0


def disco_infer_strand(cluster):
    c1 = [c.strand for c in cluster]
    c2 = [c.mstrand for c in cluster]

    if c1[0] == c2[0] and c1[-1] == c2[-1]: return '-'
    if c1[0] != c2[0] and c1[-1] != c2[-1]: return '+'

    return 'NA'


def disco_output_cluster(cluster, forest, maptrack, nonref, min_size=4, min_map=0.5):
    if len(cluster) >= min_size:
        cluster_chrom = cluster[0].chrom
        cluster_start = cluster[0].start


        if cluster_start < 0:
            cluster_start = 0

        cluster_end = cluster[-1].end
        
        map_score = 0.0

        if maptrack is not None:
            map_score = avgmap(maptrack, cluster_chrom, cluster_start, cluster_end)

        if map_score >= float(min_map) or maptrack is None:
            nr = ['NA']

            if nonref is not None:
                if cluster_chrom in nonref.contigs:
                    nr = ['|'.join(te.split()) for te in nonref.fetch(cluster_chrom, cluster_start, cluster_end)]
                else:
                    nr = ['NA']

                if not nr:
                    nr = ['NA']

            nr = ','.join(nr)

            return DiscoInsCall(cluster, cluster_chrom, cluster_start, cluster_end, disco_infer_strand(cluster), map_score, nr)


def disco_cluster(forest, bam, teindex, coords, maptrack, nonref, min_size=4, min_map=0.5, max_spacing=250):
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
                for cluster in disco_subcluster_by_label(cluster):
                    i = disco_output_cluster(cluster, forest, maptrack, nonref, min_size=min_size)

                    if not i:
                        continue

                    i.align_split(bam, teindex)
                    insertion_list.append(i)

                cluster = [c]


    for cluster in disco_subcluster_by_label(cluster):
        i = disco_output_cluster(cluster, forest, maptrack, nonref, min_size=min_size, min_map=min_map)

        # find split ends and align
        if not i:
            continue

        i.align_split(bam, teindex)

        insertion_list.append(i)

    return insertion_list


def disco_run_chunk(args, chunk):
    ''' chunk is a list of (chrom, start, end) tuples '''

    bam = pysam.AlignmentFile(args.bam, 'rb')

    maptrack = None
    if args.maptrack is not None:
        maptrack = pysam.Tabixfile(args.maptrack)

    nonref = None
    if args.nonref is not None:
        nonref = pysam.Tabixfile(args.nonref)

    logger.debug('building interval trees for %s' % args.rmsk)
    forest = interval_forest(args.rmsk)

    coords = []

    # operate over intervals in chunk
    for interval in chunk:
        chrom, start, end = interval

        logger.debug('%s:%d-%d: fetching coordinates from %s' % (chrom, start, end, args.bam))

        coords += disco_get_coords(forest, bam, chrom=chrom, start=start, end=end)

        logger.debug('%s:%d-%d: found %d anchored reads' % (chrom, start, end, len(coords)))

    return disco_cluster(forest, bam, args.teindex, coords, maptrack, nonref, min_size=int(args.minsize), min_map=float(args.minmap), max_spacing=int(args.maxspacing))


def disco_resolve_dups(ins_list):
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

    chunks = g.chunk(int(args.procs), pad=2500)

    pool = mp.Pool(processes=int(args.procs))

    reslist = []
    for chunk in chunks:
        res = res = pool.apply_async(disco_run_chunk, [args, chunk])
        reslist.append(res)

    ins_list = []
    for res in reslist:
        ins_list += res.get()

    if args.deduplicate:
        ins_list = disco_resolve_dups(ins_list)
    
    output_reads = []

    for i in ins_list:
        if i.align_count >= int(args.min_align_count):
            print(i.out())

            output_reads += [coord.rname for coord in i.coord_list]
            output_reads += i.clip_reads

    output_reads = set(output_reads)

    # bam output

    inbam  = pysam.AlignmentFile(args.bam, 'rb')
    outbam = pysam.AlignmentFile(args.outbam, 'wb', template=inbam)

    logger.output('writing to %s...' % args.outbam)

    for read in inbam.fetch():
        if read.qname in output_reads:
            outbam.write(read)

    inbam.close()
    outbam.close()


if __name__ == '__main__':
    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(format=FORMAT)

    parser = argparse.ArgumentParser(description='identify clusters of discordant read ends where one end is in BED file')
    parser.add_argument('-b', '--bam', required=True, help='Input BAM')
    parser.add_argument('-o', '--outbam', required=True, help='Output (reduced) BAM')
    parser.add_argument('-r', '--rmsk', required=True, help='locations of source locations (e.g. reference TEs) in genome: chrom, start, end, class, family, orientation')
    parser.add_argument('--teindex', required=True, help='bwa indexed te reference')
    parser.add_argument('--maptrack', default=None, help='mappability track tabix')
    parser.add_argument('--minmap', default=0.5, help='minimum region mappability score (default = 0.5)')
    parser.add_argument('--nonref', default=None, help='known nonreference element annotation')
    parser.add_argument('--minsize', default=4, help='minimum cluster size to output (default = 4)')
    parser.add_argument('--maxspacing', default=250, help='maximum spacing between support reads (default=250)')
    parser.add_argument('--min_align_count', default=2, help='minimum number of split reads mapping to te reference (default = 2)')
    parser.add_argument('-p', '--procs', default=1, help='split work over multiple processes')
    parser.add_argument('--deduplicate', action='store_true', default=False)
  
    args = parser.parse_args()
    main(args)
