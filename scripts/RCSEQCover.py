#!/usr/bin/env python

import pysam
import sys
import os
import argparse
import subprocess

from random import randint, seed
from collections import defaultdict as dd
from collections import OrderedDict as od
from shutil import move
from uuid import uuid4


def getRG(tags):
    ''' fetch RG, tags is from a pysam.AlignedRead.tags, returns RG name '''
    for tag, val in tags:
        if tag == 'RG':
            return val
    return None


def makerandom(bamfn, reffn, outdir, samples=1000000):
    assert os.path.exists(reffn + '.fai'), "reference not indexed: " + reffn
    ref = pysam.Fastafile(reffn)
    bam = pysam.Samfile(bamfn, 'rb')

    # connect chromosomes end-to-end and track offsets, otherwise sampling will be biased
    chrlen = dd(int)
    for line in open(reffn + '.fai', 'r'):
        chrom, size = line.strip().split()[:2]
        size = int(size)
        chrlen[chrom] = size

    # calculate offsets
    offset = 0
    chroffset = dd(int)
    for chrom in sorted(chrlen.iterkeys()):
        offset += int(chrlen[chrom])
        chroffset[chrom] = offset

    genomelen = offset

    with open(outdir + '/' + 'matched.random.txt', 'w') as mrout:
        n = 1
        for read in bam.fetch(until_eof=True):
            rndoffset = randint(0,genomelen)
            rndchrom = None
            lastoffset = 0
            for chrom in sorted(chrlen.iterkeys()):
                offset = int(chroffset[chrom])
                assert lastoffset < offset
                if rndoffset >= lastoffset and rndoffset < offset:
                    rndchrom = chrom
                    rndloc = rndoffset - lastoffset
                lastoffset = offset

            rndstart, rndend = rndloc, rndloc+read.rlen

            # don't extend random intervals past the end of the chromosome
            if rndloc+read.rlen > chrlen[rndchrom]:
                rndstart, rndend = rndloc-read.rlen, rndloc

            assert rndstart < rndend
            mrout.write('\t'.join((rndchrom, str(rndstart), str(rndend))) + '\n')

            if n == samples:
                break

            if n % 1000000 == 0:
                print n,"segments sampled..."

            n += 1


    tmpfn = outdir + '/' + str(uuid4()) + '.tmp'

    with open(tmpfn, 'w') as tmp:
        p = subprocess.Popen(['sort', '-k1,1', '-k2,2n', outdir + '/' + 'matched.random.txt'], stdout=subprocess.PIPE)
        for line in p.stdout:
            tmp.write(line)

    move(tmpfn, outdir + '/' + 'matched.random.txt')

    subprocess.call(['bgzip', '-f', outdir + '/' + 'matched.random.txt'])
    assert os.path.exists(outdir + '/' + 'matched.random.txt.gz'), "bgzip failed: " + outdir + '/' + 'matched.random.txt'

    subprocess.call(['tabix', '-s', '1', '-b', '2', '-e', '3', '-f', outdir + '/' + 'matched.random.txt.gz'])
    assert os.path.exists(outdir + '/' + 'matched.random.txt.gz.tbi'), "tabix failed: " + outdir + '/' + 'matched.random.txt'

    return outdir + '/' + 'matched.random.txt.gz'


def overlap(bamfn, bedfn, outdir, tabix=False):
    bam = None

    if tabix:
        bam = pysam.Tabixfile(bamfn)
    else:
        bam = pysam.Samfile(bamfn, 'rb')

    overlaps = {}

    with open(bedfn, 'r') as bed:
        for line in bed:
            lib = 'sampled'

            chrom, start, end, strand, name, eltstart, eltend = line.strip().split()
            elt5p, elt3p, eltstart, eltend = map(int, (start, end, eltstart, eltend))

            if strand == '-':
                elt5p, elt3p = elt3p, elt5p

            # count 5p junction overlaps
            for read in bam.fetch(chrom, elt5p, elt5p+1):
                info = line.strip()
                if not tabix:
                    info = '\t'.join((line.strip(), str(read.positions[0]), str(read.positions[-1])))

                if not tabix:
                    lib = getRG(read.tags)
                if lib not in overlaps:
                    overlaps[lib] = dd(list)

                if tabix or not read.is_duplicate:
                    overlaps[lib]['5pNonDup'].append(info)

                overlaps[lib]['5pAll'].append(info)

            # count 3p junction overlaps
            for read in bam.fetch(chrom, elt3p, elt3p+1):
                info = line.strip()
                if not tabix:
                    info = '\t'.join((line.strip(), str(read.positions[0]), str(read.positions[-1])))

                if not tabix:
                    lib = getRG(read.tags)
                if lib not in overlaps:
                    overlaps[lib] = dd(list)

                if tabix or not read.is_duplicate:
                    overlaps[lib]['3pNonDup'].append(info)

                overlaps[lib]['3pAll'].append(info)


    #for lib in overlaps.keys():
        #print lib
        #print "*"*30 + " debug " + "*"*30
        #print '\t'.join(overlaps[lib].keys())
        #print '\t'.join(map(str, [len(v) for v in overlaps[lib].values()]))
        #print "*"*67

    return overlaps


def readcount(bamfn):
    bam = pysam.Samfile(bamfn, 'rb')
    counts = dd(int)
    for read in bam.fetch(until_eof=True):
        rg = getRG(read.tags)
        counts[rg] += 1

    return counts


def main(args):
    assert args.bam.endswith('.bam'), "not a BAM file: " + args.bam
    seed(int(args.seed))

    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)
    assert os.path.exists(args.outdir), 'could not create output directory: ' + args.outdir

    rndfn = makerandom(args.bam, args.ref, args.outdir, samples=int(args.samples))
    real_overlaps = overlap(args.bam, args.elts, args.outdir)
    sim_overlaps  = overlap(rndfn, args.elts, args.outdir, tabix=True)

    normalised = od()
    enrichment = od()

    counts = readcount(args.bam)
    counts['sampled'] = int(args.samples)

    for lib in real_overlaps:
        normalised[lib] = {}
        for category in real_overlaps[lib]:
            normalised[lib][category] = float(len(real_overlaps[lib][category]))/float(counts[lib])

    for lib in sim_overlaps:
        normalised[lib] = {}
        for category in sim_overlaps[lib]:
            normalised[lib][category] = float(len(sim_overlaps[lib][category]))/float(counts[lib])

    for lib in normalised:
        enrichment[lib] = {}
        for category in normalised[lib]:
            enrichment[lib][category]  = normalised[lib][category] / normalised['sampled'][category]

    with open(args.outdir + '/' + 'enrichment.txt', 'w') as enrich_out:
        enrich_out.write("Library\t5p Overlaps\t5p Overlaps (without PCR dups)\t3p Overlaps\t3p Overlaps (without PCR dups)\n")
        for lib in enrichment:
            enrich_out.write('\t'.join((lib, 
                             str(enrichment[lib]['5pAll']),
                             str(enrichment[lib]['5pNonDup']),
                             str(enrichment[lib]['3pAll']),
                             str(enrichment[lib]['3pNonDup'])
                             )) + '\n')

    for lib in real_overlaps:
        for category in real_overlaps[lib]:
            with open(args.outdir + '/' + category + '.' + lib + '.txt', 'w') as libout:
                for rec in real_overlaps[lib][category]:
                    libout.write(rec + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Examine RC-Seq enrichment vs random sampling')
    parser.add_argument('-r', '--ref', required=True, help='Reference Genome (samtools indexed)')
    parser.add_argument('-b', '--bam', required=True, help='BAM file')
    parser.add_argument('-e', '--elts', required=True, help='Element coordinates (chrom, start, end, strand, name, eltstart, eltend')
    parser.add_argument('-o', '--outdir', required=True, help='output directory')
    parser.add_argument('-s', '--samples', default=10000000, help='maximum number of sites to sample for random distribution')
    parser.add_argument('--seed', default=123, help='seed for random number generator (default=123)')
    args = parser.parse_args()
    main(args)
