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

    return outdir + '/' + 'matched.random.txt.gz', n


def overlap(bamfn, bedfn, outdir, tabix=False):
    bam = None

    if tabix:
        bam = pysam.Tabixfile(bamfn)
    else:
        bam = pysam.Samfile(bamfn, 'rb')

    overlaps = {}

    with open(bedfn, 'r') as bed:
        for line in bed:
            chrom, start, end, strand, name, eltstart, eltend = line.strip().split()
            elt5p, elt3p, eltstart, eltend = map(int, (start, end, eltstart, eltend))

            if strand == '-':
                elt5p, elt3p = elt3p, elt5p

            # initialize counters
            if tabix:
                libs = ['sampled']
            else:
                assert 'RG' in bam.header
                libs = [rg['ID'] for rg in bam.header['RG']]

            for lib in libs:
                if lib not in overlaps:
                    overlaps[lib] = {}

                overlaps[lib][line.strip()] = {}
                overlaps[lib][line.strip()]['5pNonDup'] = 0
                overlaps[lib][line.strip()]['3pNonDup'] = 0
                overlaps[lib][line.strip()]['5pAll'] = 0
                overlaps[lib][line.strip()]['3pAll'] = 0


            # count 5p junction overlaps
            for read in bam.fetch(chrom, elt5p, elt5p+1):
                lib = 'sampled'

                if not tabix:
                    lib = getRG(read.tags)

                if tabix or not read.is_duplicate:
                    overlaps[lib][line.strip()]['5pNonDup'] += 1

                overlaps[lib][line.strip()]['5pAll'] += 1

            # count 3p junction overlaps
            for read in bam.fetch(chrom, elt3p, elt3p+1):
                lib='sampled'

                if not tabix:
                    lib = getRG(read.tags)

                if tabix or not read.is_duplicate:
                    overlaps[lib][line.strip()]['3pNonDup'] += 1

                overlaps[lib][line.strip()]['3pAll'] += 1

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

    rndfn, rndcount = makerandom(args.bam, args.ref, args.outdir, samples=int(args.samples))
    real_overlaps = overlap(args.bam, args.elts, args.outdir)
    sim_overlaps  = overlap(rndfn, args.elts, args.outdir, tabix=True)

    normalised = od()
    enrichment = od()

    counts = readcount(args.bam)
    counts['sampled'] = rndcount

    for lib in real_overlaps:
        normalised[lib] = {}
        for site in real_overlaps[lib]:
            normalised[lib][site] = {}
            for category in real_overlaps[lib][site]:
                normalised[lib][site][category] = float(real_overlaps[lib][site][category])/float(counts[lib])

    for lib in sim_overlaps:
        normalised[lib] = {}
        for site in sim_overlaps[lib]:
            normalised[lib][site] = {}
            for category in sim_overlaps[lib][site]:
                normalised[lib][site][category] = float(sim_overlaps[lib][site][category])/float(counts[lib])

    for lib in normalised:
        enrichment[lib] = {}
        for site in normalised[lib]:
            enrichment[lib][site] = {}
            for category in normalised[lib][site]:
                if normalised['sampled'][site][category] == 0:
                    enrichment[lib][site][category] = 'inf'
                else:
                    enrichment[lib][site][category] = normalised[lib][site][category] / normalised['sampled'][site][category]

    bam = pysam.Samfile(args.bam, 'rb')
    libs = normalised.keys()

    for lib in real_overlaps.keys():
        with open(args.outdir + '/' + lib + '.enrichment.txt', 'w') as libout:
            libout.write("Total reads in library:\t" + str(counts[lib]) + "\n")
            libout.write("Sampled reads in matched random:\t" + str(counts['sampled']) + "\n")
            libout.write("Chrom\tStart\tEnd\tStrand\tID\tTE_Start\tTE_End\t5p Overlaps\t5p Matched Random Overlaps\t5p Overlap Enrichment\t5p Overlaps (without PCR dups)\t5p NonDup Matched Random Overlaps\t5p NonDup Overlap Enrichment\t3p Overlaps\t3p Matched Random Overlaps\t3p Overlap Enrichment\t3p Overlaps (without PCR dups)\t3p NonDup Matched Random Overlaps\t3p NonDup Overlap Enrichment\n")
            with open(args.elts, 'r') as elts:
                for line in elts:
                    site = line.strip()
                    libout.write('\t'.join((site, 
                                     str(real_overlaps[lib][site]['5pAll']),
                                     str(sim_overlaps['sampled'][site]['5pAll']),
                                     str(enrichment[lib][site]['5pAll']),
                                     str(real_overlaps[lib][site]['5pNonDup']),
                                     str(sim_overlaps['sampled'][site]['5pNonDup']),
                                     str(enrichment[lib][site]['5pNonDup']),
                                     str(real_overlaps[lib][site]['3pAll']),
                                     str(sim_overlaps['sampled'][site]['3pAll']),
                                     str(enrichment[lib][site]['3pAll']),
                                     str(real_overlaps[lib][site]['3pNonDup']),
                                     str(sim_overlaps['sampled'][site]['3pNonDup']),
                                     str(enrichment[lib][site]['3pNonDup'])
                                     )) + '\n')


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
