#!/usr/bin/env python

import os
import gzip
import pysam
import subprocess
import argparse
import logging

from string import maketrans

FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def rc(dna):
    ''' reverse complement '''
    complements = maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
    return dna.translate(complements)[::-1]


def bwa_align(fq, ref, outbase, threads=1, mem='1G'):
    sam_cmd  = ['bwa', 'mem', '-t', str(threads), ref, fq]
    view_cmd = ['samtools', 'view', '-Su', '-']
    sort_cmd = ['samtools', 'sort', '-@', str(threads), '-m', mem, '-', '-T', outbase, '-o', outbase+'.bam']

    aln  = subprocess.Popen(sam_cmd, stdout=subprocess.PIPE)
    view = subprocess.Popen(view_cmd, stdin=aln.stdout, stdout=subprocess.PIPE)
    sort = subprocess.Popen(sort_cmd, stdin=view.stdout, stdout=subprocess.PIPE)

    for line in sort.stdout:
        pass

    subprocess.call(['samtools', 'index', outbase+'.bam'])

    return outbase+'.bam'


def genedict(fn):
    infile = None
    if fn.endswith('.gz'):
        infile = gzip.open(fn, 'rb')
    else:
        infile = open(fn, 'r')

    gdict = {}
    for line in infile:
        c = line.strip().split()
        name   = c[12]
        chrom  = c[2]
        strand = c[3]
        start  = int(c[4]) - 1000
        end    = int(c[5]) + 1000

        gdict[name] = ['%s:%d-%d' % (chrom, start, end), strand]

    infile.close()

    return gdict


def callmuts(bam, ref, outbase, region, gene):
    samtools_cmd = ['samtools', 'mpileup', '-r', region[0], '-ugf', ref, bam]
    bcftools_cmd = ['bcftools', 'call', '-vm']

    st = subprocess.Popen(samtools_cmd, stdout=subprocess.PIPE)
    bt = subprocess.Popen(bcftools_cmd, stdin=st.stdout, stdout=subprocess.PIPE)

    with open(outbase+'.vcf', 'w') as vcf:
        for line in bt.stdout:
            if line.startswith('#CHROM'):
                vcf.write('##GENE=%s\n' % gene)
                vcf.write('##STRAND=%s\n' % region[1])
            if not line.startswith('#'):
                chrom, pos, id, ref, alt, qual, filter, info, format, sample = line.strip().split('\t')

                #if region[1] == '-':
                #    ref = rc(ref)
                #    alt = rc(alt)

                if len(ref) == len(alt):
                    print '\t'.join((chrom, pos, ref, alt, qual, filter, info, gene, region[1]))

            vcf.write(line)

    return outbase+'.vcf'


def main(args):
    assert os.path.exists(args.path), 'path not found: %s' % args.path

    genes = genedict(args.genes)

    with open(args.tabfile, 'r') as tab:
        for i, line in enumerate(tab):
            if i == 0: # header
                header = line.strip().split('\t')
                continue

            rec = {}
            for n, field in enumerate(line.strip().split('\t')):
                rec[header[n]] = field

            bamfn = '%s/tebreak.%s.discoremap.bam' % (args.path, rec['UUID'])
            fastq = '%s/tebreak.%s.discoremap.fastq.gz' % (args.path, rec['UUID'])

            if os.path.exists(bamfn):
                bam = pysam.AlignmentFile(bamfn, 'rb')

                if bam.mapped > int(args.minmap):
                    
                    with gzip.open(fastq, 'wb') as out:
                        for i, read in enumerate(bam.fetch()):
                            name = rec['UUID'] + '.' + read.qname

                            seq = read.seq
                            if read.is_reverse:
                                seq = rc(seq)

                            out.write('@%s\n%s\n+\n%s\n' % (name, seq, read.qual))

                    logger.info('wrote %d reads to %s' % (i, fastq))

                    outbam = bwa_align(fastq, args.ref, '%s/tebreak.%s.genomic' % (args.path, rec['UUID']))

                    gene_name = '.'.join(rec['Superfamily'].split('.')[:-1])
                    if gene_name in genes:
                        outvcf = callmuts(outbam, args.ref, '%s/tebreak.%s.genomic' % (args.path, rec['UUID']), genes[gene_name], gene_name)
                        logger.info('generated callset %s on region %s' % (outvcf, genes[gene_name][0]))

                    else:
                        logger.info('no entry in %s for %s' % (args.genes, gene_name))

                else:
                    logger.info('less than %d reads for UUID: %s' % (int(args.minmap), rec['UUID']))

            else:
                logger.info('BAM not found: %s' % bamfn)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="remap intra-insertion reads to reference")
    parser.add_argument('-t', '--tabfile', required=True)
    parser.add_argument('-p', '--path', required=True, help='path to TEBreak/resolve.py sub-BAMs')
    parser.add_argument('-r', '--ref', required=True, help='reference genome')
    parser.add_argument('-g', '--genes', required=True, help='refGene.txt.gz')
    parser.add_argument('--minmap', default=10)
    args = parser.parse_args()
    main(args)
