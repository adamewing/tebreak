#!/usr/bin/env python

import sys, os
import argparse
import subprocess
import pysam
import datetime
import gzip

from re import sub
from uuid import uuid4


def now():
    return str(datetime.datetime.now())


def bamtofastq(bam, outdir, samtofastq, threads=1):
    assert os.path.exists(samtofastq)
    assert bam.endswith('.bam')

    outfq = sub('bam$', 'fastq', outdir + '/' + os.path.basename(bam))

    cmd = ['java', '-XX:ParallelGCThreads=' + str(threads), '-Xmx4g', '-jar', samtofastq, 
           'INPUT=' + bam, 'INTERLEAVE=true', 'VALIDATION_STRINGENCY=SILENT', 'FASTQ=' + outfq]
    sys.stderr.write("INFO\t" + now() + "\tconverting BAM " + bam + " to FASTQ\n")
    subprocess.call(cmd)

    assert os.path.exists(outfq) # conversion failed

    return outfq


def cramtofastq(cram, outdir, cramjar, ref, threads=1):
    assert os.path.exists(cramjar)
    assert cram.endswith('.cram')

    outfq = sub('cram$', 'fastq', outdir + '/' + os.path.basename(cram))

    cmd = ['java', '-XX:ParallelGCThreads=' + str(threads), '-Xmx4g', '-jar', cramjar, 
           'fastq', '-I', cram, '-R', ref, '--enumerate', '--skip-md5-check']
    sys.stderr.write("INFO\t" + now() + "\tconverting CRAM " + cram + " to FASTQ\n")

    with open(outfq, 'w') as fq:
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        for line in p.stdout:
            fq.write(line)

    assert os.path.exists(outfq) # conversion failed

    return outfq


def bwamem(fq, ref, outdir, threads=1, width=150, uid=None):
    ''' FIXME: add parameters to commandline '''
    assert os.path.exists(ref + '.fai')
    
    fqroot = sub('fastq$', '', fq)
    if uid is not None:
        fqroot = uid

    fqroot   = outdir + '/' + os.path.basename(fqroot)

    sam_out   = fqroot + 'capture.sam'
    bam_out   = fqroot + 'capture.bam'
    sort_out  = fqroot + 'capture.sort'

    sam_cmd  = ['bwa', 'mem', '-t', str(threads), '-M', ref, fq]
    bam_cmd  = ['samtools', 'view', '-bt', ref + '.fai', '-o', bam_out, sam_out]

    sys.stderr.write("INFO\t" + now() + "\trunning bwa-mem: " + ' '.join(sam_cmd) + "\n")

    with open(sam_out, 'w') as sam:
        p = subprocess.Popen(sam_cmd, stdout=subprocess.PIPE)
        for line in p.stdout:
            sam.write(line)

    sys.stderr.write("INFO\t" + now() + "\twriting " + sam_out + " to BAM...\n")
    subprocess.call(bam_cmd)

    os.remove(sam_out)
    return bam_out


def getmatch(read):
    ''' return number of mismatches / aligned length of TE sub-read '''
    nm = [value for (tag, value) in read.tags if tag == 'NM'][0]
    return 1.0 - (float(nm)/float(read.alen))


def parsebam(inbam, fqout, minmatch=0.0):
    # pass 1, identify hardclipped reads
    hc = {}
    bam = pysam.Samfile(inbam, 'rb')
    for read in bam.fetch(until_eof=True):
        if read.is_secondary:
            hc[read.qname] = True

    sys.stderr.write("INFO\t" + now() + "\tIdentified " + str(len(hc)) + " hardclipped reads\n")

    # pass 2, output reads
    found_hc = 0
    outcount = 0
    bam.close()

    sys.stderr.write("INFO\t" + now() + "\toutputting to " + fqout + "\n")
    with gzip.open(fqout, 'wb') as fq:
        bam = pysam.Samfile(inbam, 'rb')
        for read in bam.fetch(until_eof=True):
            if read.qname in hc and not read.is_secondary:
                fq.write(fqread(read))
                found_hc += 1
                outcount += 1
            elif getmatch(read) >= float(minmatch) and not read.is_secondary:
                fq.write(fqread(read))
                outcount += 1

    sys.stderr.write("INFO\t" + now() + "\toutput " + str(outcount) + " reads, " + str(found_hc) + " of them rescued from hardclips\n")


def fqread(read):
    salt = str(uuid4()).split('-')[0]
    return '\n'.join(('@' + read.qname + ':' + salt, read.seq, '+', read.qual)) + '\n'


def main(args):
    sys.stderr.write("INFO\t" + now() + "\tstarting " + sys.argv[0] + " called with args:\n" + ' '.join(sys.argv) + "\n")

    if not os.path.exists(args.outdir):
        try:
            os.mkdir(args.outdir)
        except:
            sys.stderr.write("ERROR\t" + now() + "\tcould not create output directory: " + args.outdir + "\n")
            sys.exit(1)

    for seq in args.seqs:
        sys.stderr.write("INFO\t" + now() + "\tprocessing " + seq + "\n")
        fastq = None
        temp  = False

        if seq.endswith('.bam'):
            fastq = bamtofastq(seq, args.outdir, args.samtofastq, threads=int(args.threads))
            temp  = True

        elif seq.endswith('.cram'):
            assert args.genome is not None and os.path.exists(args.genome), "must specify a valid reference genome for CRAM mode\n"
            fastq = cramtofastq(seq, args.outdir, args.cramjar, args.genome, threads=int(args.threads))
            temp  = True

        elif seq.endswith('.fastq') or seq.endswith('.fastq.gz'):
            fastq = seq

        else:
            sys.stderr.write("ERROR\t" + now() + "\tunrecognized file format (extension != .bam or .cram or .fastq\n")
            sys.exit(1)

        assert fastq.endswith('.fastq')

        sys.stderr.write("INFO\t" + now() + "\taligning fastq " + fastq + " to TE reference " + args.ref + "\n")
        bam = bwamem(fastq, args.ref, args.outdir, threads=args.threads)

        if not args.keepfastq and temp:
            os.remove(fastq)

        fqout = args.outdir + '/' + os.path.basename(args.fqout)
        if not fqout.endswith('.gz'):
            fqout += '.gz'
        parsebam(bam, fqout, minmatch=float(args.minmatch))
        os.remove(bam)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='find relevant reads from sequence data')
    parser.add_argument(metavar='<input file bam/cram/fastq>', dest='seqs', nargs=1, help='input files')
    parser.add_argument('--samtofastq', default=None, help='path to picard SamToFastq.jar')
    parser.add_argument('--cramjar', default=None, help='path to cramtools .jar')
    parser.add_argument('-o', '--outdir', default='.', help='path to output directory')
    parser.add_argument('-f', '--fqout', required=True, help='output FASTQ file')
    parser.add_argument('-r', '--ref', required=True, help='TE reference fasta (needs bwa index and samtools faidx')
    parser.add_argument('-g', '--genome', default=None, help='genome reference fasta (e.g. GRCh37, needs bwa index and samtools faidx')
    parser.add_argument('-t', '--threads', default=1, help='threads for alignment (default = 1)')
    parser.add_argument('-m', '--minmatch', default=0.90, help='minimum percent match to output read (default=0.9)')
    parser.add_argument('--keepfastq', action='store_true', default=False, help='keep temporary fastq file if created (default = False)')
    args = parser.parse_args()
    main(args)