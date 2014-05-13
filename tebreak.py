#!/usr/bin/env python

'''
retrosplit.py: Identify TE insertion sites from split reads


contact: Adam Ewing (adam.ewing@mater.uq.edu.au)

'''


import argparse
import subprocess
import gzip
import os
import sys
import pysam
import datetime

from string import maketrans
from collections import Counter
from collections import OrderedDict as od
from uuid import uuid4
from itertools import izip
from re import sub, search
from textwrap import dedent
from shutil import move


class SplitRead:
    ''' gread = genome read, tread = TE read '''
    def __init__(self, gread, tread, tname, chrom, loc, tlen, rescue=False, breakloc=None):
        ''' if read is obtained through hard-clipped rescue, need to set breakloc manually '''
        self.gread     = gread
        self.tread     = tread  
        self.alignloc  = int(loc)
        self.breakloc  = None
        self.teside    = None # '5' or '3' (5' or 3')
        self.chrom     = chrom
        self.tlen      = int(tlen)
        self.hardclips = [] # keep hardclipped reads for later
        self.rescued   = False
        self.embedded  = False # nonref TE is completely embedded in read

        self.original  = None # used for storing non-hardclipped read if rescued

        self.tclass, self.tname = tname.split(':')
        
        # find breakpoint
        if rescue:
            self.breakloc = breakloc
            self.rescued  = True
        else:
            if gread.qstart < gread.rlen - gread.qend:
                self.breakloc  = gread.positions[-1] # breakpoint on right
            else:
                self.breakloc  = gread.positions[0] # breakpoint on left

        self.dist_from_3p = self.tlen - max(self.tread.positions)

        if self.dist_from_3p < 50:
            self.teside = '3'
        else:
            self.teside = '5'

        if self.tclass == 'POLYA':
            self.teside = '3'

        assert self.breakloc is not None
        assert self.teside is not None

    def get_teloc(self):
        ''' return start and end location of alignment within TE '''
        return min(self.tread.positions), max(self.tread.positions)

    def get_teorient(self):
        ''' return orientation of TE alignment relative to TE reference '''
        if self.tread.is_reverse:
            return '-'
        else:
            return '+'

    def get_tematch(self):
        ''' return number of mismatches / aligned length of TE sub-read '''
        nm = [value for (tag, value) in self.tread.tags if tag == 'NM'][0]
        return 1.0 - (float(nm)/float(self.tread.alen))

    def get_gematch(self):
        ''' return number of mismatches / aligned length of TE sub-read '''
        nm = [value for (tag, value) in self.gread.tags if tag == 'NM'][0]
        return 1.0 - (float(nm)/float(self.gread.alen))

    def getbreakseq(self, flank=5):
        ''' get breakend sequence from genome alignments '''

        if self.rescued:
            origseq = self.original.seq.upper()
            assert origseq is not None

            rstart = origseq.find(self.gread.seq.upper())
            if rstart < 0:
                rstart = origseq.find(rc(self.gread.seq).upper())

            assert rstart >= 0

            rend = rstart + len(self.gread.seq)
            return lc(origseq, rstart, rend)

        else:
            leftseq  = ''
            rightseq = ''

            if self.gread.qstart < self.gread.rlen - self.gread.qend: # right
                leftseq  = self.gread.seq[self.gread.qend-flank:self.gread.qend].upper()
                rightseq = self.gread.seq[self.gread.qend:self.gread.qend+flank].lower()
            else: # left
                leftseq  = self.gread.seq[self.gread.qstart-flank:self.gread.qstart].lower()
                rightseq = self.gread.seq[self.gread.qstart:self.gread.qstart+flank].upper()

            return leftseq + rightseq

    def __gt__(self, other):
        ''' enables sorting of SplitRead objects '''
        if self.chrom == other.chrom:
            return self.breakloc > other.breakloc
        else:
            return self.chrom > other.chrom

    def __str__(self):
        return ','.join(map(str, ('SplitRead',self.chrom, self.breakloc, self.tname, self.tread.pos)))


def lc(seq, start, end):
    ''' lowercase part of a sequence '''
    assert start < end
    return seq[:start].upper() + seq[start:end].lower() + seq[end:].upper()


class Cluster:
    ''' store and manipulate groups of SplitRead objects '''
    def __init__(self, firstread=None):
        self._splitreads = []
        self._start  = 0
        self._end    = 0
        self._median = 0
        self.chrom   = None

        # data for VCF fields
        self.POS    = self._median
        self.INFO   = od()
        self.FORMAT = od()
        self.FILTER = []
        self.REF    = '.'
        self.ID     = '.'
        self.ALT    = '<INS:ME:FIXME>'
        self.QUAL   = 100

        # data about insertion point
        self.tsd = None
        self.deletion = None
        self.fraction = None

        if firstread is not None:
            self.add_splitread(firstread)

        self.INFO['KNOWN'] = '0'

    def add_splitread(self, sr):
        ''' add a SplitRead and update '''
        self._splitreads.append(sr)
        if self.chrom is None:
            self.chrom = sr.chrom

        assert self.chrom == sr.chrom # clusters can't include > 1 chromosome

        if self._median > 0 and (self._median - sr.breakloc) > 1000:
            print "WARNING: Splitread", str(sr), "more than 1000 bases from median of cluster."

        ''' update statistics '''
        self._splitreads.sort()
        self._start  = self._splitreads[0].breakloc
        self._end    = self._splitreads[-1].breakloc
        self._median = self._splitreads[len(self)/2].breakloc
        self.POS     = self._median

    def greads(self):
        ''' return list of genome-aligned pysam.AlignedRead objects '''
        return [sr.gread for sr in self._splitreads]

    def treads(self):
        ''' return list of genome-aligned pysam.AlignedRead objects '''
        return [sr.tread for sr in self._splitreads]

    def te_classes(self):
        ''' return distinct classes of TE in cluster '''
        return list(set([read.tclass for read in self._splitreads]))

    def te_names(self):
        ''' return distinct classes of TE in cluster as Counter '''
        return Counter([read.tname for read in self._splitreads])

    def te_sides(self):
        ''' return distinct classes of TE in cluster '''
        return list(set([read.teside for read in self._splitreads]))

    def guess_telen(self):
        ''' approximate TE length based on TE alignments'''
        extrema = self.te_extrema()
        return max(extrema) - min(extrema)

    def te_extrema(self):
        te_positions = []
        [[te_positions.append(tepos) for tepos in sr.tread.positions] for sr in self._splitreads]
        return min(te_positions), max(te_positions)

    def all_breakpoints(self):
        return Counter([read.breakloc for read in self._splitreads])

    def best_breakpoints(self, refbamfn): #TODO left-shift
        ''' Return the most supported breakends. If tied, choose the pair that results in the smallest TSD '''
        be   = []
        mech = 'SingleEnd'

        if len(self.te_sides()) == 1:
            return self.all_breakpoints().most_common(1)[0][0], None, mech

        n = 0
        for breakend, count in self.all_breakpoints().most_common():
            n += 1
            if len(be) < 2:
                be.append([breakend, count])
                if n == 2:
                    mech = self.guess_mechanism(refbamfn, be[0][0], breakend)
            else:
                if count == be[-1][1]: # next best breakend is tied for occurance count
                    newmech = self.guess_mechanism(refbamfn, be[0][0], breakend)
                    # prefer TSD annotations over all others
                    if mech != 'TSD' and newmech == 'TSD':
                        be[-1] = [breakend, count]
                        mech = 'TSD'
                    
                    # ... but if neither or both are TSD, choose the smallest
                    if (newmech == mech == 'TSD') or (mech != 'TSD' and newmech != 'TSD'): 
                        if abs(be[-1][0] - be[0][0]) > abs(breakend - be[0][0]):
                            be[-1] = [breakend, count]
                            mech = newmech

        if len(be) == 1:
            return be[0][0], None, mech

        if len(be) == 2:
            return be[0][0], be[1][0], mech

        raise ValueError('Cluster.best_breakpoints returned more than two breakends, panic!\n')

    def majorbreakseq(self, breakloc, flank=5):
        ''' return most common breakend sequence '''
        gbest = [sr.getbreakseq(flank=flank) for sr in self._splitreads if sr.breakloc == breakloc and not sr.rescued]

        if len(gbest) == 0:
            gbest = [sr.getbreakseq(flank=flank) for sr in self._splitreads if sr.breakloc == breakloc and sr.rescued]
            return Counter(gbest).most_common(1)[0][0] + ',HCLIP'
        else:
            return Counter(gbest).most_common(1)[0][0] + ',SCLIP'

    def guess_mechanism(self, refbamfn, bstart, bend):
        ''' check depth manually (only way to get ALL reads), decide whether it's likely a TSD, Deletion, or other '''
        if bstart > bend:
            bstart, bend = bend, bstart

        if len(self.te_sides()) == 1:
            return "SingleEnd"

        if bstart == bend:
            return "EndoMinus"

        window = 10
        bam    = pysam.Samfile(refbamfn, 'rb')
        depth  = od([(pos, 0) for pos in range(bstart-window, bend+window)])

        for read in bam.fetch(self.chrom, bstart-window, bend+window):
            for pos in read.positions:
                if pos in depth:
                    depth[pos] += 1

        cov = depth.values()

        meancov_left  = float(sum(cov[:window]))/float(window)
        meancov_right = float(sum(cov[-window:]))/float(window)
        meancov_break = float(sum(cov[window:-window]))/float(len(cov[window:-window]))

        if meancov_break > meancov_left and meancov_break > meancov_right:
            return "TSD"

        if meancov_break < meancov_left and meancov_break < meancov_right:
            return "Deletion"

        return "Unknown"

    def subcluster_by_class(self, teclass):
        ''' return a new cluster contaning only TEs of teclass '''
        new = Cluster()
        [new.add_splitread(read) for read in self._splitreads if read.tclass == teclass or read.tclass == 'POLYA']
        return new

    def subcluster_by_breakend(self, breakends):
        ''' return a new cluster containing only reads with breakpoints in passed list '''
        new = Cluster()
        [new.add_splitread(read) for read in self._splitreads if read.breakloc in breakends]
        return new

    def mean_genome_qual(self):
        quals = [sr.gread.mapq for sr in self._splitreads]
        return float(sum(quals))/float(len(quals))

    def get_hardclips(self):
        hclist = []
        for sr in self._splitreads:
            for hc in sr.hardclips:
                hclist.append(hc)
        return hclist

    def unclipfrac(self, refbamfn):
        ''' track fraction of reads over cluster region that are unclipped '''
        used   = {}
        unclip = 0
        total  = 0

        bam = pysam.Samfile(refbamfn, 'rb')
        for read in bam.fetch(self.chrom, self._start, self._end+1):
            uid = ':'.join(map(str, (read.tid,read.pos,read.rlen,read.is_reverse)))
            if uid not in used:
                if read.alen == read.rlen and not search('H', read.cigarstring):
                    unclip += 1
                total += 1
            used[uid] = True

        return float(unclip)/float(total)

    def checksnps(self, dbsnpfn, refbamfn):
        ''' dbsnpfn should be a tabix-indexed version of a dbSNP VCF '''
        dbsnp = pysam.Tabixfile(dbsnpfn)
        bam   = pysam.Samfile(refbamfn, 'rb')

        foundsnps  = []
        start, end = self.find_extrema() 
        
        snps = od()
        if self.chrom in dbsnp.contigs:
            for rec in dbsnp.fetch(self.chrom, start, end):
                snps[int(rec.strip().split()[1])] = rec

        for read in bam.fetch(self.chrom, start, end):
            for rpos, gpos in read.aligned_pairs:
                if gpos is not None and rpos is not None and int(gpos+1) in snps:
                    snp_chrom, snp_pos, snp_id, snp_ref, snp_alt = snps[int(gpos+1)].split()[:5]
                    if rpos is None:
                        print "DEBUG: rpos is None somewhere:", read.aligned_pairs
                    clusterbase = read.seq[rpos+read.qstart].upper() # rpos offsets to first _aligned_ base
                    if len(snp_ref) == len(snp_alt) == 1 and clusterbase == snp_alt.upper():
                        foundsnps.append(':'.join((snp_id, clusterbase, 'ALT')))
                    if len(snp_ref) == len(snp_alt) == 1 and clusterbase == snp_ref.upper():
                        foundsnps.append(':'.join((snp_id, clusterbase, 'REF')))

        bam.close()
        reportsnps = [snp for snp, count in Counter(foundsnps).iteritems() if count > 1]
        return reportsnps

    def find_extrema(self):
        ''' return leftmost and rightmost aligned positions in cluster vs. reference '''
        positions = []
        positions += [pos for sr in self._splitreads for pos in sr.gread.positions]
        return min(positions), max(positions)

    def __str__(self):
        ''' convert cluster into VCF record '''
        output = [self.chrom, str(self.POS), self.ID, self.REF, self.ALT, str(self.QUAL)]
        output.append(';'.join(self.FILTER))
        output.append(';'.join([key + '=' + str(val) for key,val in self.INFO.iteritems()]))
        output.append(';'.join([key for key,val in self.FORMAT.iteritems()]))
        output.append(';'.join([str(val) for key,val in self.FORMAT.iteritems()]))
        return '\t'.join(output)

    def __len__(self):
        return len(self._splitreads)


class AlignedColumn:
    ''' used by MSA class to store aligned bases '''
    def __init__(self):
        self.bases = od() # species name --> base
        self.annotations = od() # metadata

    def gap(self):
        if '-' in self.bases.values():
            return True
        return False

    def subst(self):
        if self.gap():
            return False

        if len(set(self.bases.values())) > 1:
            return True
        return False

    def cons(self):
        ''' consensus prefers bases over gaps '''
        for base, count in Counter(map(str.upper, self.bases.values())).most_common():
            if base != '-':
                return base

    def __str__(self):
        return str(self.bases)


class MSA:
    ''' multiple sequence alignment class '''
    def __init__(self, infile=None):
        self.columns = []
        self.ids     = []
        self.seqs    = od()

        if infile is not None:
            self.readFastaMSA(infile)

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
        bases = [column.cons() for column in self.columns]
        if bases is not None:
            return ''.join(bases)
        else:
            sys.stderr.write("ERROR\t" + now() + "\tNone found in consensus sequence\n")


def now():
    return str(datetime.datetime.now())


def print_vcfheader(fh, samplename):
    fh.write(sub('FIXME',samplename,dedent("""\
    ##fileformat=VCFv4.1
    ##phasing=none
    ##INDIVIDUAL=FIXME
    ##SAMPLE=<ID=FIXME,Individual="FIXME",Description="sample name">
    ##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
    ##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
    ##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
    ##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
    ##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic mutation in primary">
    ##INFO=<ID=END,Number=1,Type=Integer,Description="Position of second breakend (same as first if not known)">
    ##INFO=<ID=TESIDES,Number=.,Type=Integer,Description="Note if 3prime and 5prime ends are detected">
    ##INFO=<ID=MECH,Number=1,Type=String,Description="Hints about the insertion mechanism (TSD, Deletion, EndoMinus, etc.)">
    ##INFO=<ID=KNOWN,Number=1,Type=Integer,Description="1=Known from a previous study (in whitelist)">
    ##FORMAT=<ID=BREAKS,Number=.,Type=String,Description="Positions:Counts of all breakends detected">
    ##FORMAT=<ID=POSBREAKSEQ,Number=1,Type=String,Description="Sequence of the POS breakend">
    ##FORMAT=<ID=ENDBREAKSEQ,Number=1,Type=String,Description="Sequence of the INFO.END breakend">
    ##FORMAT=<ID=TEALIGN,Number=.,Type=String,Description="Retroelement subfamilies (or POLYA) with alignments">
    ##FORMAT=<ID=RC,Number=1,Type=Float,Description="Read Count">
    ##FORMAT=<ID=UCF,Number=1,Type=Float,Description="Unclipped Fraction">
    ##ALT=<ID=DEL,Description="Deletion">
    ##ALT=<ID=INS,Description="Insertion">
    ##ALT=<ID=INS:ME:ALU,Description="Insertion of ALU element">
    ##ALT=<ID=INS:ME:L1,Description="Insertion of L1 element">
    ##ALT=<ID=INS:ME:SVA,Description="Insertion of SVA element">
    ##ALT=<ID=INS:ME:POLYA,Description="Insertion of POLYA sequence">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tFIXME"""))+'\n')


def rc(dna):
    ''' reverse complement '''
    complements = maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
    return dna.translate(complements)[::-1]


def read_fasta(infa):
    ''' potato '''
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


def bwamem(fq, ref, threads=1, width=150, sortmem=2000000000, uid=None):
    ''' FIXME: add parameters to commandline '''
    fqroot = sub('.extendedFrags.fastq.gz$', '', fq)
    fqroot = sub('.fq$', '', fqroot)
    if uid is not None:
        fqroot = uid

    if str(sortmem).rstrip('Gg') != str(sortmem):
        sortmem = int(str(sortmem).rstrip('Gg')) * 1000000000

    if str(sortmem).rstrip('Mm') != str(sortmem):
        sortmem = int(str(sortmem).rstrip('Gg')) * 1000000

    sortmem = sortmem/int(threads) # avoid PBS killing my jobs

    sam_out  = '.'.join((fqroot, 'sam'))
    bam_out  = '.'.join((fqroot, 'bam'))
    sort_out = '.'.join((fqroot, 'sorted'))

    sam_cmd  = ['bwa', 'mem', '-t', str(threads), '-M', ref, fq]
    bam_cmd  = ['samtools', 'view', '-bt', ref + '.bai', '-o', bam_out, sam_out]
    sort_cmd = ['samtools', 'sort', '-@', str(threads), '-m', str(sortmem), bam_out, sort_out]
    idx_cmd  = ['samtools', 'index', sort_out + '.bam']

    sys.stdout.write("INFO\t" + now() + "\trunning bwa-mem: " + ' '.join(sam_cmd) + "\n")

    with open(sam_out, 'w') as sam:
        p = subprocess.Popen(sam_cmd, stdout=subprocess.PIPE)
        for line in p.stdout:
            sam.write(line)

    sys.stdout.write("INFO\t" + now() + "\twriting " + sam_out + " to BAM...\n")
    subprocess.call(bam_cmd)
    os.remove(sam_out)

    sys.stdout.write("INFO\t" + now() + "\tsorting output: " + ' '.join(sort_cmd) + "\n")
    subprocess.call(sort_cmd)
    os.remove(bam_out)

    sys.stdout.write("INFO\t" + now() + "\tindexing: " + ' '.join(idx_cmd) + "\n")
    subprocess.call(idx_cmd)

    return sort_out + '.bam'


def flash_wrapper(fq1, fq2, max_overlap, threads, uid=None):
    ''' wrapper for FLASH (http://ccb.jhu.edu/software/FLASH/) '''
    out = 'RCTMP-' + str(uuid4())
    if uid is not None:
         out = uid

    args = ['flash', '-d', os.path.dirname(out), '-o', os.path.basename(out), '-M', str(max_overlap), '-z', '-t', str(threads), fq1, fq2]
    print "INFO\t" + now() + "\tcalling FLASH:", args

    subprocess.call(args)

    return out + '.extendedFrags.fastq.gz'


def bamrec_to_fastq(read, diffseq=None, diffqual=None, diffname=None):
    ''' convert BAM record to FASTQ record: can substitute name/seq/qual using diffseq/diffqual/diffname '''
    name = read.qname
    if diffname is not None:
        name = diffname

    qual = read.qual
    if diffqual is not None:
        qual = diffqual

    seq  = read.seq
    if diffseq is not None:
        seq = diffseq

    if read.is_reverse:
        qual = qual[::-1]
        seq  = rc(seq)

    return '\n'.join(('@'+name,seq,"+",qual))


def fetch_clipped_reads(inbamfn, minclip=50, maxaltclip=2): # TODO PARAMS
    ''' FIXME: add parameters to commandline '''
    assert minclip > maxaltclip

    outfqfn = sub('.bam$', '.clipped.fq', inbamfn)
    inbam   = pysam.Samfile(inbamfn, 'rb')

    with open(outfqfn, 'w') as outfq:
        for read in inbam.fetch():
            uid = ':'.join(map(str, (read.tid,read.pos,read.rlen,read.is_reverse)))
            unmapseq = None
            unmapqua = None

            if read.qual is None:
                sys.stderr.write("ERROR\t" + now() + "\tread with no quality score:\n")
                sys.stderr.write(str(read) + "\n")
                sys.stderr.write("ERROR\t" + now() + "\tpossible FASTA alignment instead of FASTQ - please re-do alignments from FASTQ\n")
                sys.exit(1)

            if read.rlen - read.alen >= int(minclip): # 'soft' clipped?

                # length of 'minor' clip
                altclip = min(read.qstart, read.rlen-read.qend)

                if altclip <= maxaltclip:
                    # get unmapped part
                    if read.qstart <= maxaltclip:
                        # (align) AAAAAAAAAAAAAAAAAA
                        # (read)  RRRRRRRRRRRRRRRRRRRRRRRRRRRR
                        unmapseq = read.seq[read.qend:]
                        unmapqua = read.qual[read.qend:]

                    else:
                        # (align)           AAAAAAAAAAAAAAAAAA
                        # (read)  RRRRRRRRRRRRRRRRRRRRRRRRRRRR
                        unmapseq = read.seq[:read.qstart]
                        unmapqua = read.qual[:read.qstart]

                else:
                    #TODO, consider this case, important for ALU insertions w/ read length > 300 bp and very short L1s
                    # (align)           AAAAAAAAAAAAAAA
                    # (read)   RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
                    pass

            if not read.is_unmapped and unmapseq is not None and not read.is_duplicate:
                infoname = ':'.join((read.qname, str(inbam.getrname(read.tid)), str(read.pos))) # help find read again
                outfq.write(bamrec_to_fastq(read, diffseq=unmapseq, diffqual=unmapqua, diffname=infoname) + '\n')

        inbam.close()

    return outfqfn


def build_te_splitreads(refbamfn, tebamfn, teref, min_te_match = 0.9, min_ge_match = 0.98): # TODO PARAM
    ''' g* --> locations on genome; t* --> locations on TE '''
    refbam = pysam.Samfile(refbamfn, 'rb')
    tebam  = pysam.Samfile(tebamfn, 'rb')

    ''' get sorted list of split reads '''
    splitreads = []
    for tread in tebam.fetch():
        if not tread.is_unmapped and not tread.is_secondary:
            tname = tebam.getrname(tread.tid)
            gname = ':'.join(tread.qname.split(':')[:-2])
            gchrom, gloc = tread.qname.split(':')[-2:]

            gloc = int(gloc)
            sr = None
            hcreads = []
            for gread in refbam.fetch(reference=gchrom, start=gloc, end=gloc+1):
                if (gread.qname, gread.pos, refbam.getrname(gread.tid)) == (gname, gloc, gchrom):
                    tlen = len(teref[tname])
                    sr = SplitRead(gread, tread, tname, gchrom, gloc, tlen)
                    if sr.get_tematch() >= min_te_match and sr.get_gematch() >= min_ge_match:
                        splitreads.append(sr)
                if search('H', gread.cigarstring) and gread.is_secondary:
                    hcreads.append(gread)
            if sr is not None:
                [sr.hardclips.append(hc) for hc in hcreads]

    refbam.close()
    tebam.close()
    splitreads.sort()

    return splitreads


def build_te_clusters(splitreads, searchdist=100): # TODO PARAM
    ''' cluster SplitRead objects into Cluster objects and return a list of them '''
    clusters  = []

    for sr in splitreads:
        if len(clusters) == 0:
            clusters.append(Cluster(sr))

        elif clusters[-1].chrom != sr.chrom:
            clusters.append(Cluster(sr))

        else:
            if abs(clusters[-1]._median - sr.breakloc) > searchdist:
                clusters.append(Cluster(sr))

            else:
                clusters[-1].add_splitread(sr)

    return clusters


def rescue_hardclips(clusters, refbamfn, telib, width=150, threads=1):
    ''' hard-clipped reads are often mostly TE-derived, hunt down the original read and find breakpoints '''
    hc_clusters = {}
    hc_reads    = {}
    hc_original = {}
    cn = 0

    for cluster in clusters:
        for hcread in cluster.get_hardclips():
            hc_clusters[hcread.qname] = cn
            hc_reads   [hcread.qname] = hcread
        cn += 1

    bam  = pysam.Samfile(refbamfn, 'rb')
    fqfn = sub('.bam$', '.rescue.fq', refbamfn)

    with open(fqfn, 'w') as fqout:
        for read in bam.fetch(until_eof=True):
            if not read.is_secondary:
                if read.qname in hc_clusters:
                    # names match, hardclipped read should be a subsequence of read.seq
                    cliploc = read.seq.find(hc_reads[read.qname].seq)

                    if cliploc < 0:
                        cliploc = read.seq.find(rc(hc_reads[read.qname].seq))

                    assert cliploc >= 0
                    hc_original[read.qname] = read
                    
                    fqrec = bamrec_to_fastq(read, diffseq=read.seq[:cliploc], diffqual=read.qual[:cliploc])
                    fqout.write(fqrec + '\n')

    print "INFO\t" + now() + "\tealigning hardclips to TE library..."
    bamout = bwamem(fqfn, telib, width=width, threads=threads)

    tebam = pysam.Samfile(bamout, 'rb')

    te_length = dict(zip(tebam.references, tebam.lengths))

    for teread in tebam.fetch():
        if not teread.is_secondary:
            breakloc = hc_reads[teread.qname].pos
            # break location is at right end of read
            if hc_reads[teread.qname].cigarstring.endswith('H'):
                breakloc = hc_reads[teread.qname].positions[-1]

            gread  = hc_reads[teread.qname]
            tname  = tebam.getrname(teread.tid)
            gchrom = bam.getrname(gread.tid)
            gloc   = gread.pos
            tlen   = te_length[tname]

            sr = SplitRead(gread, teread, tname, gchrom, gloc, tlen, rescue=True, breakloc=breakloc)
            sr.original = hc_original[teread.qname]
            cn = hc_clusters[teread.qname]
            clusters[cn].add_splitread(sr)

    tebam.close()
    bam.close()

    return clusters


def filter_clusters(clusters, active_elts, refbamfn, minsize=4, bothends=False, maskfile=None, whitelistfile=None, minctglen=1e7, minq=1, unclip=1.0):
    ''' return only clusters meeting cutoffs '''
    ''' active_elts is a list of active transposable element names (e.g. L1Hs) '''
    ''' refbamfn is the filename of the reference-aligned BAM '''
    ''' maskfile / whitelist should be a tabix-indexed BED file of intervals to ignore (e.g. L1Hs reference locations) '''
    filtered = []
    active_elts = Counter(active_elts)

    mask = None
    if maskfile is not None:
        mask = pysam.Tabixfile(maskfile)

    whitelist = None
    if whitelistfile is not None:
        whitelist = pysam.Tabixfile(whitelistfile)

    bam = pysam.Samfile(refbamfn, 'rb')
    chromlen = dict(zip(bam.references, bam.lengths))

    for cluster in clusters:
        for teclass in cluster.te_classes(): # can consider peaks with more than one TE class
            subcluster = cluster.subcluster_by_class(teclass)
            tclasses   = Counter([read.tclass for read in subcluster._splitreads]).most_common()

            # prefer to report element type even if POLYA is the dominant annotation
            tclass = tclasses[0][0]
            if tclass == 'POLYA' and len(tclasses) > 1:
                tclass = tclasses[1][0]

            subcluster.ALT  = sub('FIXME', tclass, subcluster.ALT)
            subcluster.QUAL = int(subcluster.mean_genome_qual())

            reject = True
            
            # at least one alignment in a cluster has to be from an active element
            for tn in subcluster.te_names():
                if tn in active_elts:
                    reject = False
                    
            if reject:
                subcluster.FILTER.append('tetype')

            if subcluster.QUAL <= minq:
                reject = True
                subcluster.FILTER.append('avgmapq')

            # clusters have to be at least some minimum size
            if len(subcluster) < minsize:
                reject = True
                subcluster.FILTER.append('csize')

            # filter out clusters without both TE ends represented if user wants this
            if bothends and len(cluster.te_sides()) < 2:
                reject = True
                subcluster.FILTER.append('bothends')

            # apply position-based masking
            if mask is not None and subcluster.chrom in mask.contigs:
                if len(list(mask.fetch(subcluster.chrom, subcluster._start, subcluster._end+1))) > 0:
                    reject = True
                    subcluster.FILTER.append('masked')

            # filter out contigs below specified size cutoff
            if chromlen[subcluster.chrom] < minctglen:
                reject = True
                subcluster.FILTER.append('ctgsize')

            # save the more time-consuming cutoffs for last... check unclipped read mappings:
            if not reject and unclip < 1.0:
                subcluster.FORMAT['UCF'] = cluster.unclipfrac(refbamfn)
                if subcluster.FORMAT['UCF'] > unclip:
                    reject = True
                    subcluster.FILTER.append('unclipfrac')

            if 'masked' not in subcluster.FILTER and whitelist is not None and subcluster.chrom in whitelist.contigs:
                tclasses = [wl.strip().split()[-1] for wl in whitelist.fetch(subcluster.chrom, subcluster._start, subcluster._end+1)]
                if tclass in tclasses:
                    reject = False
                    subcluster.INFO['KNOWN'] = '1'
                    subcluster.FILTER = []

            if not reject:
                subcluster.FILTER.append('PASS')

            filtered.append(subcluster)

    bam.close()
    return filtered


def annotate(clusters, reffa, refbamfn, allclusters=False, dbsnp=None, minclip=10):
    ''' populate VCF INFO field with information about the insertion '''

    ref = pysam.Fastafile(reffa)

    for cluster in clusters:
        cluster.INFO['SVTYPE'] = 'INS'

        breaklocs = []
        for breakloc, count in cluster.all_breakpoints().iteritems():
            breaklocs.append(str(breakloc) + ':' + str(count))

        cluster.FORMAT['BREAKS']  = ','.join(breaklocs)
        cluster.FORMAT['TESIDES'] = ','.join(cluster.te_sides())

        tetypes = []
        for tetype, count in cluster.te_names().iteritems():
            tetypes.append(tetype + ':' + str(count))
        cluster.FORMAT['TEALIGN'] = ','.join(tetypes)
        cluster.FORMAT['TEMINPOS'], cluster.FORMAT['TEMAXPOS'] = map(str, cluster.te_extrema())

        if cluster.FILTER[0] == 'PASS' or allclusters: 
            leftbreak, rightbreak, mech = cluster.best_breakpoints(refbamfn)
            if rightbreak is not None and leftbreak > rightbreak:
                leftbreak, rightbreak = rightbreak, leftbreak

            cluster.POS = leftbreak
            cluster.INFO['END'] = leftbreak

            if rightbreak is not None:
                cluster.INFO['END']   = rightbreak
                cluster.INFO['MECH']  = mech
                cluster.INFO['TELEN'] = cluster.guess_telen()

            if dbsnp is not None:
                snps = cluster.checksnps(dbsnp, refbamfn)
                if snps:
                    cluster.INFO['SNPS'] = ','.join(snps)

        refbase = ref.fetch(cluster.chrom, cluster.POS, cluster.POS+1)
        if refbase != '':
            cluster.REF = refbase
        else:
            cluster.REF = '.'

        cluster.FORMAT['POSBREAKSEQ'] = cluster.majorbreakseq(cluster.POS, flank=int(minclip))

        if 'END' in cluster.INFO and cluster.INFO['END'] != cluster.POS:
            cluster.FORMAT['ENDBREAKSEQ'] = cluster.majorbreakseq(cluster.INFO['END'], flank=int(minclip))

        cluster.FORMAT['RC'] = len(cluster)

    return clusters


def consensus(cluster, breakend):
    ''' create consensus sequence for a breakend '''
    subcluster = cluster.subcluster_by_breakend([breakend])
    greads = [sr.gread for sr in subcluster._splitreads]

    # don't try to make a consensus of just one read...
    if len(greads) == 1:
        return greads[0].seq

    tmpfa = str(uuid4()) + '.fa'
    with open(tmpfa, 'w') as fa:
        [fa.write('>' + gread.qname + '\n' + gread.seq + '\n') for gread in greads]

    alnfa = mafft(tmpfa)
    msa = MSA(alnfa)

    os.remove(tmpfa)
    os.remove(alnfa)

    return msa.consensus()


def mafft(infafn, iterations=100, threads=1):
    ''' use MAFFT to create MSA '''

    outfafn = str(uuid4()) + '.aln.fa'

    args = ['mafft', '--localpair', '--maxiterate', str(iterations), '--thread', str(threads), infafn]

    FNULL = open(os.devnull, 'w')

    with open(outfafn, 'w') as outfa:
        p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=FNULL)
        for line in p.stdout:
            outfa.write(line)

    return outfafn


def bamtofq(inbam, outfq):
    assert inbam.endswith('.bam')
    with gzip.open(outfq, 'wb') as fq:
        bam = pysam.Samfile(inbam, 'rb')
        for read in bam.fetch(until_eof=True):
            seq  = read.seq
            qual = read.qual
            salt = str(uuid4()).split('-')[0]

            if read.is_reverse:
                seq  = rc(seq)
                qual = qual[::-1]

            fq.write('\n'.join(('@' + read.qname + ':' + salt, seq, '+', qual)) + '\n')
    
    bam.close()
    return outfq


def vcfoutput(clusters, outfile, samplename):
    ''' output results in VCF format '''
    with open(outfile, 'w') as out:
        print_vcfheader(out, samplename)
        for cluster in clusters:
            out.write(str(cluster) + '\n')


def bamoutput(clusters, refbamfn, tebamfn, prefix, passonly=False):
    refbam = pysam.Samfile(refbamfn, 'rb')
    tebam  = pysam.Samfile(tebamfn, 'rb')
    refout = pysam.Samfile(prefix + ".ref.bam", 'wb', template=refbam)
    teout  = pysam.Samfile(prefix + ".te.bam", 'wb', template=tebam)
    if passonly:
        refout = pysam.Samfile(prefix + ".pass.ref.bam", 'wb', template=refbam)
        teout  = pysam.Samfile(prefix + ".pass.te.bam", 'wb', template=tebam)

    for cluster in clusters:
        if passonly and cluster.FILTER[0] == 'PASS':
            [refout.write(read) for read in cluster.greads()]
            [teout.write(read) for read in cluster.treads()]
        else:
            [refout.write(read) for read in cluster.greads()]
            [teout.write(read) for read in cluster.treads()]

    refbam.close()
    tebam.close()
    refout.close()
    teout.close()


def consensus_fasta(clusters, outfile, passonly=True):
    with open(outfile, 'w') as cons:
        for cluster in clusters:
            if passonly and cluster.FILTER[0] == 'PASS':
                outname = cluster.chrom + ':' + str(cluster.POS)
                cons.write('>' + outname + ':1\n' + consensus(cluster, cluster.POS) + '\n')
                if cluster.INFO.get('END') and cluster.POS != cluster.INFO['END']:
                    cons.write('>' + outname + ':2\n' + consensus(cluster, cluster.INFO['END']) + '\n')


def markdups(inbam, picard):
    ''' call MarkDuplicates.jar from Picard (requires java) '''
    md = picard + '/MarkDuplicates.jar'
    assert os.path.exists(md)

    outtmp  = str(uuid4()) + '.mdtmp.bam'
    metrics = inbam + '.markdups.metrics.txt'
    args    = ['java', '-Xmx4g', '-jar', md, 'I=' + inbam, 'O=' + outtmp, 'M=' + metrics]
    subprocess.call(args)

    assert os.path.exists(metrics)

    os.remove(inbam)
    move(outtmp, inbam)

    print "INFO:\t" + now() +"\trebuilding index for", inbam
    subprocess.call(['samtools', 'index', inbam])

    return inbam


def main(args):
    sys.stderr.write("INFO\t" + now() + "\tstarting " + sys.argv[0] + " called with args:\n" + ' '.join(sys.argv) + "\n")

    mergefq = None
    basename = args.samplename

    if args.outdir is not None:
        if not os.path.exists(args.outdir):
            try:
                os.mkdir(args.outdir)
            except:
                sys.stderr.write("ERROR: " + now() + " failed to create output directory: " + args.outdir + ", exiting.\n")

        basename = args.outdir + '/' + basename

    if args.pair2 is not None:
        sys.stderr.write("INFO: " + now() + " overlapping paired ends\n")
        mergefq = flash_wrapper(args.pair1, args.pair2, args.maxoverlap, args.threads, uid=basename)

    else:
        if args.pair1.endswith('.bam') and not args.premapped:
            sys.stderr.write("INFO: " + now() + " converting BAM " + args.pair1 + " to FASTQ " + basename + ".fq.gz\n")
            mergefq = bamtofq(args.pair1, basename + '.fq.gz')
        else: # assume fastq
            mergefq = args.pair1

    refbamfn = None
    if args.premapped:
        if args.pair1.endswith('.bam'):
            sys.stderr.write("INFO: " + now() + " using pre-mapped BAM: " + args.pair1 + "\n")
            refbamfn = args.pair1
        else:
            sys.stderr.write("ERROR: " + now() + " flag to use premapped bam (--premapped) called but " + args.pair1 + " is not a .bam file\n")
            sys.exit(1)
    else:
        sys.stderr.write("INFO: " + now() + " mapping fastq " + mergefq + " to genome " + args.ref + " using " +  str(args.threads) + " threads\n")
        refbamfn = bwamem(mergefq, args.ref, width=int(args.width), threads=args.threads, uid=basename, sortmem=args.sortmem)

    assert refbamfn is not None

    sys.stderr.write("INFO: " + now() + " marking likely PCR duplicates\n")
    markdups(refbamfn, args.picard)

    sys.stderr.write("INFO: " + now() + " finding clipped reads from genome alignment\n")
    clipfastq = fetch_clipped_reads(refbamfn, minclip=int(args.minclip))
    assert clipfastq

    sys.stderr.write("INFO: " + now() + " realigning clipped ends to TE reference library\n")
    tebamfn = bwamem(clipfastq, args.telib, width=int(args.width), threads=args.threads, sortmem=args.sortmem)
    assert tebamfn 

    sys.stderr.write("INFO: " + now() + " identifying usable split reads from alignments\n")
    splitreads = build_te_splitreads(refbamfn, tebamfn, read_fasta(args.telib))
    assert splitreads

    sys.stderr.write("INFO: " + now() + " clustering split reads on genome coordinates\n")
    clusters = build_te_clusters(splitreads)
    sys.stderr.write("INFO: " + now() + " cluster count: " + str(len(clusters)) + "\n")

    sys.stderr.write("INFO: " + now() + " further investigation of hard-clipped reads in breakend regions\n")
    clusters = rescue_hardclips(clusters, refbamfn, args.telib, width=int(args.width), threads=args.threads)

    sys.stderr.write("INFO: " + now() + " filtering clusters\n")
    clusters = filter_clusters(clusters, args.active.split(','), refbamfn, 
                               minsize=int(args.mincluster), maskfile=args.mask,
                               whitelistfile=args.whitelist, unclip=float(args.unclip))
    sys.stderr.write("INFO: " + now() + " passing cluster count: " + str(len([c for c in clusters if c.FILTER[0] == 'PASS'])) + "\n")

    sys.stderr.write("INFO: " + now() + " annotating clusters\n")
    clusters = annotate(clusters, args.ref, refbamfn, allclusters=args.processfiltered, dbsnp=args.snps, minclip=args.minclip)
    assert clusters

    sys.stderr.write("INFO: " + now() + " writing VCF\n")
    vcfoutput(clusters, args.outvcf, args.samplename)

    if args.clusterbam:
        sys.stderr.write("INFO: " + now() + " writing BAMs\n")
        bamoutput(clusters, refbamfn, tebamfn, basename, passonly=True)

    if args.consensus is not None:
        sys.stderr.write("INFO: " + now() + " compiling consensus sequences for breakends and outputting to " + args.consensus + "\n")
        consensus_fasta(clusters, args.consensus)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyse RC-seq data')
    parser.add_argument('-1', dest='pair1', required=True, 
                        help='fastq(.gz) containing first end reads or BAM to process and re-map')
    parser.add_argument('-2', dest='pair2', default=None,
                        help='fastq(.gz) containing second end reads (optional, second read is assumed to overlap the first and will be joined via FLASH)')

    parser.add_argument('-r', '--ref', dest='ref', required=True, 
                        help='reference genome for bwa-mem, also expects .fai index (samtools faidx ref.fa)')
    parser.add_argument('-l', '--telib', dest='telib', required=True, 
                        help='TE library (BWA indexed FASTA), seq names must be CLASS:NAME')
    parser.add_argument('-v', '--outvcf', dest='outvcf', required=True, 
                        help='output VCF file')
    parser.add_argument('-o', '--outdir', dest='outdir', default=None, 
                        help='output directory')
    parser.add_argument('-p', '--picard', dest='picard', required=True,
                        help='Picard install directory, needed for MarkDuplicates.jar')

    parser.add_argument('-a', '--active', dest='active', default='L1Hs',
                        help='Comma-delimited list of relevant (i.e. active) subfamilies to target (default=L1Hs)')

    parser.add_argument('-n', '--samplename', dest='samplename', default=str(uuid4()), 
                        help='unique sample name (default = generated UUID4)')
    parser.add_argument('-m', '--mask', dest='mask', default=None, 
                        help='genome coordinate mask (recommended!!) - expects tabix-indexed BED-3')
    parser.add_argument('-w', '--whitelist', dest='whitelist', default=None,
                        help='PASS any non-mask insertions found tabix-indexed BED-3 plus a fourth column corresponding to L1/ALU/SVA')
    parser.add_argument('-s', '--snps', dest='snps', default=None,
                        help='dbSNP VCF (tabix-indexed) to link SNPs with insertions')
    parser.add_argument('-t', '--threads', dest='threads', default=1, 
                        help='number of threads (default = 1)')
    parser.add_argument('--sortmem', dest='sortmem', default='8G',
                        help='amount of memory for sorting (default = 8G)')

    parser.add_argument('--width', dest='width', default=150,
                         help='bandwidth parameter for bwa-mem: determines max size of indels in reads (see bwa docs, default=150)')
    parser.add_argument('--max-overlap', dest='maxoverlap', default=100, 
                        help='Maximum overlap used for joining paired reads with FLASH')
    parser.add_argument('--mincluster', dest='mincluster', default=4, 
                        help='minimum number of reads in a cluster')
    parser.add_argument('--minclip', dest='minclip', default=50, 
                        help='minimum clipped bases for adding a read to a cluster (default = 50)')
    parser.add_argument('--minq', dest='minq', default=1, 
                        help='minimum mean mapping quality per cluster (default = 1)')
    parser.add_argument('--unclipfrac', dest='unclip', default=1.0, 
                        help='maximum fraction of unclipped reads in cluster region (default = 1.0)')
    parser.add_argument('--consensus', dest='consensus', default=None,
                        help='build consensus sequences from breakends and output as FASTA to specified file')

    parser.add_argument('--processfiltered', action='store_true', default=False, 
                        help='perform post-processing steps on all clusters, even filtered ones')
    parser.add_argument('--premapped',  action='store_true', default=False, 
                        help='use BAM specified by -1 (must be .bam) directly instead of remapping')
    parser.add_argument('--clusterbam',  action='store_true', default=False, 
                        help='output genome and TE clusters to BAMs')

    args = parser.parse_args()
    main(args)


