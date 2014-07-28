#!/usr/bin/env python

import sys
import argparse 

from pysam import Tabixfile
from re import search
from collections import defaultdict as dd
from collections import Counter


columns = ['Chr',
           'Left_Position',
           'Right_Position',
           'Strand',
           'Strand_Conf',
           'Class',
           'Family',
           'Family_Conf',
           'Found_Ends',
           'Elt_Length',
           'TSD',
           'TSD_Length',
           'DEL',
           'DEL_Length',
           'Left_Snip',
           'Right_Snip',
           'Left_Consensus',
           'Right_Consensus',
           'Support',
           'Libraries',
           'Library_List',
           'Mapscore',
           'Known_Nonref']


def getlibs(invcf):
    rg = {}
    with open(invcf, 'r') as vcf:
        for line in vcf:
            if not line.startswith('#'):
                chrom, pos, id, ref, alt, qual, filter, info, format, sample = line.strip().split('\t')
                rg[sample.split(':')[-1].split(',')] = True
    return rg.keys() 


def main(args):
    print '\t'.join(columns)

    # read consensus
    cons = []
    seq  = []

    nr = Tabixfile(args.nrtabix)

    with open(args.fasta, 'r') as consfa:
        for line in consfa:
            if line.startswith('>'):
                cons.append(line.strip().lstrip('>'))
            else:
                seq.append(line.strip())

    consdict = dict(zip(cons, seq))

    sys.stderr.write("build consensus dict of " + str(len(consdict)) + " sequences\n")

    with open(args.vcf, 'r') as vcf:
        for line in vcf:
            if not line.startswith('#'):
                chrom, pos, id, ref, alt, qual, filter, info, format, sample = line.strip().split('\t') 

                infodict = {}
                for i in info.split(';'):
                    if search('=', i):
                        key, val = i.split('=')
                        infodict[key] = val
                    else:
                        infodict[i] = True

                formatdict = dict([(key, val) for key, val in zip(format.split(':'), sample.split(':'))])

                if filter == 'PASS':

                    start = int(pos)
                    end   = 'NA'
                    if 'END' in infodict:
                        end = int(infodict['END'])
                        if start > end:
                            start, end = end, start

                    alt = alt.split(':')[-1].strip('>')
                    strand = 'NA'
                    strandconf = 0.0

                    if int(formatdict['FWD']) > int(formatdict['REV']):
                        strand = '+'
                        strandconf = float(formatdict['FWD']) / (float(formatdict['FWD']) + float(formatdict['REV']))

                    if int(formatdict['FWD']) < int(formatdict['REV']):
                        strand = '-'
                        strandconf = float(formatdict['REV']) / (float(formatdict['FWD']) + float(formatdict['REV']))

                    # build family info
                    famvec  = []
                    for rawfam in formatdict['TEALIGN'].split(','):
                        fam, count = rawfam.split('|')
                        for _ in range(int(count)):
                            famvec.append(fam)

                    famtotal = len([f for f in famvec if f != 'POLYA'])
                    family, famcount = Counter(famvec).most_common(1)[0]

                    if family == 'POLYA':
                        family, famcount = Counter(famvec).most_common(2)[1]

                    famconf = float(famcount)/float(famtotal)

                    eltends = formatdict['TESIDES']

                    eltlen  = 'NA'
                    if 'TELEN' in infodict:
                        eltlen = infodict['TELEN']

                    support = formatdict['RC']

                    TSD = 'N'
                    DEL = 'N'
                    tlen = 0
                    dlen = 0

                    if 'MECH' in infodict:
                        if infodict['MECH'] == 'TSD':
                            TSD = 'Y'
                            tlen = infodict['TSDLEN']

                        if infodict['MECH'] == 'Deletion':
                            DEL = 'Y'
                            dlen = infodict['TSDLEN']

                    rightsnip = formatdict['POSBREAKSEQ']
                    leftsnip = 'NA'
                    if 'ENDBREAKSEQ' in formatdict:
                        leftsnip = formatdict['ENDBREAKSEQ']

                    leftconskey  = ':'.join((chrom, pos, alt, '1'))
                    rightconskey = ':'.join((chrom, pos, alt, '2'))

                    leftcons  = consdict[leftconskey]
                    rightcons = 'NA'
                    if rightconskey in consdict:
                        rightcons = consdict[rightconskey]

                    liblist  = formatdict['RG']
                    libcount = len(liblist.split(','))
                    mapscore = infodict['MAPSCORE']

                    searchwindow = 100
                    nonref = None
                    if args.nrtabix is not None and chrom in nr.contigs:
                        nonreflist = [elt.split()[3] + '|' + elt.split()[-1] for elt in nr.fetch(chrom, start-searchwindow, end+searchwindow)]
                        if nonreflist:
                            nonref = ','.join(nonreflist)

                    data = [chrom, start, end, strand, strandconf, alt, family, famconf, eltends, eltlen,
                            TSD, tlen, DEL, dlen, rightsnip, leftsnip, leftcons, rightcons, support, 
                            libcount, liblist, mapscore, nonref]

                    assert len(data) == len(columns), "column count: " + str(len(columns)) + " data count: " + str(len(data))

                    print '\t'.join(map(str, data))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='make a human-readable table from tebreak VCF')
    parser.add_argument('-v', '--vcf', required=True, help='Input VCF')
    parser.add_argument('-f', '--fasta', required=True, help='consensus FASTA (not reference fasta!)')
    parser.add_argument('-n', '--nrtabix', default=None, help='known non-ref element tabix')
    args = parser.parse_args()
    main(args)
