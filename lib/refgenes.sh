#!/bin/sh

wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFaMasked.tar.gz
tar xvzf chromFaMasked.tar.gz
cat *.fa.masked > ucsc.hg19.masked.fa
rm *.masked
samtools faidx ucsc.hg19.masked.fa 
./make_refGene_fa.py ucsc.hg19.masked.fa refGene.txt.gz > refGenes.fa
