#!/bin/sh

git clone https://github.com/samtools/tabix.git
cd tabix
make
cd ..

if [ ! -e tabix/bgzip ]
then
    echo "tabix/bgzip did not build properly"
    exit 1
fi

if [ ! -e tabix/tabix ]
then
    echo "tabix/tabix did not build properly"
    exit 1 
fi

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign100mer.bigWig

wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigToWig
chmod +x bigWigToWig

./bigWigToWig wgEncodeCrgMapabilityAlign100mer.bigWig wgEncodeCrgMapabilityAlign100mer.wig

if [ ! -e wgEncodeCrgMapabilityAlign100mer.wig ]
then
    echo "could not create wig file"
    exit 1
fi

grep -v '^#' wgEncodeCrgMapabilityAlign100mer.wig | sed -e s'/^chr//' > wgEncodeCrgMapabilityAlign100mer.bed

tabix/bgzip wgEncodeCrgMapabilityAlign100mer.bed
tabix/tabix -s 1 -b 2 -e 3 wgEncodeCrgMapabilityAlign100mer.bed.gz

if [ ! -e wgEncodeCrgMapabilityAlign100mer.bed.gz.tbi ]
then
    echo "could not create tabix index"
    exit 1
else
    echo "mappability index created successfully"
    rm wgEncodeCrgMapabilityAlign100mer.bigWig
    rm wgEncodeCrgMapabilityAlign100mer.wig
fi
