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

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign50mer.bigWig

wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigToWig
chmod +x bigWigToWig

./bigWigToWig wgEncodeCrgMapabilityAlign50mer.bigWig wgEncodeCrgMapabilityAlign50mer.wig

if [ ! -e wgEncodeCrgMapabilityAlign50mer.wig ]
then
    echo "could not create wig file"
    exit 1
fi

grep -v '^#' wgEncodeCrgMapabilityAlign50mer.wig | sed -e s'/^chr//' > wgEncodeCrgMapabilityAlign50mer.bed

tabix/bgzip wgEncodeCrgMapabilityAlign50mer.bed
tabix/tabix -s 1 -b 2 -e 3 wgEncodeCrgMapabilityAlign50mer.bed.gz

if [ ! -e wgEncodeCrgMapabilityAlign50mer.bed.gz.tbi ]
then
    echo "could not create tabix index"
    exit 1
else
    echo "mappability index created successfully"
    rm wgEncodeCrgMapabilityAlign50mer.bigWig
    rm wgEncodeCrgMapabilityAlign50mer.wig
fi
