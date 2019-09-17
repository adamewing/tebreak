#!/bin/bash

if [ ! -d /tmp/rmsk ]
then
    mkdir /tmp/rmsk
fi

if [ ! -d /tmp/rmsk ]
then
    echo "cannot create directory: /tmp/rmsk"
    exit 1
fi

wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.out.gz -O /tmp/rmsk/hg38.fa.out.gz

zcat /tmp/rmsk/hg38.fa.out.gz | grep 'L1HS' | awk '{print $5"\t"$6"\t"$7"\t"$11"\t"$10"\t"$9}' | sed -e 's/^chr//' -e 's/C$/-/' > hg38.te.disctgt.txt
zcat /tmp/rmsk/hg38.fa.out.gz | grep 'L1PA2' | awk '{print $5"\t"$6"\t"$7"\t"$11"\t"$10"\t"$9}' | sed -e 's/^chr//' -e 's/C$/-/' >> hg38.te.disctgt.txt
zcat /tmp/rmsk/hg38.fa.out.gz | grep 'L1PA3' | awk '{print $5"\t"$6"\t"$7"\t"$11"\t"$10"\t"$9}' | sed -e 's/^chr//' -e 's/C$/-/' >> hg38.te.disctgt.txt
zcat /tmp/rmsk/hg38.fa.out.gz | grep 'AluY' | awk '{print $5"\t"$6"\t"$7"\t"$11"\t"$10"\t"$9}' | sed -e 's/^chr//' -e 's/C$/-/' >> hg38.te.disctgt.txt
zcat /tmp/rmsk/hg38.fa.out.gz | grep 'AluS' | awk '{print $5"\t"$6"\t"$7"\t"$11"\t"$10"\t"$9}' | sed -e 's/^chr//' -e 's/C$/-/' >> hg38.te.disctgt.txt
zcat /tmp/rmsk/hg38.fa.out.gz | grep 'SVA_' | awk '{print $5"\t"$6"\t"$7"\t"$11"\t"$10"\t"$9}' | sed -e 's/^chr//' -e 's/C$/-/' -e 's/Other/SVA/' >> hg38.te.disctgt.txt


