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

wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromOut.tar.gz -O /tmp/rmsk/chromOut.tar.gz
tar xvzf /tmp/rmsk/chromOut.tar.gz -C /tmp/rmsk

cat /tmp/rmsk/*/*.out | grep 'L1HS' | awk '{print $5"\t"$6"\t"$7"\t"$11"\t"$10"\t"$9}' | sed -e 's/^chr//' -e 's/C$/-/' > hg19.te.disctgt.txt
cat /tmp/rmsk/*/*.out | grep 'L1PA2' | awk '{print $5"\t"$6"\t"$7"\t"$11"\t"$10"\t"$9}' | sed -e 's/^chr//' -e 's/C$/-/' >> hg19.te.disctgt.txt
cat /tmp/rmsk/*/*.out | grep 'L1PA3' | awk '{print $5"\t"$6"\t"$7"\t"$11"\t"$10"\t"$9}' | sed -e 's/^chr//' -e 's/C$/-/' >> hg19.te.disctgt.txt
cat /tmp/rmsk/*/*.out | grep 'AluY' | awk '{print $5"\t"$6"\t"$7"\t"$11"\t"$10"\t"$9}' | sed -e 's/^chr//' -e 's/C$/-/' >> hg19.te.disctgt.txt
cat /tmp/rmsk/*/*.out | grep 'AluS' | awk '{print $5"\t"$6"\t"$7"\t"$11"\t"$10"\t"$9}' | sed -e 's/^chr//' -e 's/C$/-/' >> hg19.te.disctgt.txt
cat /tmp/rmsk/*/*.out | grep 'SVA_' | awk '{print $5"\t"$6"\t"$7"\t"$11"\t"$10"\t"$9}' | sed -e 's/^chr//' -e 's/C$/-/' -e 's/Other/SVA/' >> hg19.te.disctgt.txt


