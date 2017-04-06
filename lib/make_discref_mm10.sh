#!/bin/bash

if [ ! -d /tmp/rmsk.mm10 ]
then
    mkdir /tmp/rmsk.mm10
fi

if [ ! -d /tmp/rmsk.mm10 ]
then
    echo "cannot create directory: /tmp/rmsk.mm10"
    exit 1
fi

wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/chromOut.tar.gz -O /tmp/rmsk.mm10/chromOut.tar.gz
tar xvzf /tmp/rmsk.mm10/chromOut.tar.gz -C /tmp/rmsk.mm10

cat /tmp/rmsk.mm10/*/*.out | grep 'L1Md' | awk '$2 < 10.0 {print $5"\t"$6"\t"$7"\t"$11"\t"$10"\t"$9}' | sed -e 's/C$/-/' > mm10.te.disctgt.txt
cat /tmp/rmsk.mm10/*/*.out | grep 'LTR' | awk '$2 < 10.0 {print $5"\t"$6"\t"$7"\t"$11"\t"$10"\t"$9}' | sed -e 's/C$/-/' >> mm10.te.disctgt.txt
cat /tmp/rmsk.mm10/*/*.out | grep 'B[12]_M' | awk '$2 < 10.0 {print $5"\t"$6"\t"$7"\t"$11"\t"$10"\t"$9}' | sed -e 's/C$/-/' >> mm10.te.disctgt.txt


