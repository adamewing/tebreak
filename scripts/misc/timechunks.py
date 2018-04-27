#!/usr/bin/env python

'''
Tool for assessing runtimes over genome chunks, useful for troubleshooting

Matches these:
2016-04-29 11:03:44,695 Processing chunk: 1:51177340-51264889 ...


With these:
2016-04-29 11:03:49,902 Finished chunk: 1:51177340-51264889

2016-10-01 07:44:18,591 Chunk chr1:711835-755970: Mapping 0 breakends ...

2016-10-02 02:17:37,840 Chunk chr16:16353178-16397313: Summarising insertions ...

'''

import sys
import re
import time
import datetime

start  = {}
finish = {}

if len(sys.argv) == 2:
    with open(sys.argv[1], 'r') as log:
        for line in log:
            if re.search('Processing chunk:', line):
                c = line.strip().split()

                date  = c[0]
                hms   = c[1]
                chunk = c[-2]

                t = '%s %s' % (date, hms)
                start[chunk] = time.mktime(datetime.datetime.strptime(t, "%Y-%m-%d %H:%M:%S,%f").timetuple())

            if re.search('Summarising', line) or re.search('Mapping 0 breakends', line):
                c = line.strip().split()
                date  = c[0]
                hms   = c[1]
                chunk = c[3].rstrip(':')

                t = '%s %s' % (date, hms)
                finish[chunk] = time.mktime(datetime.datetime.strptime(t, "%Y-%m-%d %H:%M:%S,%f").timetuple())

            if re.search('Finished chunk:', line):
                c = line.strip().split()

                date  = c[0]
                hms   = c[1]
                chunk = c[-1]

                t = '%s %s' % (date, hms)
                finish[chunk] = time.mktime(datetime.datetime.strptime(t, "%Y-%m-%d %H:%M:%S,%f").timetuple())


    for chunk in start:
        if chunk in finish:
            print '%s\t%d' % (chunk, finish[chunk] - start[chunk])
        else:
            print '%s\tUNFINISHED' % chunk

else:
    sys.exit('usage: %s <tebreak log (stderr output)>' % sys.argv[0])
