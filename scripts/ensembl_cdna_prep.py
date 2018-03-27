#!/usr/bin/env python

import sys

if len(sys.argv) == 2:
    seq = []
    with open(sys.argv[1], 'r') as _:
        for line in _:
            if line.startswith('>'):
                if seq:
                    print seq[0] + '\n' + ''.join(seq[1:]) + 'A'*50

                seq = [line.strip().split()[0]]
                continue

            seq.append(line.strip())

    if seq:
        print seq[0] + '\n' + ''.join(seq[1:]) + 'A'*50

else:
    sys.exit('usage %s <ensembl.fa>' % sys.argv[0])
