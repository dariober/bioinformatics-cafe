#!/usr/bin/env python

import sys

docstring= """
DESCRIPTION
===========
Parse the output of samtools mpileup to sum the number of bases in Plasmodium
berghei chromosomes and in mouse chromosomes.

See code for how counts are assigned.

USAGE
=====
samtools mpileup <args> | pileupParserPbergMouse.py -

OUTPUT
======
Tab separated line to stdout:
# Count of bases in P.berghei
# Count of bases in mouse
# Count of bases unassigned to the above
"""

if len(sys.argv) != 2 or sys.argv[1] in ['-h', '--help']:
    print(docstring)
    sys.exit(1)

if sys.argv[1] == '-':
    fin= sys.stdin
else:
    fin= open(sys.argv[1])

pb_cnt= 0
mm_cnt= 0
unk_cnt= 0
for line in fin:
    line= line.strip().split('\t')
    chrom= line[0]
    cnt= int(line[3])
    if chrom.startswith('berg') or chrom.startswith('PBANKA'):
        pb_cnt += cnt
    elif chrom.startswith('chr'):
        mm_cnt += cnt
    else:
        unk_cnt += cnt
print('\t'.join([str(pb_cnt), str(mm_cnt), str(unk_cnt)]))
fin.close()
