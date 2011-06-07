#!/usr/bin/python

"""
Extract single mappers from SAM file produced by BWA (i.e. look at the X0:i: tag)
"""

import re

sam= open('/exports/work/vet_roslin_nextgen/dario/bwa/output/20101216_cage_050810_iterate/cage_050810_bwa_20101216.clean.sam')
sam_single= open('/exports/work/vet_roslin_nextgen/dario/bwa/output/20101216_cage_050810_iterate/cage_050810_bwa_20101216.clean.single.sam', 'w')

n= 0
for line in sam:
    if line.startswith('@'):
        """ Output header """
        sam_single.write(line)
    mm= re.findall('.?\tX0:i:(\d+)\t.?', line)
    if mm == [] or int(mm[0]) > 1:
        continue
    sam_single.write(line)
    n += 1

sam.close()
sam_single.close()
print(str(n) + ' single-mapping reads found.')


        


