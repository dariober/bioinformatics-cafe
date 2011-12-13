#!/usr/bin/python

import os
import re

os.chdir('/exports/work/vet_roslin_nextgen/dario/bwa/output/20101216_cage_050810_iterate')
sam_list= ['aln21.sam', 'aln22.sam', 'aln23.sam', 'aln24.sam', 'aln25.sam', 'aln26.sam', 'aln27.sam']
dataset_id= [21, 22, 23, 24, 25, 26, 27]

sam_out= open('cage_050810_bwa_20101216.sam', 'w')


for s in sam_list:
    n_map= 0
    n_unm= 0
    n_mult= 0
    sam= open(s)
    for line in sam:
        if line.startswith('@'):
            continue
        mm= re.findall('.?\tX0:i:(\d+)\t.?', line)
        if mm == []:
            sam_out.write(line)
            n_unm += 1
            continue
        if int(mm[0]) == 1:
            sam_out.write(line)
	    n_map += 1
            continue
        else:
            n_mult += 1
            if sam_list[-1] == s:
                """ If you are processing the last file, output also the multimappers """
                sam_out.write(line)
    sam.close()
    print(s + '\t' + str(n_map) + '\tsingle mappers;\t' + str(n_unm) + '\tunmapped;\t' + str(n_mult) + '\tmultimappers.')
sam_out.close()

    
