#!/usr/bin/python

import os
import re

# --------------------------[ Input/Output ]-----------------------------------

# Current working dir:
os.chdir('/exports/work/vet_roslin_nextgen/dario/bwa/output/20110106_fantom5_macro_monocyte_vs_ss9')

# Sam files to be merged (they have to be in current dir set above)
sam_list= ['aln21.sam', 'aln22.sam', 'aln23.sam', 'aln24.sam', 'aln25.sam', 'aln26.sam', 'aln27.sam', 'aln28.sam', 'aln29.sam', 'aln30.sam']

# Output file, in current dir
sam_out= open('mphage_mcyte_f5_hg19_21bp.sam', 'w')

# This is not used:
dataset_id= [21, 22, 23, 24, 25, 26, 27, 28, 29, 30]
# -----------------------------------------------------------------------------

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

    
