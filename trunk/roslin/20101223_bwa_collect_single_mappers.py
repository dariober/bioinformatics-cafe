#!/usr/bin/python

#
# Collect in a single file all the single mapping reads found in a
# list of sam files produced by BWA
#

# -------------------------[ Input/Output ]------------------------------------

sam_dir= '/exports/work/vet_roslin_nextgen/dario/bwa/output/20101222_blackface_dge'

## Note: This name shouldn't end in '.sam' otherwise it will be processed itself in the loop.
sam_out= 'LN_all.single.sam'

# -----------------------------------------------------------------------------

import re
import os

sam_files= os.listdir(sam_dir)
sam_files= [s for s in sam_files if s.endswith('.sam')]
outfile= open(sam_out, 'w')

for file in sam_files:
    print('Processing file: ' + file)
    sam_id= file.replace('.fq.sam', '') ## Here is where the SAM ID is defined
    n= 0           # Number of reads processed
    n_single= 0    # Number of single-mappers
    sam= open(file)
    for line in sam:
        if line.startswith('@'):
            continue
        n += 1
        mm= re.findall('.?\tX0:i:(\d+)\t.?', line)   ## Get the BWA tag X0:i giving the numebr of alignments (or none if unmapped)
        if mm == []:
            continue
        elif int(mm[0]) == 1:
            outfile.write(sam_id + '\t' + line)
            n_single += 1
        else:
            continue
    print(str(n_single) + '/' + str(n) + ' lines in output.')
    sam.close()
outfile.close()
