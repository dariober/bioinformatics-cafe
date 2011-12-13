#!/usr/bin/python

"""
Statistics aboout alignment produced by bwa (i.e. one line per read.)
"""
import re
import os

# ---------------------------[ Input/Output ]----------------------------------

## Directory with SAM files to read (everything ending with '.sam' will be precessed)
insam= '/exports/work/vet_roslin_nextgen/dario/bwa/output/20101222_blackface_dge'

## Output summary file
stats_out= 'bwa_samstats.txt'

# -----------------------------------------------------------------------------

sam_files= os.listdir(insam)
sam_files= [s for s in sam_files if s.endswith('.sam')]


out= open(stats_out, 'w')
out.write('\t'.join(['SAM', 'Tot. reads', 'Single-mappers\t', 'Multi-mappers\t', 'Unmapped\t\n']))
## List containing summary stats
summary= []

for file in sam_files:
    print('Processing file: ' + file)
    n= 0           # Number of reads processed
    n_single= 0    # Number of single-mappers
    n_multi= 0     # Number of multimappers
    n_unmap= 0     # Number unmapped
    sam= open(file)
    for line in sam:
        if line.startswith('@'):
            continue
        n += 1
        mm= re.findall('.?\tX0:i:(\d+)\t.?', line)   ## Get the BWA tag X0:i giving the numebr of alignments (or none if unmapped)
        if mm == []:
            n_unmap += 1
        elif int(mm[0]) == 1:
            n_single += 1
        elif int(mm[0]) > 1:
            n_multi += 1
        else:
            sys.exit('Unrecognized X0:i tag at alignment ' + str(n) + ', line:\n' + line)
    file_summary= '\t'.join([file, str(n), str(n_single), str( round( (n_single/float(n))*100 ,2)  )+'%' , str(n_multi), str( round( (n_multi/float(n))*100 ,2)  )+'%', str(n_unmap), str( round( (n_unmap/float(n))*100 ,2)  )+'%'])
    out.write(file_summary + '\n')
    sam.close()
out.close()