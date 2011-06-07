#!/usr/bin/python

"""
Combines different outputs of HTSeq.scripts.count in a single
file, in database format. The summary lines at the bottom of
each file are stripped.

EXAMPLE:

------ File_1:

ENSSSCG00000000001	0
ENSSSCG00000000002	4
ENSSSCG00000000003	187
...
ENSSSCG00000000010	100
no_feature	3522369
ambiguous	69638
too_low_aQual	0
not_aligned	0
alignment_not_unique	3924823

------ File_2:

ENSSSCG00000000001	0
ENSSSCG00000000002	4
ENSSSCG00000000003	187
...
ENSSSCG00000000010	100
no_feature	3522369
ambiguous	69638
too_low_aQual	0
not_aligned	0
alignment_not_unique	3924823

------- Output:
file_1	ENSSSCG00000000001	0
file_1	ENSSSCG00000000002	4
file_1	ENSSSCG00000000003	187 
...
file_1	ENSSSCG00000000010	100
file_2	ENSSSCG00000000001	0
file_2	ENSSSCG00000000002	4
file_2	ENSSSCG00000000003	187
...
file_1	ENSSSCG00000000010	100

"""

import sys

# -------------------------------[ Input/Output ]------------------------------

## List of input files:
htseq_in= ['/exports/work/vet_roslin_nextgen/dario/tophat/output/20110202_rnaseq_am_ctrl/feature_count.htseq',
           '/exports/work/vet_roslin_nextgen/dario/tophat/output/20110202_rnaseq_am_lps/feature_count.htseq',
           '/exports/work/vet_roslin_nextgen/dario/tophat/output/20110202_rnaseq_bmdm_ctrl/feature_count.htseq',
           '/exports/work/vet_roslin_nextgen/dario/tophat/output/20110202_rnaseq_bmdm_lps/feature_count.htseq']

## List of file identifier (first column in output). Same order as htseq_in !!
file_id= ["20110202_am_ctrl",
          "20110202_am_lps",
          "20110207_bmdm_ctrl",
          "20110207_bmdm_lps"]

## Output_file:
htseq_out= '/exports/work/vet_roslin_nextgen/dario/Tritume/feature_count.htseq_all'

## Should features with zero counts be removed? True/False
remove_zero= True

# -----------------------------------------------------------------------------

if len(htseq_in) != len(file_id):
    sys.exit('Number of files is not equal to the number of identifiers.')

outfile= open(htseq_out, 'w')

for n in range(0, len(htseq_in)):
    fin= open(htseq_in[n])
    id= file_id[n]
    print('Processing file ' + htseq_in[n])
    for line in fin:
        if line.startswith('no_feature\t'):
            break
        sline= line.strip().split('\t')
        if int(sline[1]) == 0 and remove_zero is True:
            continue
        outline= id + '\t' + line
        outfile.write(outline)
    fin.close()
outfile.close()



