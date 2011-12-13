#!/usr/bin/python

"""
Extract single mappers from SAM file produced by BWA (i.e. look at the X0:i: tag)
"""

# ----------------------------------------[ I/O ]------------------------------

# Sam file from which to extract: 
sam= open('/exports/work/vet_roslin_nextgen/dario/bwa/output/20110106_fantom5_macro_monocyte_vs_ss9/mphage_mcyte_f5_hg19_21bp.nordna.sam')

# Output file
sam_single= open('/exports/work/vet_roslin_nextgen/dario/bwa/output/20110106_fantom5_macro_monocyte_vs_ss9/mphage_mcyte_f5_hg19_21bp.nordna.single.sam', 'w')

# -----------------------------------------------------------------------------

import re

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


        


