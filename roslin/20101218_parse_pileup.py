#!/usr/bin/python

"""
Extract defined columns from pileup file produced by "samtools pileup -c [...]". 
Two pileups are processed and splice together: forward (+ strand) and reverse (- strand). 
Ouput will look like this

rname	pos	ref_base	read_count	strand
1	10	A		100		+
1	11	C		120		+
1	12	G		140		+
...

"""

import os

# ---------------------------------[ Input/Output ]----------------------------

os.chdir('/exports/work/vet_roslin_nextgen/dario/bwa/output/20101218_cage_050810_increment')

pile_in= ['cage_050810_bwa_20101218.clean.single.forward.pileup', 'cage_050810_bwa_20101218.clean.single.reverse.pileup']

## Strand for the first and second file in pile_in (make sure strands and files match!)
strands= ['+', '-']

## Columns in pileup to keep in output
keep_col= [0, 1, 2, 7]

## Header line
header= 'rname\tpos\tref_base\tread_count\tstrand\n'

## Output file
pile_out= open('/exports/work/vet_roslin_nextgen/dario/bwa/output/20101218_cage_050810_increment/cage_050810_bwa_20101218.clean.single.pileup', 'w')

# -----------------------------------------------------------------------------

## Add header to 
pile_out.write(header)

for pile in pile_in:
    inpile= open(pile)
    strand= strands.pop(0)
    print('Processing: "' + pile + '" with strand "' + strand + '"...')
    for line in inpile:
        line= line.split('\t')
        out= []
        [out.append(line[i]) for i in keep_col]
        out.append(strand)
        out= '\t'.join(out)
        pile_out.write(out + '\n')
    inpile.close()
pile_out.close()