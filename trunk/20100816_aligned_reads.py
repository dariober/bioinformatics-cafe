#!/usr/bin/python

"""
   Extract from a SAM file the aligned reads (i.e. exclude reads flagged with 4) 
   
"""

" ---------------------------------[ Input ]--------------------------------------- "

infile= open('/exports/work/vet_roslin_nextgen/dario/bowtie/output/20100816_human_CAGE_vs_sscrofa9/20100816_human_cage_vs_sscrofa9.sam', 'r')

outfile= open('/exports/work/vet_roslin_nextgen/dario/bowtie/output/20100816_human_CAGE_vs_sscrofa9/20100816_human_cage_vs_sscrofa9.sam.aln', 'w')

"--------------------------------------[]------------------------------------------ "

counter= 0
counter_out= 0

for line in infile:
    if line.startswith('@'):
        continue

    counter += 1
    if counter % 1000000 == 0:
        print('Processed line number: ' + str(counter))

    line= line.split('\t')
    if line[1] == '4' or line[1] == '8	':
        continue
    line= '\t'.join(line)
    outfile.write(line)
    counter_out += 1

infile.close()
outfile.close()

print('\nNumber of reads aligned: ' + str(counter_out) + '\n')