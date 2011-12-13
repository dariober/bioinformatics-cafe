#!/usr/bin/python

"""

Removes reads mapped in one SAM file from another SAM file.
This is used to remove tags mapped against rDNA from the SAM
file produced from the whole genome alignment.

"""

# -----------------------[ Input/Output ]--------------------------------------

## SAM file containing the reads to remove (e.g. rDNA)
sam_rdna= open('/exports/work/vet_roslin_nextgen/dario/bwa/output/20101216_CAGE_050810_hg_rdna/20101216_cage_vs_hg_rDNA_21bp.sam')

## SAM file to be filtered 
sam= open('/exports/work/vet_roslin_nextgen/dario/bwa/output/20101216_cage_050810_iterate/cage_050810_bwa_20101216.sam')

## SAM output
sam_filtered= open('/exports/work/vet_roslin_nextgen/dario/bwa/output/20101216_cage_050810_iterate/cage_050810_bwa_20101216.clean.sam', 'w')

# -----------------------------------------------------------------------------

## Dictionary of reads to filter out

filter= {}
n= 0
for line in sam_rdna:
    if line.startswith('@'):
        continue
    n += 1
    if n % 500000 == 0:
        print(str(n) + ' reads processed.')
    line= line.split('\t')
    if line[2] == '*':
        continue
    filter[line[0]]= filter.get(line[0], 0)
    
print(str(len(filter)) + ' reads to be filtered out.\n\n')
sam_rdna.close()

n= 0
n_out= 0
for line in sam:
    if line.startswith('@'):
        continue
    n += 1
    if n % 500000 == 0:
        print(str(n) + ' reads processed.')
    qname= line.split('\t')[0]
    if filter.has_key(qname):
        continue
    sam_filtered.write(line)
    n_out += 1
print(str(n_out) + ' reads retained.')
sam.close()
sam_filtered.close()



