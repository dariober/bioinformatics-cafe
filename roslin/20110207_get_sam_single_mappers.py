##!/usr/bin/python

import os
import sys

def get_flag(line):
    """ Return the SAM flag of the reads (line as string) """
    return( int(line.split('\t')[1]) )

def count_proper_pairs(readgroup):
    """ Takes a list of sam lines and returns a list of sam flags saying for each
    line whether the alignment is properly paired (1) or not (0) """
    pp= []
    for line in readgroup:
        line= line.split('\t')
        if int(line[1]) & 2 == 2:
            pp.append(1)
        else:
            pp.append(0)
    return(pp)

samtools= '/exports/work/vet_roslin_nextgen/dario/samtools/samtools-0.1.10/samtools '

"""
## Convert bam to sam
print('Converting BAM to SAM')
os.system(samtools + ' view accepted_hits.bam > accepted_hits.sam')
"""
## Sort by read name (actually by entire line)
sort_py= "python /exports/work/vet_roslin_nextgen/dario/python_scripts/sort_file-2.4.py"
os.system(sort_py + ' -b 500000 accepted_hits.sam accepted_hits.qnamesorted.sam')
"""

## Iterate through lines to get groups of reads:
n_lines= 0   ## Number of lines in accepted_hits.sam
n_qname= 0   ## Number of unique qnames found
flagstat= {} ## Dictionary with the counts of each sam flag

sam= open(sys.argv[1], 'r')

cur_line= sam.readline().rstrip()
cur_read= cur_line.split('\t')[0]
readgroup= [cur_line]   ## A list to store the same group of reads
while cur_line != '':
    while cur_line != '':
        next_line= sam.readline().rstrip()
        next_read= next_line.split('\t')[0]
        if cur_read == next_read:
            readgroup.append(next_line)
            cur_line= next_line
            cur_read= cur_line.split('\t')[0]
        else:
            cur_line= next_line
            cur_read= cur_line.split('\t')[0]
            break
    if len(readgroup) == 1:
        """ If only one valid alignment is present send it to output
        (read mate is unmapped but this read has only one position)"""
        print(readgroup[0])
    else:
        """ If more valid alignments are reported see if there is only one
        that is properly paired. """
        pp= count_proper_pairs(readgroup)
        if pp.count(1) > 2:
            """ If more than 2 valid pairs are found don't output this read  """
            pass
        else:
            """ Determine which reads are properly paired and output them  """
            for i in range(0, len(pp)):
                if pp[i] == 1:
                    print(readgroup[i])
    readgroup= [cur_line] ## Restart new readgroup
sam.close()