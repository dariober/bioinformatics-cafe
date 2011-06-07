#!/usr/bin/python

## ---------------------------[ pairend_filter.py ]----------------------------
#  
# Reads the two fastq files of a pair end sequence library (Solexa) to filter
# out low quality reads and trim low quality 3'ends.
#
# Filter applied: Tolerate n bases with quality below a given threshold. At the n+1  
# low-quality base (reading the reads 5' -> 3'), chop away the
# read from n+1 onwardsthe 2nd low quality base onward.

# The longer mate of a pair is trimmed to the length of the shorter (as required 
#  by tophat).
#

import sys

## --------------------------------[ User's input ]----------------------------

pair_end_file_1 = '/exports/work/vet_roslin_nextgen/dario/fastq/20091217_RNAseq_AM_LPS_CTRL/s_7_1_sequence.txt' ## Raw fastq files
pair_end_file_2 = '/exports/work/vet_roslin_nextgen/dario/fastq/20091217_RNAseq_AM_LPS_CTRL/s_7_2_sequence.txt'

window_size= 5  ## Size of the sliding window
min_trim_q = 20 ## Minimum average quality of the window (PHRED score)
min_tail = 15   ## Trim right-hand side bases below this PHRED score

min_read_length = 20 ## Discard (filtered) reads when shorter than min_read_length bases

## Output files are in the same dir as the input.
## Output name: input filename + '.filtered' (e.g. s_7_1_sequence.txt.filtered)


## --------------------------[ Define trim_sequence filterer]------------------

def average(values):
    """Computes the arithmetic mean of a list of numbers.

    >>> print average([20, 30, 70])
    40.0
    """
    return sum(values, 0.0) / len(values)

    
def trim_window(seq_quality, window_size, min_q):
    """
    Takes a list of phred scores and
    applies a sliding window (of size 'window_size') starting from the right-hand side and trims one base at a 
    time until the average quality of the window is above the phred score in 'min_q'.
    
    Returns the length of the trimmed sequence. 
    
    """
    
    while len(seq_quality) > 0:
        window_q= seq_quality[(len(seq_quality) - window_size):]
        if average(window_q) < min_q:
            seq_quality.pop()
##            print(seq_quality)
        else:
            break
    return(len(seq_quality))

    
## ----------------------------------------------------------------------------

## Open pair-end files
pe1= open(pair_end_file_1, 'r')
pe2= open(pair_end_file_2, 'r')

## Open output files
pe1_out= open(pair_end_file_1 + '.filtered', 'w')
pe2_out= open(pair_end_file_2 + '.filtered', 'w')

counter = 1 ## This counter used to convert each read (4 lines) to a list

read_pe1 = list() ## List to contain a 4-lines bloc (one read)
read_pe2 = list()

i = 0
j = 0 ## Count number of output reads
k = 0 ## Average length of the reads 
m = 0 ## Count input reads
line = 'file_iterator'

while line != '':
       
    ## --------------------[ Convert reads (4 line bloc) to list ]-------------
    while counter <= 4:
        line = (pe1.readline()).rstrip('\n')
        read_pe1.append(line)
        
        line = (pe2.readline()).rstrip('\n')
        read_pe2.append(line)

        i += 1         
        counter += 1
    counter = 1
    m += 1
    ## ---------------------------[ End of converter ]-------------------------

    if m % (1000000) == 0:
        print('Processing ' + str(m) + ' reads')

    seq_quality_1 = [ord(x) - 64 for x in read_pe1[3]]  ## check -64 is correct!
    seq_quality_2 = [ord(x) - 64 for x in read_pe2[3]]  ## check -64 is correct!

    
    ## Find where to trim the sequences using
    pe1_trim = trim_window(seq_quality_1, window_size, min_trim_q)
    pe2_trim = trim_window(seq_quality_2, window_size, min_trim_q)
    
    ## Remove low quality tail
    pe1_trim = trim_window(seq_quality_1, 1, min_tail)
    pe2_trim = trim_window(seq_quality_2, 1, min_tail)
    
#    if m > 3000:
#        break

    ## Check the length of the trimmed reads.
    ## If one of the two in the pair is too short, drop both (i.e. continue to next read)
    if pe1_trim < min_read_length or pe2_trim < min_read_length:
        read_pe1 = list()
        read_pe2 = list()
        continue

    ## Reset the length of the two reads to the shortest one (tophat doesn't like pairs with different length)
    min_length= min([pe1_trim, pe2_trim])
    
    ## Return trimmed reads
                        # Name       # Sequence                   # Comment    # Quality
    read_pe1_trimmed = [read_pe1[0], read_pe1[1][0 : min_length], read_pe1[2], read_pe1[3][0 : min_length]]
    read_pe2_trimmed = [read_pe2[0], read_pe2[1][0 : min_length], read_pe2[2], read_pe2[3][0 : min_length]]

    ## Check consistency of names bewteen pair reads (just in case)
    if read_pe1[0][:-2] != read_pe2[0][:-2]:
        print('Inconsistency of names for reads at line ' + str(i))
        print('Read from file 1:\n' + read_pe1_trimmed)
        print('Read from file 2:\n' + read_pe2_trimmed)
        break

    ## Write out reads    
    pe1_out.write('\n'.join(read_pe1_trimmed) + '\n')
    pe2_out.write('\n'.join(read_pe2_trimmed) + '\n')
    j += 1
    k = 1.0 * ((k * (j-1)) + min_length) / j
    
    ## Reset     
    read_pe1 = list()
    read_pe2 = list()
    
  
pe1.close()
pe2.close()
pe1_out.close()
pe2_out.close()

print('\n----------------------[ Filtering reads ]------------------------')
print('\nFASTQ files:\n' + pair_end_file_1 + '\n' + pair_end_file_2)
print('\nNumber of input reads: ' + str(m))
print('\nNumber of output reads: ' + str(j))
print('\nAverage read length: ' + str(k))
print('')



## -----------------------------[ Tritume ]------------------------------------
def line_counter(file):
    f = open(file) 
    i = 0
    for line in f:
        i += 1
    f.close()
    return(i)