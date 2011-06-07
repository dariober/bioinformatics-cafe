#!/usr/bin/python

import time

# ---------------------------------[ User's input ]----------------------------
pair_end_file_1 = '/exports/work/vet_roslin_nextgen/dario/fastq/20091217_RNAseq_AM_LPS_CTRL/s_7_1_sequence.txt' ## Raw fastq files
pair_end_file_2 = '/exports/work/vet_roslin_nextgen/dario/fastq/20091217_RNAseq_AM_LPS_CTRL/s_7_2_sequence.txt'

## Output: Same location as current script. Same filename as input with '.avgq' appended

## Minumum avg quality to retain a read
min_avg_quality= 20

# -------------------------------[ Define functions ]--------------------------

def average(values):
    """Computes the arithmetic mean of a list of numbers.

    >>> print average([20, 30, 70])
    40.0
    """
    return sum(values, 0.0) / len(values)

def fastq_tab(fastq_handle):
    """
    Reads four lines from an opened fastq file 
    and returns a list of four elements each 
    corresponding to a line ([read name, sequence, comment, quality])
    
    Example:
    fastq= open(fastq_file)
    xv= fastq_tab(fastq)
    --> ['@name', 'ACTG', '+comment', 'abCD']

    """
    fq_counter = 1
    tabq= []
    while fq_counter <= 4:
        line = (fastq_handle.readline()).rstrip('\n')
        tabq.append(line) 
        fq_counter += 1
    return(tabq)

def avq_fastq_q(tab_fastq):
    """
    Returns the average quality of a fastq read. Input is a list of four elements
    corresponding to the four lines of a FASTQ block (E.g. ouput of fastq_tab)
    """
    offset= 64 ## To convert ASCII to PHRED
    quality= [ord(x) - offset for x in tab_fastq[3]]

    return(average(quality))
     

# -------------------------------[ Output files ]------------------------------

## Output: Same location as current script. Same filename as input with '.avgq' appended
## 
basen= pair_end_file_1[::-1].find('/')    ## Length of file name (excluding path)
outPe_name_1 = pair_end_file_1[-basen :] + '.avgq'

## File to output single reads (for one a read is rejected in one library but not in the other)
outSe_name_1 = pair_end_file_1[-basen :] + '.single'

basen= pair_end_file_2[::-1].find('/')    ## Length of file name (excluding path)

outPe_name_2 = pair_end_file_2[-basen :] + '.avgq'
outSe_name_2 = pair_end_file_2[-basen :] + '.single'


# -----------------------------[ Start iterating ]-----------------------------
print('\nStart filtering reads...\n')
t0= time.time()

fastq_1= open(pair_end_file_1)
fastq_2= open(pair_end_file_2)

outPe_1= open(outPe_name_1, 'w')
outSe_1= open(outSe_name_1, 'w')

outPe_2= open(outPe_name_2, 'w')
outSe_2= open(outSe_name_2, 'w')

counter= 0
counterPE= 0
counterSE_1= 0
counterSE_2= 0
while True:
    readq_1 = fastq_tab(fastq_1)

    if readq_1[0] == '':
        break

    readq_1.append('')
    avg_q_1 = avq_fastq_q(readq_1)
    
    readq_2 = fastq_tab(fastq_2)
    readq_2.append('')
    avg_q_2 = avq_fastq_q(readq_2)
    
    if avg_q_1 >= min_avg_quality and avg_q_2 >= min_avg_quality:
        outPe_1.write('\n'.join(readq_1))
        outPe_2.write('\n'.join(readq_2))
        counterPE += 1

    elif avg_q_1 >= min_avg_quality:
        outSe_1.write('\n'.join(readq_1))
        counterSE_1 += 1

    elif avg_q_2 >= min_avg_quality:
        outSe_2.write('\n'.join(readq_2))
        counterSE_2 += 1

    counter += 1
    if counter % 1000000 == 0:
        print(str(counter) + ' reads processed')
##        break


fastq_1.close(); fastq_2.close(); 
outPe_1.close(); outPe_2.close()
outSe_1.close(); outSe_2.close()

t1= time.time()

print('\nNumber of reads in raw FASTQ files: ' + str(counter))
print('... of which, after filtering:')
print('    ' + str(counterPE) + ' pair-ended')
print('    ' + str(counterSE_1) + ' single-end (from ' + pair_end_file_1 + ')')
print('    ' + str(counterSE_2) + ' single-end (from ' + pair_end_file_2 + ')')
print('    ' + str(counter - (counterPE + counterSE_1 + counterSE_2)) + ' discarded')
print('\nRun-time: ' + str(round(t1-t0,2)) + ' seconds\n')



            