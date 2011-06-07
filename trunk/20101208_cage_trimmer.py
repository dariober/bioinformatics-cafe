#!/usr/bin/python

"""
Trim CAGE reads.
Trims CAGE reads starting from the 3'end (right). Bases are chopped off 
if they fall below a givn threshold. Reads too short are also rejected
"""

import time

# -------------------------------[ Input/Output ]------------------------------

## Raw fastq files
fastq_file = '/exports/work/vet_roslin_nextgen/dario/fastq/20100805_CAGE_AM/Beraldi_CAGE_50810_101202_EBRI093151_0060_s_8_sequence.txt'

## Minimum quality below which a base is chopped off the read: 
## Reminder: p-error= 10^(-Q/10) where Q is the phred score.
min_trim_q = 20 

## Discard reads entirely when shorter than this:
min_length = 15 

## Output file;
fastq_out= 'cage_050810_right_trimmed.fq'

# -----------------------------------------------------------------------------

t0= time.time()

fastq= open(fastq_file)
out= open(fastq_out, 'w')
read_counter= 0
n_out= 0
n_rej= 0
r_length= {}
while True:
    qname= fastq.readline().rstrip('\n')
    if qname == '':
        break
    seq= fastq.readline().rstrip('\n')
    comment= fastq.readline().rstrip('\n')
    qual= fastq.readline().rstrip('\n')
    
    read_counter += 1
    if read_counter % 500000 == 0:
        print(str(read_counter) + ' reads processed.')
    seq_length= len(seq)
    n_reject= 0                  ## Number of bases to reject, from the right (3' end)
    for b in qual[::-1]:         ## Start reading from 3' end
        q= ord(b) - 64
        if (seq_length - n_reject) < min_length: ## Check the read hasn't become too short
            break
        if q < min_trim_q:       
            n_reject += 1
        else:
            if n_reject == 0:    ## As soon as a base has quality above the threshold, return the read from the beginning to this base.
                tseq= seq
                tqual= qual
                break
            else:
                tseq= seq[: -n_reject]
                tqual= seq[: -n_reject]
                break
    seq_len= len(tseq)
    if seq_len < min_length:
        n_rej += 1
        continue
    r_length[seq_len]= r_length.get(seq_len, 0) + 1
    out.write(qname + '\n')
    out.write(tseq + '\n')
    out.write(comment + '\n')
    out.write(tqual + '\n')
    n_out += 1
t1= time.time()
print('\nNumber of reads processed: ' + str(read_counter))
print('Number of reads retained:  ' + str(n_out))
print('Number of reads rejected:  ' + str(n_rej))
print('Histogram of read length (read_length:count):\n' + str(r_length))
print('\nRun time: ' + str(round(t1-t0, 2)) + ' sec')