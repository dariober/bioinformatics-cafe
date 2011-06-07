#!/usr/bin/python

#
#  Removes n bases from the 3' (left) of FASTQ reads
#

# ---------------------------------[ I/O ]---------------------------------------

fastq_file= '/exports/work/vet_roslin_nextgen/dario/fastq/20100805_CAGE_AM/CAGE_050810_no_adapter.fq'
fastq_out= 'CAGE_050810_21bp.fq'
seq_len= 21 ## Length of the sequence after trimming

# -------------------------------------------------------------------------------

fastq= open(fastq_file)
out= open(fastq_out, 'w')

n_reads= 0
while True:
    qname= fastq.readline()
    if qname == '':
        break
    seq= fastq.readline().rstrip('\n')
    comment= fastq.readline()
    qual= fastq.readline().rstrip('\n')
    n_reads += 1
    out.write(qname)
    out.write(seq[0:seq_len] + '\n')
    out.write(comment)
    out.write(qual[0:seq_len] + '\n')
fastq.close()
out.close()
print('Reads processed: ' + str(n_reads))

    
   