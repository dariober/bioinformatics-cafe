#!/usr/bin/python

#
#  Removes n bases from the 3' (left) of FASTQ reads
#

# ---------------------------------[ I/O ]---------------------------------------

fastq_file= '/exports/work/vet_roslin_nextgen/dario/bwa/output/20110106_fantom5_macro_monocyte_vs_ss9/Macrophage%20-%20monocyte%20derived%20donor1.CNhs10861.11232-116C8.hg19.ctss.fq'
fastq_out= 'mphage_mcyte_f5_hg19_21bp.fq'
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

    
   