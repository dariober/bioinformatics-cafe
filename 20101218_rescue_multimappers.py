#!/usr/bin/python

"""
Extract reads froma fastq file that appear to be multimappers
in a BWA sam file. Such reads have one more base added to the 5' end
(sam file being produced with reads trimmed n bases from the right)
"""

import re
import sys
# -----------------------------[ Input/Output ]--------------------------------

# This is the file without adapters but with complete tags (27 bp)
fastq= open('/exports/work/vet_roslin_nextgen/dario/fastq/20100805_CAGE_AM/CAGE_050810_no_adapter.fq', 'r')

# Latest output of bwa
sam= open('/exports/work/vet_roslin_nextgen/dario/bwa/output/20101218_cage_050810_increment/bwa_output.sam')

# Latest FASTQ input that was used for bwa:
multi_fastq= open('/exports/work/vet_roslin_nextgen/dario/bwa/output/20101218_cage_050810_increment/bwa_input.fq', 'w')

# ------------------------------[ Dictionary of multimappers ]-----------------

mmdict= {}
# n= 0
for line in sam:
    if line.startswith('@'):
        continue
    """ Extract the digits in the X0:i tag """
    mm= re.findall('.?\tX0:i:(\d+)\t.?', line)
#    n += 1
#    print(line)
#    print(mm)
#    if n > 10:
#        sys.exit()
    if mm == []:
        """ If X0 tag is not found, the read is unmapped """
        continue
    mm= int(mm[0])
    if mm > 1:
        line= line.split('\t')
        qname= line[0]
        mmdict[qname]= mm
        seq_len= len(line[9])
print('Number of multimapping reads: ' + str(len(mmdict)))
print('Sequence lenght in sam: ' + str(seq_len))
sam.close()

n= 0
while True:
    qname= fastq.readline().rstrip('/1\n').lstrip('@')
    seq= fastq.readline()
    comment= fastq.readline()
    qual= fastq.readline()
    if qname == '':
        break
    if mmdict.has_key(qname):
        multi_fastq.write('@' + qname + '/1\n')
        multi_fastq.write(seq[0:(seq_len+1)] + '\n')
        multi_fastq.write(comment)
        multi_fastq.write(qual[0:(seq_len+1)] + '\n')
        n += 1
multi_fastq.close()
fastq.close()
print('Reads sent to output: ' + str(n))
