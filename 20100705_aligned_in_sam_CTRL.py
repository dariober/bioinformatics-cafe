#!/usr/bin/python

"""
  Extract the lines from a SAM files not flagged as unmapped (i.e. flag does not contain 4 or 8).
  The output includes the same SAM header lines found in input but without eventual contig names
  not present in output.

-----------------[ SAM header looks like this (tab separated)]---------------

@HD     VN:1.0  SO:unsorted
@SQ     SN:IS1#ARTEFACT_        LN:768
@SQ     SN:IS2#ARTEFACT_        LN:1331
@SQ     SN:IS3#ARTEFACT_        LN:1258
@SQ     SN:IS5#ARTEFACT_        LN:1195
...
@PG     ID:Bowtie       VN:0.12.5       CL:"/exports/work/vet_roslin_nextgen/dario/bowtie/current/bowtie /exports/work/vet_roslin_nextgen/dario/bowtie/indexes/pig_repbase_20090604_renamed -1 /exports/work/vet_roslin_nextgen/dario/fastq/20091217_RNAseq_AM_LPS_CTRL/s_7_1_sequence.txt -2 /exports/work/vet_roslin_nextgen/dario/fastq/20091217_RNAseq_AM_LPS_CTRL/s_7_2_sequence.txt -q --solexa1.3-quals -a --seedmms 3 --best --strata --sam /exports/work/vet_roslin_nextgen/dario/bowtie/output/20100629_RNAseq_CTRL_repbase/20100629_RNAseq_CTRL_repbase.sam"

-----------------[ SAM alingment looks like (tab separated)]------------------

EBRI093151_0001:7:1:52:880#NNNNNN       77      *       0       0       *       *       0       0       GCGACCCGTGGGCCCCGGTGCAGANNNNNNNNNNN     ABBBBBBBA@BAAAABAA#################     XM:i:0
EBRI093151_0001:7:1:52:880#NNNNNN       141     *       0       0       *       *       0       0       NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN     ###################################     XM:i:0
EBRI093151_0001:7:1:53:1002#NNNNNN      77      *       0       0       *       *       0       0       AAGGGCTGTGTCGCGGATGACATCNNNNNNNNNNN     BBB@BBBB=@6=<6=####################     XM:i:0

"""

import time
import os

# ----------------------------[ Input ]----------------------------------------

infilename_sam= '/exports/work/vet_roslin_nextgen/dario/bowtie/output/20100629_RNAseq_CTRL_repbase/20100629_RNAseq_CTRL_repbase.sam'

outfilename_sam= '/exports/work/vet_roslin_nextgen/dario/miscellanea/20100705_repbase/20100629_RNAseq_CTRL_repbase.aligned.sam'

# -----------------------------------------------------------------------------

t0= time.time()

infile= open(infilename_sam)
outfile= open(outfilename_sam + '.aln.tmp', 'w')    ## Temp file to output alignment lines
outfile_header= open(outfilename_sam + '.tmp', 'w') ## Temp file to output header lines

## Store sam headers here. 
sam_header= []

## Store contig names found in SAM here.
contigs= []

counter_aln= 0
counter_header= 0
counter_out= 0

print('\nStart reading SAM input...\n')

for line in infile:
    if line.startswith('@'):
        ## Process header lines
        sam_header.append(line)
        counter_header += 1
    else:
        line= line.split('\t')
        if (int(line[1]) & 4 == 0) or (int(line[1]) & 8 == 0):
            ## Output reads properly aligned
            contigs.append(line[2])
            line= '\t'.join(line)
            outfile.write(line)
            counter_out += 1
        counter_aln += 1
        if counter_aln % 1000000 == 0:
            print('Read SAM alignment line number: ' + str(counter_aln))

infile.close()
outfile.close()

sam_header_aln= []

for i in xrange(0, len(sam_header)):
    ## Retain header lines where the contig name is present in output (plus any the other headers)
    if sam_header[i].startswith('@SQ'):
        line= sam_header[i].split('\t')
        contig_name= line[1]
        contig_name= contig_name[3:]
        if contig_name not in contigs:
            continue
    sam_header_aln.append(sam_header[i])

for i in sam_header_aln:
    outfile_header.write(i)
outfile_header.close()

print('\nOutputting SAM file complete with header...')
os.system('cat ' + outfilename_sam + '.tmp' + ' ' + outfilename_sam + '.aln.tmp' ' > ' + outfilename_sam)
os.remove(outfilename_sam + '.tmp')
os.remove(outfilename_sam + '.aln.tmp')

t1= time.time()

print('\nNumber of header lines: ' + str(counter_header))
print('Total number of lines in SAM input (excluding header): ' + str(counter_aln))
print('Total number of lines sent to output (excluding header): ' + str(counter_out) + '\n')

print('Total run time: ' + str(round(t1-t0, 2)) + ' sec')