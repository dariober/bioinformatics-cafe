import sys
import datetime
import inspect
import re

"""
  Combines a GTF file with the coverage from RNAseq data in pileup format
  Pileup produced by something like 'samtools pileup -s -c -f'
  
  Pileup format should look like (tab separated):
  MT      1       C       C       36      0       60      3       ^~.^~.^~.       `a`     ~~~
  MT      2       A       W       1       1       60      12      ...^~,^~,^~,^~,^~,^~.^~T^~,^~.  \abaYZVaXaa`    ~~~~~~~~~~~~

"""

# ------------------------[ User's input ]-------------------------------------

## No header in pileup from samtools. Check this!
pileup_file= '/exports/work/vet_roslin_nextgen/dario/samtools/output/20100409_RNAseq_LPS_tophat/20100409_RNAseq_LPS_sscrofa9.56.pileup'

## Assumes NO header in GTF
gtf_file= '/exports/work/vet_roslin_nextgen/dario/GTF/Sus_scrofa.Sscrofa9.56.gtf'

## Output file name
pileup_ann= '/exports/work/vet_roslin_nextgen/dario/samtools/output/20100409_RNAseq_LPS_tophat/20100409_RNAseq_LPS_sscrofa9.56.cov'

## Single line comment (First line of the output file commented by ## and followed by timestamp and script name)
## Set to None for no comment line
comment= 'GTF file annotated with coverage from RNAseq. Input from pileup file (samtools pileup -s -c -f) from RNAseq_LPS from tophat (aided by GFF file)'

# -------------------------[ Generate pilup dictionary ]-----------------------

"""
This part transforms the pileup format to a dictionary of chromosomes (key) containing 
a dictionary of positions (key) and coverage (value). Something like:
pileup_dict= {'chr1': {1:10, 2:12, 5:20}, 'chr2': {1:5, 2:20, 10:6}}
"""

print('\nPreparing pilup dictionary...\n')

pileup_f= open(pileup_file, 'r')

## Discard header line if present
#pileup_f.readline()

## Initialize dict
pileup={}

line= pileup_f.readline()
line= line.split('\t')

chr= line[0]
pileup[chr]= {}
pos= int(line[1])
cov= int(line[7])
pileup[chr][pos]= cov

counter = 1
while True:
    line= pileup_f.readline()
    if line == '':
        break
    line= line.split('\t')
    while line[0] == chr:
        line= pileup_f.readline()
        if line == '':
            break
        line= line.split('\t')
        pos= int(line[1])
        cov= int(line[7])
        pileup[chr][pos]= cov
        counter += 1
        if counter % 250000 == 0:
            print('Processed line: ' + str(counter))
    if line == '':
        break
    chr= line[0]
    pileup[chr]= {}
    pos= int(line[1])
    cov= int(line[7])
    pileup[chr][pos]= cov
    if line == '':
        break
    counter += 1
    if counter % 250000 == 0:
        print('Processed line: ' + str(counter))            
        
print('\nPositions in pileup file ' + pileup_file + ': ' + str(counter+1) + '\n')

pileup_f.close()

#------------------------[ Annotation ]----------------------------------------

gtf_fopen= open(gtf_file, 'r')
pileup_out= open(pileup_ann, 'w')

## Add comment line
if comment is not None:
    tstamp= (datetime.datetime.now()).isoformat(' ')
    script= inspect.getfile( inspect.currentframe() )
    comment= comment.replace('\n', ' ')
    comment= '## ' + 'Produced by : "' + script + '" on ' + tstamp + ' ' + comment + '\n'
    pileup_out.write(comment)

header= ['transcript_id', 'rname', 'exon_number', 'f_start', 'f_end', 'strand', 'sum_coverage', 'coverage']
header= '\t'.join(header)
pileup_out.write(header + '\n')

prev_chr= '5'
print('Processing chromosome: ' + prev_chr)

for line in gtf_fopen:
    " Parse GTF line "
    line= line.split('\t')
    if line[2] != 'exon':
        continue
    rname= line[0]

    t= re.compile(r"""transcript_id "([a-zA-Z0-9]*?)";""")
    transcript_id= re.findall(t, line[8])[0]

    ex= re.compile(r"""exon_number "([0-9]*?)";""")
    exon_number= re.findall(ex, line[8])[0]

    f_start= line[3]
    f_end= line[4]
    " Retrieve chr from pileup dict"
    chr= pileup[line[0]]    
    if line[0] != prev_chr:
        print('Processing chromosome: ' + rname)
    prev_chr= line[0]
    f_start= int(line[3])
    f_end= int(line[4])
    strand= line[6]
    
    coverage_line=[]
    for pos in range(f_start, f_end+1):
        try:
            coverage_line.append(chr[pos])
        except:
            coverage_line.append(0)
    if set(coverage_line) == set([0]):
        " Ouput a blank if coverage is 0 across the whole exon "
        coverage_line=''
        sum_coverage= '0'
    else:
        sum_coverage= str(sum(coverage_line))
        coverage_line= [str(x) for x in coverage_line]
        if strand == '-':
            coverage_line.reverse()
        coverage_line= ','.join(coverage_line)
    gtf_line= [transcript_id, rname, exon_number, str(f_start), str(f_end), strand, sum_coverage]
    gtf_line.append(coverage_line)
    gtf_line= '\t'.join(gtf_line) + '\n'
    pileup_out.write(gtf_line)

gtf_fopen.close()
pileup_out.close()

sys.exit()