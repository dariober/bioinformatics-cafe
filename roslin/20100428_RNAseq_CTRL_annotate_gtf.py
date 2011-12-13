import sys
import datetime
import re

"""

  Combines a GTF file with the coverage from RNAseq data in pileup format
  Pileup produced by something like '20100428_RNAseq_LPS_pileup_to_genotypes.py'
  using 'samtools pileup -s -c -f'

"""

# ------------------------[ User's input ]-------------------------------------

## Assumes header in pileup
pileup_file= 'C:/Tritume/20100409_RNAseq_CTRL_sscrofa9.56.pileup.genot'

## Assumes NO header in GTF
gtf_file= 'C:/Tritume/Sus_scrofa.Sscrofa9.56.gtf'

## Output file name
pileup_ann= 'C:/Tritume/20100428_Ss9.56_gtf_RNAseq_CTRL.cov'

## Output file comment
## This will be added as first line as comment preceeded by timestamp
comment= 'Created by 20100428_RNAseq_CTRL_annotate_gtf.py. GTF file completed with RNAseq coverage from RNAseq_CTRL lib'

# -----------------------------------------------------------------------------

pileup_f= open(pileup_file, 'r')
pileup_f.readline() ## Discard header line

pileup={}
line= pileup_f.readline()
line= line.split('\t')

chr= line[1]
pileup[chr]= {}
pos= int(line[2])
cov= int(line[8])
pileup[chr][pos]= cov

counter = 1
while True:
    line= pileup_f.readline()
    if line == '':
        break
    line= line.split('\t')
    while line[1] == chr:
        line= pileup_f.readline()
        if line == '':
            break
        line= line.split('\t')
        pos= int(line[2])
        cov= int(line[8])
        pileup[chr][pos]= cov
        counter += 1
        if counter % 250000 == 0:
            print('Processed line: ' + str(counter))            
    if line == '':
        break
    chr= line[1]
    pileup[chr]= {}
    pos= int(line[2])
    cov= int(line[8])
    pileup[chr][pos]= cov
    if line == '':
        break
    counter += 1
    if counter % 250000 == 0:
        print('Processed line: ' + str(counter))            
        
print('Positions in pileup file ' + pileup_file + ': ' + str(counter))

pileup_f.close()

#------------------------[ Annotation ]----------------------------------------

gtf_fopen= open(gtf_file, 'r')
pileup_out= open(pileup_ann, 'w')

comment= '## ' + (datetime.datetime.now()).isoformat(' ') + ' ' + comment + '\n'
pileup_out.write(comment)

header= ['transcript_id', 'rname', 'exon_number', 'f_start', 'f_end', 'coverage']
header= '\t'.join(header)
pileup_out.write(header + '\n')

prev_chr= '5'
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
    
    coverage_line=[]
    for pos in range(f_start, f_end+1):
        try:
            coverage_line.append(str(chr[pos]))
        except:
            coverage_line.append(str(0))
    if set(coverage_line) == set(['0']):
        coverage_line=''
    else:
        coverage_line= ','.join(coverage_line)
    gtf_line= [transcript_id, rname, exon_number, str(f_start), str(f_end)]
    gtf_line.append(coverage_line)
    gtf_line= '\t'.join(gtf_line) + '\n'
    pileup_out.write(gtf_line)

gtf_fopen.close()
pileup_out.close()

sys.exit()