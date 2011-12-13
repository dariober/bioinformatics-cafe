#!/usr/bin/python

#
#  Extract DGE tags from raw FASTQ files and convert to Sanger format
#  

## Solexa to sanger command
sol2std= 'perl /exports/work/vet_roslin_nextgen/dario/miscellanea/fq_all2std.pl sol2std '

## Inpout FASTQ dir
fastq_dir= '/exports/work/vet_roslin_nextgen/dario/fastq/20101223_Blackface'

## Input FASTQ files:
fastq= ['LN_109.txt', 'LN124.txt', 'LN_146.txt', 'LN_183.txt', 'LN_21.txt', 'LN_47.txt',
        'LN_57.txt', 'LN_92.txt', 'LN114.txt', 'LN_130.txt', 'LN_173.txt', 'LN_20B.txt',
        'LN_38.txt', 'LN_50.txt', 'LN_58.txt']

# ----------------------------------------

import os
os.chdir(fastq_dir)

for file in fastq:
    print('Processing file: ' + file)
    fin= open(file)
    fastq_tmp= file.replace('.txt', '.tmp')
    fout= open(fastq_tmp, 'w')
    n= 0
    while True:
        line= fin.readline() ## qname
        if line == '':
            break
        fout.write(line)
        
        line= fin.readline()  ## Sequence
        line= 'CATG' + line[0:17]
        fout.write(line + '\n')

        line= fin.readline()  ## Comment
        fout.write(line)
        
        line= fin.readline()  ## Quality
        line= 'hhhh' + line[0:17]
        fout.write(line + '\n')
        n += 1
    print(str(n) + ' reads spliced.')
    fout.close()
    fin.close()
    fastq_out= file.replace('.txt', '.fq')
    print('Converting quality Solexa to Sanger...\n')
    os.system(sol2std + ' ' + fastq_tmp + ' > ' + fastq_out)
    os.remove(fastq_tmp)
        



