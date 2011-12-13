#!/usr/bin/python

"""
  Extract heterozygous SNPs from a pileup file. Pileup in the format:

MT      1       C       C       33      0       60      2       ^~.^~.  =B      ~~
MT      2       A       A       54      0       60      9       ..^~,^~,^~,^~.^~,^~,^~, BB9?>;67=       ~~~~~~~~~
MT      3       A       A       127     0       60      40      ..,,,.,,,^~.^~,^~.^~,^~.^~,^~,^~.^~.^~,^~,^~,^~,^~.^~,^~,^~,^~.^~.^~,^~c^~,^~.^~.^~.^~.^~.^~,^~,^~,^~,  @C,1/'215B#>:B98@B48/0B2+59A/0<BB#BB#75;        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
MT      3223    A       W       187     221     60      192     ,$,$,$,t,tTTTTTGTTTTTTGTGTT.t.tt.TgTTt.,,t,t,T,tttt,t.tt,gttgtt,,tg,ttg,,tttt.,t,ttT,tTg,,t,tgT,,,gTgttttg,T,t,,,,,tgg.tt,,ttTGt,tt..t,tTTt.tTgtttTt,tttTtgtTttt,TTG,TTT,C,,ttGt,gt.t.T..T..T,^~t^~,^~T^~t^~t	BBBCAB@@B?BBBCBB@BA@BCB<?>>B=2<CBBC@CCCBBCAC@ACBBAA8B?C6BC5??C?A7BBBCCBCB=ABB8ACCBBB>ACBCBBBB8@B:>9??BACCBBBB?BCBABAAAB@@@CC@A:?;C<4:BB9B8@7B;9B;@<@)C60?8@<C@BBA?B<BBB@A>?B<B37?=ABABCBBB@-@B:*	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Output is the first 8 columns of the original pileup plus a column for tyhe count of allele 1 and a column
  for the count of allele 2. (Allele 1: the first in alphabetical order)

"""

# --------------------------------[ Input ]------------------------------------

infile_name= '/exports/work/vet_roslin_nextgen/dario/samtools/output/20100630_RNAseq_CTRL/20100630_RNAseq_CTRL.pileup'

outfile_name= '/exports/work/vet_roslin_nextgen/dario/samtools/output/20100630_RNAseq_CTRL/20100728_CTRL_snp.pileup'

phred_snp= 20 ## PHRED score for SNP calling. SNPs less than this will be discarded

# -----------------------------------------------------------------------------

import time

t0= time.time()

# Heterozygous bases as IUPAC:
# R (A or G), Y (C or T), S (G or C), W (A or T), K (G or T), M (A or C)
snp_codes= {'R': ['A', 'G'], 'Y' : ['C', 'T'], 'S' : ['C', 'G'], 'W' : ['A', 'T'], 'K' : ['G', 'T'], 'M' : ['A', 'C']}

infile= open(infile_name, 'r')
outfile= open(outfile_name, 'w')

counter= 0
counter_snp= 0
for line in infile:
    line= line.split('\t')
    counter += 1
    if counter % 1000000 == 0:
        print('Reading line number: ' + str(counter))
    if (snp_codes.has_key(line[3]) == False) or int(line[5]) < phred_snp:
        continue
    outline= line[:8]
    
    alleles= snp_codes[line[3]]
    allele_1= alleles[0]
    allele_2= alleles[1]

    if line[2] == allele_1:
        """ If one of the alleles is in the reference, count the number of ',' and '.' """
        count_a1= line[8].count('.') + line[8].count(',')
    else:
        count_a1= line[8].count(allele_1) + line[8].count(allele_1.lower())

    if line[2] == allele_2:
        """ If one of the alleles is in the reference, count the number of ',' and '.' """
        count_a2= line[8].count('.') + line[8].count(',')
    else:
        count_a2= line[8].count(allele_2) + line[8].count(allele_2.lower())

    outline.append(allele_1)
    outline.append(str(count_a1))
    outline.append(allele_2)
    outline.append(str(count_a2))

    outline='\t'.join(outline)
    outfile.write(outline + '\n')
    counter_snp += 1
    if counter_snp % 10000 == 0:
        print('Output snp number: ' + str(counter_snp))
outfile.close()
infile.close()
t1= time.time()

print('\nLines in input pileup: ' + str(counter))
print('SNPs sent to output: ' + str(counter_snp))
print('Run time: ' + str(round(t1-t0, 2)) + ' sec\n')
