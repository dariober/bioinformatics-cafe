"""
This script extract the query name and full reference name from the alignment
file produced by BLAST (see the end of the script for example of input)

Output file has the same dir as the input and its name has '.parsed' added to
the name of the input

"""

##------------------------------[ User Input ]---------------------------------

blastf= 'C:/Tritume/blast_insert_100.txt'

##-----------------------------------------------------------------------------

infile= open(blastf, 'r')
outfile= open(blastf + '.parsed', 'w')
outfile.write('qname\trname\n')

rDNA_count= 0
q_count= 0

for line in infile:
#    if line.startswith('Query=') or line.startswith('>'):
#        outfile.write(line)
    if 'ribosomal' in line or 'Ribosomal' in line or 'rDNA' in line or 'rRNA' in line:
        ## print(line)
        rDNA_count += 1
    if line.startswith('Query='):
        q_count +=1
        
print('Matches to some rDNA-like sequence: ' + str(rDNA_count))
print('Number of queries: ' +str(q_count))

infile= open(blastf, 'r')

line= infile.readline()
while 1<2:
    if line.startswith('Query='):
        query_name= ''
        line= line.rstrip('\n')
        line= line.lstrip('Query=  ')
        query_name= line + '\t'
    description_line = ''
    if line.startswith('>'):
        line= line.rstrip('\n')
        description_line = description_line + line
        
        line= infile.readline()
        while 'Length=' not in line:
            line= line.rstrip('\n')
            description_line = description_line + line
            line= infile.readline()
        outfile.write(str(query_name) + str(description_line) + '\n')
    line= infile.readline()
    if line == '':
        break
    
outfile.close()
infile.close()

## ------------------Example of input and output ------------------------------
"""
BLASTN 2.2.22+
Reference: Stephen F. Altschul, Thomas L. Madden, Alejandro
A. Schaffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and
David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new
generation of protein database search programs", Nucleic
Acids Res. 25:3389-3402.


RID: SF9W9H7V01S


Database: All GenBank+EMBL+DDBJ+PDB sequences (but no EST, STS,
GSS,environmental samples or phase 0, 1 or 2 HTGS sequences)
           10,955,361 sequences; 30,277,851,791 total letters
Query=  01
Length=27


                                                                   Score     E
Sequences producing significant alignments:                       (Bits)  Value

gb|EU852965.1|  HIV-1 clone 120FIc11 from Uganda envelope glyc...  38.2    0.78 
gb|AC140298.3|  Mus musculus BAC clone RP24-319I6 from chromos...  38.2    0.78 
gb|EF506888.1|  Homo sapiens fibroblast growth factor 2 (basic...  36.2    3.1  
gb|AC134248.3|  Mus musculus BAC clone RP23-338G21 from chromo...  36.2    3.1  
gb|AC021205.6|  Homo sapiens BAC clone RP11-170N16 from 4, com...  36.2    3.1  
gb|AC157217.2|  Pan troglodytes BAC clone CH251-305B5 from chr...  36.2    3.1  
emb|CT033793.10|  Mouse DNA sequence from clone RP23-21O10 on ...  36.2    3.1  
gb|AC156502.2|  Mus musculus BAC clone RP23-355N1 from 9, comp...  36.2    3.1  
gb|AC133650.4|  Mus musculus BAC clone RP23-453B15 from 9, com...  36.2    3.1  
ref|XM_002534027.1|  Ricinus communis RNA-binding protein, put...  34.2       12

ALIGNMENTS
>gb|EU852965.1| HIV-1 clone 120FIc11 from Uganda envelope glycoprotein (env) 
gene, complete cds
Length=2592

 Score = 38.2 bits (19),  Expect = 0.78
 Identities = 19/19 (100%), Gaps = 0/19 (0%)
 Strand=Plus/Minus

Query  4     GGTTGGGGGTATGGGTCTG  22
             |||||||||||||||||||
Sbjct  2170  GGTTGGGGGTATGGGTCTG  2152


>gb|AC140298.3| Mus musculus BAC clone RP24-319I6 from chromosome 18, complete 
sequence
Length=148703

 Score = 38.2 bits (19),  Expect = 0.78
 Identities = 19/19 (100%), Gaps = 0/19 (0%)
 Strand=Plus/Minus

Query  4      GGTTGGGGGTATGGGTCTG  22
              |||||||||||||||||||
Sbjct  78838  GGTTGGGGGTATGGGTCTG  78820

etc. etc...
"""

## -------------------------[ Ouput ]------------------------------------------