#!/cygdrive/c/Python26/python
"""
Produce dummy FASTQ file given sequence and read name.
This is to convert FANTOM5 transcription start sites back
to CAGE tags to be mapped against S.scrofa
"""

infile= open('F:/data/20110105_FANTOM5/Macrophage%20-%20monocyte%20derived%20donor1.CNhs10861.11232-116C8.hg19.ctss.seq.txt')
outfile= open('F:/data/20110105_FANTOM5/Macrophage%20-%20monocyte%20derived%20donor1.CNhs10861.11232-116C8.hg19.ctss.fq', 'w')

n= 0
for line in infile:
    if n == 0:
        n += 1
        continue
    line= line.strip()
    line= line.split('\t')
    outfile.write('@' + line[0] + '\n')    ## qname
    outfile.write(line[6] + '\n')          ## sequence
    outfile.write('+\n')                   ## comment
    outfile.write('I'*len(line[6]) + '\n') ## quality (sanger)
    if n % 100000 == 0:
        print('%d read processed.' %n)
    n += 1
n -= 1
outfile.close()
infile.close()
print('%d reads sent to output'%n)
    
    

