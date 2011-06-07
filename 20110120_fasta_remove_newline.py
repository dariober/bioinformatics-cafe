"""
   Remove newline from sequences in fasta file
"""

import os
os.chdir('D:/Tritume')
fasta= open('sequence_upstream_hs.fa')
fasta_new= open('sequence_upstream_hs.new.fa', 'w')

n= 0
for line in fasta:
    if line.startswith('>') and n != 0:
        line= '\n' + line
    elif n != 0:
        line= line.strip()
    fasta_new.write(line)
    n += 1
fasta_new.close()
fasta.close()