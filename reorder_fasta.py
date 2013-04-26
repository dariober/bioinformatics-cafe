#!/usr/bin/env python

import sys

docstring= """DESCRIPTION
    Reorder the sequences in a FASTA file according to the order given in a reference
    file. The reference file has one sequence name per line.
USAGE
    reorder_fasta.py <file.fasta> <file.reference>

----------- EXAMPLE -------------
## fasta file
echo '>second_seq
AAAAAAAAAAAA
AAAAAAAAAAAA
AAAA
>first_seq
TTTTTTTTTTTTTT
TTTTTTTTTTTTTT
TTTTTTTTT
>third_seq
CCCCCCCCCCCCCCCCCC' > seq.fasta

## reference file
echo 'first_seq
second_seq
third_seq' > ref.txt

## Reorder fasta according to reference:
reorder_fasta.py seq.fasta ref.txt
>first_seq
TTTTTTTTTTTTTT
TTTTTTTTTTTTTT
TTTTTTTTT
>second_seq
AAAAAAAAAAAA
AAAAAAAAAAAA
AAAA
>third_seq
CCCCCCCCCCCCCCCCCC
"""

if len(sys.argv) != 3:
    sys.exit(docstring)

fasta= open(sys.argv[1])
ref= open(sys.argv[2])

seq_dict= {}
while True:
    line= fasta.readline()
    if line == '':
        break
    if line.strip().startswith('>'):
        seq_name= line.strip()[1:]
        seq_dict[seq_name]= []
    else:
        seq_dict[seq_name].append(line.strip())
fasta.close()
for seq_name in ref:
    seq_name= seq_name.strip()
    print('>' + seq_name)
    print('\n'.join(seq_dict[seq_name]))
ref.close()
sys.exit()
