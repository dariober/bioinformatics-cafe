#!/usr/bin/python

import re
import sys
import string
import argparse

parser = argparse.ArgumentParser(description="""

DESCRIPTION
    
    Search for matches to a regex in a fasta file and return a bed file with the
    coordinates of the match and the matched sequence itself. 
    
    With defaults, quadparser.py searches for putative quadruplexes on forward
    and reverse strand using the quadruplex rule described at
    http://en.wikipedia.org/wiki/G-quadruplex#Quadruplex_prediction_techniques.
    
    The defualt regex is '([gG]{3}\w{1,7}){3,}[gG]{3}' and its complement
    produce the same output as in http://www.quadruplex.org/?view=quadbaseDownload
    
    Output bed file has columns:
    1. Name of fasta sequence (e.g. chromosome)
    2. Start of the match
    3. End of the match
    4. ID of the match
    5. Length of the match
    6. Strand 
    7. Matched sequence
    
    Note: Fasta sequences are read in memory one at a time. Also the bed file
    of of the matches are kept in memeory.

EXAMPLE:
    ## Test data:
    echo '>mychr' > /tmp/mychr.fa
    echo 'ACTGnACTGnACTGnTGAC' >> /tmp/mychr.fa
    
    quadparser.py -f /tmp/mychr.fa -r 'ACTG'
        mychr	0	4	mychr_0_4_for	4	+	ACTG
        mychr	5	9	mychr_5_9_for	4	+	ACTG
        mychr	10	14	mychr_10_14_for	4	+	ACTG
        mychr	15	19	mychr_15_19_rev	4	-	TGAC

    
    ls /tmp/mychr.fa | quadparser.py -f - -r 'A\w\wGn'
        mychr	0	5	mychr_0_5_for	5	+	ACTGn
        mychr	5	10	mychr_5_10_for	5	+	ACTGn
        mychr	10	15	mychr_10_15_for	5	+	ACTGn

DOWNLOAD
    quadparser.py is hosted at http://code.google.com/p/bioinformatics-misc/

TODO

    - Better handling of forward and reverse matches (i.e. other than complementing
      the forward regex?).
    - Add support for protein search in reverse (currently --regex is translated
      using the nucleotide dictionary so the protein reverse match should be
      provided separately to the --regexrev/-R arg.)
    - Read sequence from stdin ('less myseq.fa | quadparser.py -f - ...').
    """, formatter_class= argparse.RawTextHelpFormatter)


parser.add_argument('--regex', '-r',
                   type= str,
                   help='''Regex to be searched in the fasta input.
Matches to this regex will have + strand. This string passed to python
re.compile(). The default regex is '([gG]{3}\w{1,7}){3,}[gG]{3}' which searches
for G-quadruplexes.
                                   
                   ''',
                   default= '([gG]{3,}\w{1,7}){3,}[gG]{3,}')

parser.add_argument('--regexrev', '-R',
                   type= str,
                   help='''The second regex to be searched in fasta input.
Matches to this regex will have - strand.
By default (None), --regexrev will be --regex complemented by replacing
'actguACTGU' with 'tgacaTGACA'. NB: This means that [a-zA-Z] will be translated
to [t-zT-Z] and proteins are not correctly translated. 
                                   
                   ''',
                   default= None)


parser.add_argument('--fasta', '-f',
                   type= str,
                   help='''Input file in fasta format containing one or more
sequences. Use '-' to get the name of the file from stdin
                                   
                   ''',
                   required= True)

args = parser.parse_args()

" --------------------------[ Check and pare arguments ]---------------------- "

""" Reverse forward match """
intab=  'actguACTGU'
outtab= 'tgacaTGACA'
if args.regexrev is None:
    transtab = string.maketrans(intab, outtab)
    regexrev= args.regex.translate(transtab)
else:
    regexrev= args.regex
    
if args.fasta == '-':
    args.fasta= sys.stdin.readlines()
    if len(args.fasta) > 1:
        sys.exit('\nquadpareser.py: Only one input file at a time can be processed:\n--fasta/-f: %s\n' %(args.fasta))
    args.fasta= args.fasta[0].strip()
    
" ------------------------------[  Functions ]--------------------------------- "

"""                               LIST SORTER
Code to sort list of lists
see http://www.saltycrane.com/blog/2007/12/how-to-sort-table-by-columns-in-python/
"""
import operator

def sort_table(table, cols):
    """ sort a table by multiple columns
        table: a list of lists (or tuple of tuples) where each inner list 
               represents a row
        cols:  a list (or tuple) specifying the column numbers to sort by
               e.g. (1,0) would sort by column 1, then by column 0
    """
    for col in reversed(cols):
        table = sorted(table, key=operator.itemgetter(col))
    return(table)

#if __name__ == '__main__':
#    mytable = (
#        ('Joe', 'Clark', '1989'),
#        ('Charlie', 'Babbitt', '1988'),
#        ('Frank', 'Abagnale', '2002'),
#        ('Bill', 'Clark', '2009'),
#        ('Alan', 'Clark', '1804'),
#        )
#    for row in sort_table(mytable, (1,0)):
#        print(row)

"""                           END of SORTER
-----------------------------------------------------------------------------
"""


psq_re_f= re.compile(args.regex)
psq_re_r= re.compile(regexrev)

ref_seq_fh= open(args.fasta)

ref_seq=[]
line= (ref_seq_fh.readline()).strip()
chr= re.sub('^>', '', line)
line= (ref_seq_fh.readline()).strip()
gquad_list= []
while True:
    while line.startswith('>') is False:
        ref_seq.append(line)
        line= (ref_seq_fh.readline()).strip()
        if line == '':
            break
    ref_seq= ''.join(ref_seq)
    for m in re.finditer(psq_re_f, ref_seq):
        quad_id= str(chr) + '_' + str(m.start()) + '_' + str(m.end()) + '_for'
        gquad_list.append([chr, m.start(), m.end(), quad_id, len(m.group(0)), '+', m.group(0)])
    for m in re.finditer(psq_re_r, ref_seq):
        quad_id= str(chr) + '_' + str(m.start()) + '_' + str(m.end()) + '_rev'
        gquad_list.append([chr, m.start(), m.end(), quad_id, len(m.group(0)), '-', m.group(0)])
    chr= re.sub('^>', '', line)
    ref_seq= []
    line= (ref_seq_fh.readline()).strip()
    if line == '':
        break

gquad_sorted= sort_table(gquad_list, (0,1,2,3))

for line in gquad_sorted:
    line= '\t'.join([str(x) for x in line])
    print(line)
sys.exit()
