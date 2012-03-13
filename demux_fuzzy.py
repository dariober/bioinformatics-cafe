#!/home/berald01/.local/bin/python

import Levenshtein
import sys
import argparse
import gzip

parser = argparse.ArgumentParser(description="""

DESCRIPTION
    De-multiplex a FASTQ file into separate files on the bases of the barcode sequence
    found on the read name. The output files are gzipped.
    
    demux_fuzzy.py rescues imperfect matches provided that
    1) A unique best match can be found
    2) The edit distance between sequence and barcodes is less than a given threshold
       (1 by default).

    The edit distance is the Levenshtein distance as implemented in the python package
    Levenshtein
    
    The sample sheet gives the link between barcode and output (like demux).
    It has no header and the first two columns (1st barcode, 2nd file name) must be separated
    by space (*not* TAB). Additional columns are ignored.
    
-------sample sheet example ---------
TAGCTTA demultiplex_file-1.fq
ACTTGAA demultiplex_file-2.fq
-------------------------------------
    
USAGE:
    demux_fuzzy.py -f <fastqfile> -s <sample sheet>
    
    ## Read from stdin a zipped file. Use -f -:
    gunzip -c fastq.fq.gz | demux_fuzzy.py -f - -s samplesheet.txt

REQUIREMENTS:
    Python with package Levenshtein (http://pypi.python.org/pypi/python-Levenshtein/)
    and argparse (http://pypi.python.org/pypi/argparse)

TODO
   - Allow to output in uncompressed format.
     ...
    """, formatter_class= argparse.RawTextHelpFormatter)


parser.add_argument('--fastq', '-f',
                   type= str,
                   help='''Input FASTQ file to de-multiplex. Use - to read from stdin
                                   
                   ''',
                   required= True)

parser.add_argument('--samplesheet', '-s',
                   type= str,
                   help='''Sample sheet with barcodes (1st colulmn) and output
files (2nd column). Additional columns are ignored.
Memo: Columns are space (not tab) separated.
                                   
                   ''',
                   required= True)


parser.add_argument('--distance', '-d',
                   type= int,
                   help='''Maximum edit distance between barcodes found in the sample
sheet and sequence read from the fastq file. Default 1. (Memo: ambiguos barcodes
are always discarded)
                                   
                   ''',
                   required= False,
                   default= 1)

args = parser.parse_args()

#------------------------[ Functions ]-----------------------------------------

def read_fastq_line(fastq_fh):
    """
    Given the file handle fastq_fh (e.g. fastq_fh= open('myfastqfile')), reads
    a chunk of 4 lines.
    """
    fqline= [fastq_fh.readline().strip()]
    fqline.append(fastq_fh.readline().strip())
    fqline.append(fastq_fh.readline().strip())
    fqline.append(fastq_fh.readline().strip())
    return(fqline)
    
def read_samplesheet(sample_sheet):
    """
    Read the file sample_sheet containing. Format is space or tab separated with barcode
    sequence (first column) and output file (second column). Additional columns
    ignored.
    Returns a dictionary of codes and files:
    {'ACTACT': 'file1.fq', 'ACGACG': 'file2.fq', ...}
    """
    ssheet= open(sample_sheet).readlines()
    ssheet= [x.strip() for x in ssheet if x.strip() != ''] ## Strip leading and trailing blanks and blank lines
    ssheet= [x.split(' ') for x in ssheet]
    codes= [x[0] for x in ssheet]
    if len(codes) != len(set(codes)):
        sys.exit('Duplicate barcodes found in sample sheet')
    code_dict= {}
    print('\nSample sheet: %s' %(sample_sheet))
    for line in ssheet:
        print('\t'.join(line))
        ## For each barcode/file store the following
        if line[1].endswith('.gz'):
            filename= line[1]
        else:
            filename= line[1] + '.gz'
        code_dict[line[0]]= [
                  filename,                  ## Output file name
                  gzip.open(filename, 'wb'),  ## Output file handle
                  0                             ## Counter for number of reads in this file
                  ]
    print('')
    return(code_dict)

# -----------------------------------------------------------------------------

if args.fastq == '-':
    fh= sys.stdin
else:
    fh= open(args.fastq)
    
code_dict= read_samplesheet(args.samplesheet)
codes= []
for k in code_dict:
    codes.append(k)
    
n= 0
n_lost= 0
while True:
    fqline= read_fastq_line(fh)
    if fqline[0] == '':
        break
    n += 1
    if n % 1000000 == 0:
        print(n)
    bcode= fqline[0][-9:-2]
    if bcode in code_dict:
        ## Test for perfect match
        code_dict[bcode][1].write('\n'.join(fqline) + '\n')
        code_dict[bcode][2] += 1
    elif args.distance == 0 or bcode.count('N') > args.distance:
        ## If there are more Ns than allowed by args.distance, don't go to Levenshtein distance at all 
        n_lost += 1
        continue
    else:
        dists= [Levenshtein.distance(bcode, x) for x in codes]
        best_dist= min(dists)
        if best_dist > args.distance or len([x for x in dists if x == best_dist]) > 1:
            n_lost += 1
        else:
            bcode= codes[dists.index(best_dist)]
            code_dict[bcode][1].write('\n'.join(fqline) + '\n')
            code_dict[bcode][2] += 1
fh.close()
perc_lost= round(100*(n_lost/float(n)),2)
print('\nTotal reads: %s' %(n))
print('Lost:        %s (%s%%)' %(n_lost, perc_lost))
print('\nReads in:')
for k in code_dict:
    perc= round(100*(code_dict[k][2]/float(n)),2)
    print('%s\t%s\t%s%%' %(code_dict[k][0], code_dict[k][2], perc,))
    code_dict[k][1].close()
print('')
sys.exit()