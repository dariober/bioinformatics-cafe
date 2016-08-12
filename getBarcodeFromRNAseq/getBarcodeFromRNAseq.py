#!/usr/bin/env python

import sys
import argparse
import gzip

VERSION= '0.1.0'

parser = argparse.ArgumentParser(description= """
DESCRIPTION
    Extract the rightmost x bases from fastq read 2 and put them in the read name of both read 1 and read 2
    The extracted sequence is placed at the *start* of the read names just after the leading '@'.

USAGE:
    getBarcodeFromRNAseq.py -inR1 r1.fq -inR2 r2.fq -outR1 out1.fq -outR2 out2.fq

Gzip input and output:
    getBarcodeFromRNAseq.py -inR1 <(zcat r1.fq.gz) \\
                            -inR2 <(zcat r2.fq.gz) \\
                            -outR1 >(gzip > out1.fq.gz) \\
                            -outR2 >(gzip > out2.fq.gz)

Interleaved output 
    getBarcodeFromRNAseq.py -inR1 r1.fq -inR2 r2.fq

""", formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--inR1', '-inR1',
                   required= True,
                   help='''First-in-pair input fastq file. Use - for stdin
                   ''')

parser.add_argument('--inR2', '-inR2',
                   required= True,
                   help='''Second-in-pair input fastq file. Use - for stdin
                   ''')

parser.add_argument('--outR1', '-outR1',
                   default= '-',
                   help='''First-in-pair output fastq file. Default - (stdout).
                   ''')

parser.add_argument('--outR2', '-outR2', 
                   default= '-',
                   help='''Second-in-pair output fastq file. Default - (stdout).
                   ''')

parser.add_argument('--len', '-l',
                   type= int,
                   default= 12,
                   help='''Length of the barcode to extract from the right-end of read 2. Default 12.
                   ''')

parser.add_argument('--sep', '-s',
                   default= ':',
                   help='''String to separate the barcode from the read name. Default is ':'. 
                   ''')

parser.add_argument('--minlen', '-ml',
                   type= int,
                   default= -1,
                   help='''Skip reads that are shorter than MINLEN *before* the barcode is extracted.
Use -1 (default) to skip reads shorter than the barcode length (--len).
                   ''')

parser.add_argument('--version', action='version', version='%s' %(VERSION))
   
args= parser.parse_args()

def checkPair(r1, r2):
    """Some checks that list r1 and r2 have 4 elements, same name, etc
    """
    r1= [x for x in r1 if x != '']
    r2= [x for x in r2 if x != '']

    if len(r1) != 4 or len(r2) != 4:
        sys.stderr.write("Invalid read:\n")
        sys.stderr.write(str(r1))
        sys.stderr.write(str(r2))
        sys.exit(1)

    name1= r1[0].strip().split()[0]
    name2= r2[0].strip().split()[0]
    if name1 != name2:
        sys.stderr.write("Read names are not the same for\n")
        sys.stderr.write(str(r1))
        sys.stderr.write(str(r2))
        sys.exit(1)

    if not name1.startswith('@') or not name2.startswith('@'):
        sys.stderr.write("Invalid read name(s)\n")
        sys.stderr.write(str(r1))
        sys.stderr.write(str(r2))
        sys.exit(1)

    if (len(r1[1]) != len(r1[3])) or (len(r2[1]) != len(r2[3])):
        sys.stderr.write("Inconsistent sequence and quality lengths:\n")
        sys.stderr.write(str(r1))
        sys.stderr.write(str(r2))
        sys.exit(1)

def extractBarcode(r2, r1, xlen, sep):
    """Parse list r2 to extract barcode of length xlen and put it in the read name of r1 and r2.
    Lists r1 and r2 are modified in place
    r1= ['@name1', 'ABCDEF', '+', '123456']
    r2= ['@name2', 'abcdef', '+', '123456']
    extractBarcode(r2, r1, 2, ':')

    r1 == ['@ef:name1', 'ABCDEF', '+', '123456']
    r2 == ['@ef:name2', 'abcd', '+', '1234']
    """
    xfrom= len(r2[1]) - xlen
    if xfrom < 0:
        xfrom= 0
    barcode= r2[1][xfrom:] 
    r2[1]= r2[1][0:xfrom] # Trim read
    r2[3]= r2[3][0:xfrom] # Trim quality

    if not r2[0].startswith('@'):
        sys.stderr.write("Unexpected read name %s \n" %(r2))
        sys.exit(1)
    if not r1[0].startswith('@'):
        sys.stderr.write("Unexpected read name %s \n" %(r1))
        sys.exit(1)

    r2[0]=  '@' + barcode + sep + r2[0][1:] ## Add barcode to read name 1 and 2
    r1[0]=  '@' + barcode + sep + r1[0][1:]

if args.inR1 == '-':
    fin1= sys.stdin
else:
    fin1= open(args.inR1)

if args.inR2 == '-':
    fin2= sys.stdin
else:
    fin2= open(args.inR2)

if args.outR1 == '-':
    fout1= sys.stdout
else:
    fout1= open(args.outR1, 'w')

if args.outR2 == args.outR1:
    fout2= fout1
elif args.outR2 == '-':
    fout2= sys.stdout
else:
    fout2= open(args.outR2, 'w')

if args.len <= 0:
    sys.stderr.write("Lenght of barcode must be >0. Got: %s \n" %(args.len))
    sys.exit(1)

minlen= args.minlen
if args.minlen < 0:
    minlen= args.len

while True:
    r1= [fin1.readline(), fin1.readline(), fin1.readline(), fin1.readline()]
    r2= [fin2.readline(), fin2.readline(), fin2.readline(), fin2.readline()]
    if r1[0] == '' and r2[0] == '':
        break
    checkPair(r1, r2)

    r1= [x.strip() for x in r1]
    r2= [x.strip() for x in r2]

    if len(r2[1].strip()) < minlen:
        continue

    extractBarcode(r2, r1, args.len, args.sep)

    fout1.write('\n'.join(r1) + '\n')
    fout2.write('\n'.join(r2) + '\n')

fout1.close()
fout2.close()
fin1.close()
fin2.close()
sys.exit()