#!/usr/bin/env python

import sys
import argparse
import os
import gzip 
import atexit

parser = argparse.ArgumentParser(description= """
""", formatter_class= argparse.RawTextHelpFormatter, prog= os.path.basename(__file__))

parser.add_argument('input', nargs= '?', default= '-',
    help='''Input file or - for stdin. Default %(default)s''')

parser.add_argument('--comment', '-c', action= 'store_true',
    help='''Print also the comment line. Default is to skip it''')

parser.add_argument('--no_newline', '-n', action= 'store_true',
    help='''Do not separate records with a newline''')

parser.add_argument('--skip', '-s', type= int, default= 0,
    help='''Skip this many records before start print. Default %(default)s''')

parser.add_argument('--offset', type= int, default= 33,
    help='''Offset to interpret ASCII character as phred score, e.g. 33 for Sanger encoding
64 for old Illumina. Default %(default)s''')

parser.add_argument('--version', action='version', version='%(prog)s 0.1.0')    

# -----------------------------------------------------------------------------

NT_COLS= {'A': '12',
         'C': '9', 
         'G': '2',
         'T': '11',
         'NA': '0'}

def get_reader(x):
    if x == '-':
        return sys.stdin
    elif x.endswith('.gz'):
        return gzip.open(x, 'rt')
    else:
        return open(x)

def format_nt(x):
    """Add formatting to nucleotide x
    """
    X= x.upper()
    if(X in NT_COLS):
        return '\033[48;5;231;38;5;%sm%s\033[0m' % (NT_COLS[X], x)
    else:
        return '\033[48;5;231;38;5;%sm%s\033[0m' % (NT_COLS['NA'], x)
        
def format_phred(x, offset= 33):
    """Format the character x interpreted as phred score
    """
    phred= ord(x) - offset
    bold= ''
    if phred <= 2:
        col= 196     # 196:red1
        bold= '1;'
    elif phred <= 8:
        col= 166     # 166:darkorange3
        bold= '1;'
    elif phred <= 12:
        col= 178     # 178:gold3
        bold= '1;'
    elif phred <= 22:
        col= 10      # 13:fuchsia
    elif phred <= 27:
        col= 57      # 57:blueviolet
    elif phred <= 32:
        col= 27      # 27:dodgerblue
    elif phred <= 37:
        col= 39      # 39:deepskyblue1
    elif phred <= 250:
        col= 21      # 21:blue
    else:
        col= 9       # 9:red
    return '\033[%s48;5;231;38;5;%sm%s\033[0m' % (bold, col, x)

def terminal_width():
    try:
        rows, columns = os.popen('stty size', 'r').read().split()
    except:
        columns= 0
    return columns

def line_filler(line, width):
    fill= ' ' * (width - len(line))
    return '\033[48;5;231m' + fill + '\033[0m'

def clear():
    """On exit do this.
    """
    sys.stderr.write('\033[0m') # Clear formatting

if __name__ == '__main__':
    args= parser.parse_args()
    fin= get_reader(args.input)
    sys.stderr.write('\033[48;5;231m')
    twd= int(terminal_width())
    try:
        i= 1
        nrecs= 0
        for line in fin:
            if nrecs >= args.skip:
                line= line.strip()
                if i == 2:
                    fmt= ''.join([format_nt(x) for x in line])
                elif i == 3:
                    if args.comment: 
                        fmt= '\033[48;5;231;38;5;250m' + line + '\033[0m'
                    else:
                        fmt= ''
                elif i == 4:
                    fmt= ''.join([format_phred(x) for x in line])
                else:
                    fmt= '\033[48;5;231;38;5;250m' + line + '\033[0m'
                sep= ''
                if not args.no_newline and i == 4:
                    sep= '\n'
                if i != 3 or (i == 3 and args.comment):
                    print((fmt + line_filler(line, twd) + sep))
            i += 1
            if i > 4:
                i= 1
                nrecs += 1
    except IOError as e:
        sys.exit()

    atexit.register(clear)

