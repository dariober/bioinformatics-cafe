#!/usr/local/bin/python

import sys
import os
import socket
import hashlib
import datetime
import argparse

parser = argparse.ArgumentParser(description= """

DESCRIPTION:

    Get file stats for the file passed as first positional argument.
    Returns a string to stdout formatted like a dict.
    By default, the md5sum is not computed, unless --md5sum is specified, and the
    md5sum value in output is an empty string.
    The order of the keys os not always the same.

OUTPUT:
    Printed on a single line:
    
        "{'filename': os.path.split(args.infile)[1],
          'path':     os.path.split(args.infile)[0],
          'md5sum':   md5sum(args.infile),
          'hostname': socket.gethostname(),
          'ctime':    datetime.datetime.fromtimestamp(os.path.getctime()),
          'mtime':    os.path.getmtime(),
          'fsize':     os.path.getsize()}"

    IMPORTANT: If for any reasons any of the file stats fail (e.g. the file doesn't
    exist), get_file_stats.py exits and prints nothing, no error is returned

EXAMPLES:
        

TODO:
    
""", formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--infile', '-i',
                   required= True,
                   help='''File to stat
                   ''')

parser.add_argument('--md5sum',
                   action= 'store_true',
                   help='''Compute md5sum in addition to the other stats.
                   ''')

args = parser.parse_args()

# -----------------------------[ Functions ]-----------------------------------

def sumfile(fobj):
    '''Returns an md5 hash for an object with read() method.'''
    m = hashlib.md5()
    while True:
        d = fobj.read(8096)
        if not d:
            break
        m.update(d)
    return(m.hexdigest())

def md5sum(fname):
    '''Returns an md5 hash for file fname, or stdin if fname is "-".'''
    if fname == '-':
        ret = sumfile(sys.stdin)
    else:
        f = open(fname, 'rb')
        ret = sumfile(f)
        f.close()
    return(ret)
# -----------------------------------------------------------------------------

try:
    filename= os.path.split(args.infile)[1]
    path= os.path.abspath(os.path.split(args.infile)[0])
    if args.md5sum is True:
        md5sum= md5sum(args.infile)
    else:
        md5sum= ''
    hostname= socket.gethostname()
    ctime= datetime.datetime.fromtimestamp(os.path.getctime(args.infile)).isoformat()
    mtime= datetime.datetime.fromtimestamp(os.path.getmtime(args.infile)).isoformat()
    fsize= os.path.getsize(args.infile)
except:
    sys.exit()    

filestats= {'filename': filename,
            'path':     path,
            'md5sum':   md5sum,
            'hostname': hostname,
            'ctime':    ctime,
            'mtime':    mtime,
            'fsize':    fsize}

print(filestats)
sys.exit()