#!/usr/local/bin/python

import sys
import os
import socket
import hashlib
import datetime

if len(sys.argv) != 2:
    sys.exit("""
    Get file stats for the file passed as first positional argument.
    Returns a string to stdout formatted like a dict:

    USAGE:
        get_file_stats.py myfile.txt
    
    OUTPUT:
    
        Printed on a single line:
        "{'filename': os.path.split(sys.argv[1])[1],
          'path':     os.path.split(sys.argv[1])[0],
          'md5sum':   md5sum(sys.argv[1]),
          'hostname': socket.gethostname(),
          'ctime':    datetime.datetime.fromtimestamp(os.path.getctime()),
          'mtime':    os.path.getmtime(),
          'fsize':     os.path.getsize()}"

    IMPORTANT: If for any reasons any of the file stats fail (e.g. the file doesn't
    exist), get_file_stats.py exits and prints nothing, no error is returned
         """)

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

try:
    filename= os.path.split(sys.argv[1])[1]
    path= os.path.abspath(os.path.split(sys.argv[1])[0])
    md5sum= md5sum(sys.argv[1])
    hostname= socket.gethostname()
    ctime= datetime.datetime.fromtimestamp(os.path.getctime(sys.argv[1]))
    mtime= datetime.datetime.fromtimestamp(os.path.getmtime(sys.argv[1]))
    fsize= os.path.getsize(sys.argv[1])
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