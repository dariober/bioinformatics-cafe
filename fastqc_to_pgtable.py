#!/usr/bin/python

import os
import sys
import psycopg2
import argparse

parser = argparse.ArgumentParser(description= """
DESCRIPTION
    Convert the FastQC report (file: fastqc_data.txt) to a single line with columns
    tab-separated. The output is uploaded to the postgres database given in the
    connection string. 
        
EXAMPLE
    fastqc_to_pgtable.py -i fastqc/fastqc_data.txt
    
    ## Process all fastqc_data.txt files in all subdirs in current dir:
    for fastqc in `ls .`
    do
    fastqc_to_pgtable.py -i $fastqc/fastqc_data.txt
    done

TODO:
    - Create a template fastqc_line (which holds the output of fastqc) to assign
      correct default values to each item.
    - Read zip files

""", formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--infile', '-i',
                   required= True,
                   help='''Fastqc output file containing the file stats as text.
Typically this is fastqc_data.txt
                   ''')

parser.add_argument('--md5sum',
                   required= False,
                   default= None,
                   help='''Optional md5sum of the fastq/bam file passed to fastqc
default is None
                   ''')

parser.add_argument('--fsize',
                   required= False,
                   default= None,
                   type= int,
                   help='''Optional file size (as int) of the fastq/bam file passed to fastqc
default is None
                   ''')

parser.add_argument('--mtime',
                   required= False,
                   default= None,
                   help='''Optional modification time of the fastq/bam file passed to fastqc
default is None. Typically this is returned by os.path.getmtime().
                   ''')

args = parser.parse_args()

# -----------------------------------------------------------------------------
def get_psycopgpass():
    """
    Read file ~/.psycopgpass to get the connection string to pass to
    psycopg2.connect()
    """
    conn= open(os.path.join(os.getenv("HOME"), '.psycopgpass'))
    conn= conn.readlines()
    conn_args= [x.strip() for x in conn if x.strip() != '' and x.strip().startswith('#') is False][0]
    return(conn_args)
    
def parse_module(fastqc_module):
    """
    Parse a fastqc module from the table format to a line format (list).
    Input is list containing the module. One list-item per line. E.g.:

    fastqc_module= [
        '>>Per base sequence quality	pass',
        '#Base	Mean	Median	Lower Quartile	Upper Quartile	10th Percentile	90th Percentile',
        '1	36.34	38.0	35.0	39.0	31.0	40.0',
        '2	35.64	38.0	35.0	39.0	28.0	40.0',
        '3	35.50	38.0	34.0	39.0	28.0	40.0',
        ...
        ]
    
    Return a list like this where each sublist after 1st is a column:
    ['pass', ['1', '2', '3', ...], ['36,34', '35.64', '35.50', ...], ['40.0', '40.0', '40.0', ...], ...]
    """
    row_list= []
    module_header= fastqc_module[0]
    module_name= module_header.split('\t')[0]
    row_list.append(module_header.split('\t')[1]) ## Line with module name >>Per base ...    pass/warn/fail    
    
    # Handle odd cases:
    # Where no table is returned:
    if len(fastqc_module) == 1 and module_name == '>>Overrepresented sequences':
        return(row_list + [[]]*4)
    if len(fastqc_module) == 1 and module_name == '>>Kmer Content':
        return(row_list + ['']*5)
    # Table is not the secod row:
    if module_name == '>>Sequence Duplication Levels':
        tot_dupl= fastqc_module[1].split('\t')[1]
        row_list.append(tot_dupl)
        del(fastqc_module[1])
   
    ## Conevrt table to list of lists:
    tbl= []
    for line in fastqc_module[2:]:
        tbl.append(line.split('\t'))
    ## Put each column in a list:
    nrows= len(tbl)
    ncols= len(tbl[0])
    for i in range(0, ncols):
        col= []
        for j in range(0, nrows):
            col.append(tbl[j][i])
        row_list.append(col)
    return(row_list)

def list_to_pgcolumns(lst):
    """
    Convert the nested lists in lst to strings formatted as array type suitable
    for postgres. E.g.
    ['a', [1,2,3]] -> ['a', '{1,2,3}']
    ['x', ['a', 'b', 'c']]-> ['x', "{'a', 'b', 'c'}"]
    """
    pgrow= []
    for x in lst:
        if type(x) == list:
            try:
                x= [int(z) for z in x]
                # pgrow.append('{' + str(x).strip('[]') + '}')
                pgrow.append(x)
            except ValueError:
                try:
                    x= [float(z) for z in x]
                    pgrow.append(x)
                    # pgrow.append('{' + str(x).strip('[]') + '}')
                except:
                    pgrow.append(x)
                    # pgrow.append('{' + str(x).strip('[]') + '}')
        else:
            pgrow.append(x)
    return(pgrow)

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
 
" ---------------------------------------------------------------------------- "

fq= open(args.infile).readlines()
fq= [x.strip() for x in fq]

fastqc_line= [] ## The entire report will be in this list. Each item is a column in the database table 

fastqc_line.append(fq[0].split('\t')[1]) ## Get Fastqc version
fastqc_line.append(fq[1].split('\t')[1]) ## base stats

## Get start and end position of all modules:
mod_start= []
mod_end= []
for i in range(0, len(fq)):
    line= fq[i]
    if line == '>>END_MODULE':
        mod_end.append(i)
    elif line.startswith('>>'):
        mod_start.append(i)
    else:
        pass

## Process First module:
module= fq[mod_start[0]:mod_end[0]]
for line in module[2:]:
    line= line.split('\t')
    if line[0] == 'Filename':
        filename= line[1]
    fastqc_line.append(line[1])

## Start processing modules. First one (Basic statitics) is apart:
for s, e in zip(mod_start[1:], mod_end[1:]):
    module= fq[s:e]
    module_name= module[0].split('\t')[0]
    row= parse_module(module)
    # For these modules remove the column with base position (1st column, 2nd item in list):
    if module_name in ['>>Per base sequence quality', '>>Per base sequence content', '>>Per base GC content', '>>Per base N content']:
        del(row[1])
    if module_name == '>>Sequence Duplication Levels':
        row[2][9]= '10' ## Replace 10++ with 10
    row= list_to_pgcolumns(row) 
    fastqc_line= fastqc_line + row

## Get md5sum, fsize and mtime of fastq/bam file
fastqc_line.append(args.md5sum)
fastqc_line.append(args.fsize)
fastqc_line.append(args.mtime)

# ---------------------[ Send to postres ]--------------------------------------
conn= psycopg2.connect(get_psycopgpass())
cur= conn.cursor()
## Memo: get columns with query like this:
## select * from information_schema.columns where table_name = 'fastqc' order by ordinal_position;
cur.execute("CREATE TEMP TABLE fastqc_tmp AS (SELECT * FROM fastqc WHERE 1=2)") ## Where to put results
cur.execute("INSERT INTO fastqc_tmp VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)", fastqc_line)
cur.execute("INSERT INTO fastqc SELECT DISTINCT * FROM fastqc_tmp EXCEPT SELECT * FROM fastqc") ## Import only rows not already present 
cur.execute("DROP TABLE fastqc_tmp")
cur.close()
conn.commit()
conn.close()
sys.exit()


