#!/usr/local/bin/python

import sys

if len(sys.argv) == 1:
    sys.exit("""
Convert FastQC report (fastqc_data.txt) to table format for importing to
postgres.
""")

"""
##FastQC	0.10.0
>>Basic Statistics	pass
#Measure	Value	
Filename	ds001_sample1.bam	
File type	Conventional base calls	
Encoding	Sanger / Illumina 1.9	
Total Sequences	999556	
Filtered Sequences	0	
Sequence length	36	
%GC	42	
>>END_MODULE
>>Per base sequence quality	pass
#Base	Mean	Median	Lower Quartile	Upper Quartile	10th Percentile	90th Percentile
1	36.34938412655219	38.0	35.0	39.0	31.0	40.0
2	35.648754046796775	38.0	35.0	39.0	28.0	40.0
3	35.502683191336956	38.0	34.0	39.0	28.0	40.0
...
>>END_MODULE
>>Per sequence quality scores	pass
#Quality	Count
2	3636.0
3	203.0
...
>>END_MODULE
>>Per base sequence content	warn
#Base	G	A	T	C
1	13.47777203256334	29.391581791144084	31.86158976763948	25.269056408653096
...
"""

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
                pgrow.append('{' + str(x).strip('[]') + '}')
            except ValueError:
                try:
                    x= [float(z) for z in x]
                    pgrow.append('{' + str(x).strip('[]') + '}')
                except:
                    pgrow.append('{' + str(x).strip('[]') + '}')
        else:
            pgrow.append(x)
    return(pgrow)
   
" ---------------------------------------------------------------------------- "

fq= open(sys.argv[1]).readlines()
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
    fastqc_line.append(line.split('\t')[1])

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
print('\t'.join(fastqc_line))
sys.exit()


