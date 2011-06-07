#!/cygdrive/c/Python26/python

"""
  Add a file identifier and concatenate files. 
  081_mapping.tbl.txt.bz2 etc are files with mappings of CAGE tags from FANTOM4 (mouse)
  see http://fantom.gsc.riken.jp/4/download/Tables/mouse/CAGE/mapping/
"""
import bz2
import os
" ----------------------------[ Input/output ] ------------------------------ "

## Where files are:
wdir= 'F:/data/FANTOM4/Tables/mouse/CAGE/mapping'

## File names:
files= ['081_mapping.tbl.txt.bz2', '082_mapping.tbl.txt.bz2', '083_mapping.tbl.txt.bz2', '097_mapping.tbl.txt.bz2', '119_mapping.tbl.txt.bz2', '144_mapping.tbl.txt.bz2']

## Concatenated ouput file
outfile= open(wdir + '/con_mapping.tbl.txt', 'w')

## Header: make this the first line in the output:
header= '\t'.join(['dataset_id', 'id', 'library_count', 'edit_string', 'chrom', 'strand', 'start', 'end', 'percentage', 'map_pos', 'ribo_flag', 'refseq_flag', 'rna', 'rna_count', 'tpm_in_ribo', 'tpm_ex_ribo']) ## These columns left out: 'weight', 'rescue_weight'

## Identifiers are the file names stripped of extensions (e.g. 081_mapping.tbl.txt.bz2 -> 081_mapping)
## See function get_file_id()
" ---------------------------------------------------------------------------- "

os.chdir(wdir)

def get_file_id(filename):
    fid= filename.split('.')[0]
    return(fid)

outfile.write(header + '\n')

for file in files:
    file_id= get_file_id(file)
    fopen= bz2.BZ2File(file)
    n= 0
    print('Processing ' + file)
    for line in fopen:
        if line.startswith('#') or line.startswith('id ') or 'macrophage' not in line:
            " Skip header/comment lines and skip lines that have nothing to do with macrophages "
            continue
        else:
            line= line.strip()
            line= line.split('\t')[0:15]
            line= '\t'.join(line)
            line= file_id + '\t' + line
            n += 1
        outfile.write(line + '\n')
    fopen.close()
    print(str(n) + ' lines sent to output')
outfile.close()