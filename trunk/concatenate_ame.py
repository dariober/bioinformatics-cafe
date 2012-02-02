#!/usr/local/bin/python

"""
Walks through a dir containing the output of meme-chip (each a subdir),
concatenate the output from the AME module and returns a cross-tabulated file of
motifs (rows) and libraries (columns) with values being the corrected p-value for enrichment.

Like this:

motif   ds001.macs_combined     ds002.macs_combined     ds003.macs_combined     ds005.macs_combined     ds006.macs_combined     ds007.macs_combined     ds008.macs_combined     ds011.macs_combined
MA0026.1 Eip74EF        NA      NA      NA      4.619e-10       4.273e-10       9.624e-07       3.464e-08       NA
MA0028.1 ELK1   NA      NA      NA      NA      0.0004774       9.01e-05        4.526e-05       NA
MA0030.1 FOXF2  1.083e-11       2.182e-11       NA      NA      NA      NA      NA      NA
MA0031.1 FOXD1  1.324e-05       4.01e-07        NA      NA      NA      NA      NA      NA


USAGE

concatenate_ame.py <meme-dir> <output>

TODO:
    - Options to switch between normlized and cross-tab output.
    - Value to output in crosstab from (e.g. p-value or 0/1 for presence/absence)

MEMO:
Input looks like this:

meme                     <- Top mem-dir
|-- ds001.macs_combined
|   |-- ame_out
|       |-- ame.html
|       `-- ame.txt      <- This is the file searched and concatenated  
|   |-- dreme_ama_out
|   |-- dreme_mast_out
|   |-- dreme_out
|   |-- dreme_tomtom_out
|   |-- meme_ama_out
|   |-- meme_mast_out
|   |-- meme_out
|   `-- meme_tomtom_out
|-- ds002.macs_combined
|   |-- ame_out
|       |-- ame.html
|       `-- ame.txt      <- This is the file searched and concatenated  
|   |-- dreme_ama_out
|   |-- dreme_mast_out
|   |-- dreme_out
|   |-- dreme_tomtom_out
|   |-- meme_ama_out
|   |-- meme_mast_out
|   |-- meme_out
|   `-- meme_tomtom_out

                        ame.txt looks like this:
-----------------------------[ cut here ]---------------------------------------
AME (Analysis of Motif Enrichment): Compiled on Nov 20 2011 at 17:47:34
------------------------------
Copyright (c) Robert McLeay <r.mcleay@imb.uq.edu.au> & Timothy Bailey <t.bailey@imb.uq.edu.au>, 2009.

Command line
ame --verbose 1 --oc meme/ds005.macs_combined/ame_out --fix-partition 986 --bgformat 0 meme/ds005.macs_combined/seqs-centered_w_bg /home/berald01/applications/meme/db/motif_databases/JASPAR_CORE_2009.meme

Warning: Not in partition maximisation mode. Fixing partition at 986.

Threshold p-value for reporting results: 0.001
Motif p-values are corrected by #Motifs * #ThresholdsTested - (476 * 1 = 476)

1. Fisher-exact p-value of motif MA0316.1 HAP5 top 986 seqs: 4.622e-22 (Corrected p-value: 2.2e-19)
2. Fisher-exact p-value of motif MA0060.1 NFYA top 986 seqs: 1.703e-21 (Corrected p-value: 8.106e-19)
3. Fisher-exact p-value of motif MA0314.1 HAP3 top 986 seqs: 1.044e-17 (Corrected p-value: 4.969e-15)
-----------------------------[ cut here ]---------------------------------------

"""

import sys
import os
import re
def ame_line_to_tsv(ameline):
    """
    Parse a line of ame output to have it tab separated.
    Input and output is string (the ame line)
    """
    ameline= ameline.replace(' Fisher-exact p-value of motif ',  '\tFisher-exact p-value of motif\t')
    ameline= ameline.replace(' top ', '\ttop')
    ameline= ameline.replace(' seqs: ', ' seqs:\t')
    ameline=ameline.replace(' (Corrected p-value: ',  '\tCorrected p-value:\t')
    ameline= re.sub('\)$', '', ameline)
    return(ameline)

def crosstab(normtab, rowheader_index, colheader_index, value_index, missing= '', first_column_name= 'rowheader'):
    """
    Reshape a nomalized table given as list of lists (each inner list is a row) to
    a cross-tab format having as rows the items in column rowheader, columns are
    the elements in colheader and populated by the elements in values.
    
    normtab:
        List of lists
    rowheader_index, colheader_index, value_index:
        Index of the relevant columns
    missing:
        A value for missing variables 

    EXAMPLE:
    
    normtab=  [['a', 'x', 1],
               ['a', 'y', 0],
               ['a', 'z', 2],
               ['b', 'x', 3],
               ['b', 'z', 4]]
    ct= crosstab(testdata, 0, 1, 2, first_column_name= 'F')
    ct
    [['F',  'x', 'y', 'z'],
     ['a',   1,   0,   2 ],
     ['b',   3,  '',   4 ]]

    """
    rowheaders= sorted(set([x[rowheader_index] for x in normtab]))
    colheaders= sorted(set([x[colheader_index] for x in normtab]))
    ## Dictionary of dictionaries as:
    ## ddict= {'a' : {'x':1, 'y':0,  'z':2},
    ##               {'b':3, 'y':'', 'z':4}}
    ddict= {}

    for i in range(0, len(normtab)):
        """Traverse normalized table and get all the 'rowheaders', put them in
        the dictionary. At the end you get:
        ddict= {'a': {'y': 0, 'x': 1, 'z': 2}, 'b': {'x': 3, 'z': 4}}
        """
        row= normtab[i]
        rowheader= row[rowheader_index] 
        if rowheader not in ddict:
            ddict[rowheader]= {}
        "Get the column name and value for this row"
        colheader= row[colheader_index]
        value= row[value_index]
        if colheader not in ddict[rowheader]:
            ddict[rowheader][colheader]= value
    
    "Produce list of lists of the cross-tab"
    crosstab= [[first_column_name] + colheaders]
    for r in rowheaders:
        row= [r]
        dictrow= ddict[r]
        for c in colheaders:
            if c in dictrow:
                row.append(dictrow[c])
            else:
                row.append(missing)
        crosstab.append(row)
    return(crosstab)
    
subdirs= sorted(os.listdir(sys.argv[1]))

amelist= []
for d in subdirs:
    if os.path.isdir(os.path.join(sys.argv[1], d)):
        ameout= os.path.join(sys.argv[1], d, 'ame_out', 'ame.txt')
        try:
            ametxt= open(ameout).readlines()
        except:
            continue
        ametxt= [x.rstrip('\n\r') for x in ametxt]
        "Get last empty line, after that the output begins"
        amestart= len(ametxt) - ametxt[::-1].index('')
        for i in range(amestart, len(ametxt)):
            ameline= ame_line_to_tsv(ametxt[i])
            ameline= ameline.split('\t') + [d]
            amelist.append(ameline)

amecrosstab= crosstab(amelist, rowheader_index= 2, colheader_index= 7, value_index= 6, missing= 'NA', first_column_name= 'motif')
fout= open(sys.argv[2], 'w')
for line in amecrosstab:
    line= [str(x) for x in line]
    fout.write('\t'.join(line) + '\n')
fout.close()
sys.exit()