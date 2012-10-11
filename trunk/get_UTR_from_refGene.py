#!/usr/bin/env python

import sys
import os
import argparse
import operator

from interval import Interval, IntervalSet ## See http://pypi.python.org/pypi/interval/1.0.0 and http://stackoverflow.com/questions/6462272/subtract-overlaps-between-two-ranges-without-sets


parser = argparse.ArgumentParser(description= """

DESCRIPTION:
    Extract UTRs from UCSC table refGene.
    Output is sorted by chrom, start, end; without header line
    
    Non-coding transcripts are classified as completely 5'UTR if on + strand,
    3'UTR if on -.
        
    refGene file can be produce with
mysql --user=genome --host=genome-mysql.cse.ucsc.edu --disable-auto-rehash -e \
    "select * from hg19.refGene order by chrom, txStart, cdsStart, cdsEnd, txEnd;" > hg19.refgene.txt

EXAMPLE:
    get_UTR_from_refGene.py -i hg19.refgene.txt > hg19.refgene.utr.bed
    ## Separate 3'- from 5'-UTRs
    awk '$4 ~ /_3utr$/ {print $0}' hg19.refgene.utr.bed > hg19.refgene.3utr.bed
    awk '$4 ~ /_5utr$/ {print $0}' hg19.refgene.utr.bed > hg19.refgene.5utr.bed
    
""", formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('-i', '--input',
                    type= str,
                    required= True,
                    help="""RefGene file from UCSC 
                    
                    """)

args = parser.parse_args()
# -----------------------------------------------------------------------------

def sort_table(table, cols):
    """ sort a table by multiple columns
        table: a list of lists (or tuple of tuples) where each inner list 
               represents a row
        cols:  a list (or tuple) specifying the column numbers to sort by
               e.g. (1,0) would sort by column 1, then by column 0
    Found http://www.saltycrane.com/blog/2007/12/how-to-sort-table-by-columns-in-python/
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
#        print row

# -----------------------------------------------------------------------------
fin= open(args.input)

utr_list= []
header= False
for line in fin:
    line= line.strip().split('\t')
    if header is False:
        chromIdx= line.index('chrom')
        txStartIdx= line.index('txStart')
        txEndIdx= line.index('txEnd')
        cdsStartIdx= line.index('cdsStart')
        cdsEndIdx= line.index('cdsEnd')
        exStartsIdx= line.index('exonStarts')
        exEndsIdx= line.index('exonEnds')
        strandIdx= line.index('strand')
        nameIdx= line.index('name')
        name2Idx= line.index('name2')
        header= True
        continue
    strand= line[strandIdx]
    cdsStart= int(line[cdsStartIdx])
    cdsEnd= int(line[cdsEndIdx])
#    if cdsStart == cdsEnd:
#         "No coding sequence in this transcripts, skip it"
#        continue
    txStart= int(line[txStartIdx])
    txEnd= int(line[txEndIdx])
    exonStarts= line[exStartsIdx].strip(',').split(',')
    exonEnds= line[exEndsIdx].strip(',').split(',')
    exon_intervals= []
    for s,e in zip(exonStarts, exonEnds):
        exon_intervals.append(Interval(int(s), int(e)))
    exon_set= IntervalSet(exon_intervals) ## EXAMPLE= exon_set= IntervalSet([Interval(1, 50), Interval(100, 200), Interval(800, 1300)])

    txStart_cdsStart= IntervalSet([Interval(txStart, cdsStart)]) ## EXAMPLE= IntervalSet([Interval(1, 1000, lower_closed=True, upper_closed=True)])
    txEnd_cdsEnd= IntervalSet([Interval(cdsEnd, txEnd)])
    if strand == '+':
        utr5= txStart_cdsStart - (txStart_cdsStart - exon_set)
        utr3= txEnd_cdsEnd - (txEnd_cdsEnd - exon_set)
    if strand == '-':
        utr3= txStart_cdsStart - (txStart_cdsStart - exon_set)
        utr5= txEnd_cdsEnd - (txEnd_cdsEnd - exon_set)
    for x in utr5:
        if type(x) == int:
            "Skip transcripts which don't have UTR"
            continue
        name5utr= '_'.join([line[chromIdx], str(x.lower_bound), str(x.upper_bound), line[nameIdx], '5utr'])
        utr_list.append(([line[chromIdx], x.lower_bound, x.upper_bound, name5utr, line[nameIdx], line[strandIdx], line[name2Idx]]))
    for x in utr3:
        if type(x) == int:
            continue
        name3utr= '_'.join([line[chromIdx], str(x.lower_bound), str(x.upper_bound), line[nameIdx], '3utr'])
        utr_list.append([line[chromIdx], x.lower_bound, x.upper_bound, name3utr, line[nameIdx], line[strandIdx], line[name2Idx]])
fin.close()
utr_list_sorted= sort_table(utr_list, [0,1,2])
for line in utr_list_sorted:
    print('\t'.join([str(x) for x in line]))
sys.exit()
