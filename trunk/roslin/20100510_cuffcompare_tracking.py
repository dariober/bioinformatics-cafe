"""

Reformat the output .tracking from cuffcompare. Stack the two compared files one on the other and remove
expression data (this is alread in cufflinks output transcripts.gtf)
Note that currently only two files can be parsed

TCONS_00000001	XLOC_000001	ENSSSCG00000000005|ENSSSCT00000000006	i	q1:LPS.144792|LPS.144792.0|100|5.562231|1.020688|10.103774|1.618497	-
TCONS_00000002	XLOC_000001	ENSSSCG00000000005|ENSSSCT00000000006	e	q1:LPS.144795|LPS.144795.0|100|5.677086|1.662780|9.691392|2.323009	-
TCONS_00000003	XLOC_000001	ENSSSCG00000000005|ENSSSCT00000000006	c	q1:LPS.144798|LPS.144798.0|100|4.284900|1.045820|7.523979|1.870229	q2:CTRL.152491|CTRL.152491.0|100|1.382326|0.000000|3.337231|0.536398
TCONS_00000004	XLOC_000002	Q6QAP4_PIG|ENSSSCT00000000009	j	q1:LPS.144816|LPS.144816.0|100|29.603764|23.769873|35.437655|12.356631	-
TCONS_00000005	XLOC_000003	-        	u	q1:LPS.144849|LPS.144849.0|100|11.037749|7.690534|14.384965|4.014703	-
"""

infile_name= 'C:/Tritume/20100505_RNAseq_cuffcompare/RNAseq.tracking'
outfile_name= 'C:/Tritume/20100505_RNAseq_cuffcompare/RNAseq.tracking2'


def flatten(lst):
    " Returns a flat list from a list containing nested lists "
    for elem in lst:
        if type(elem) in (tuple, list):
            for i in flatten(elem):
                yield(i)
        else:
            yield(elem)

infile= open(infile_name, 'r')
outfile= open(outfile_name, 'w')
header= ['cuffcompare_transcript_id', 'cuffcompare_locus_id', 'ref_gene_id', 'ref_transcript_id', 'class_code', 'cufflinks_gene_id', 'cufflinks_transcript_id']
header= '\t'.join(header)
outfile.write(header + '\n')

for line in infile:
    line= line.split('\t')
    if line[2].strip() == '-':
        line[2]= ['', '']
    else:
        line[2]= line[2].split('|')
    if line[4] == '-':
        line[4]= ['']*7
    else:
        line[4]= line[4].replace('q1:', '')
        line[4]= line[4].split('|')
    if line[5] == '-\n':
        line[5]= ['']*7
    else:
        line[5]= line[5].replace('q2:', '')
        line[5]= line[5].split('|')
    line= list(flatten(line))

    if line[5] != '':
        line_1= line[0:7]
        line_1= '\t'.join(line_1)
        outfile.write(line_1 + '\n')

    if line[12] != '':
        line_2= line[0:5] + line[12:14]
        line_2= '\t'.join(line_2)
        outfile.write(line_2 + '\n')

infile.close()
outfile.close()
