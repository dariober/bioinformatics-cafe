#!/usr/bin/env python

"""
Get gene positions from GTF file
E.g.
/lustre/sblab/berald01/reference_data/genomes/iGenomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf

chr1    unknown exon    11874   12227   .       +       .       gene_id "DDX11L1"; transcript_id "NR_046018_1"; gene_name "DDX11L1"; tss_id "TSS14523";
chr1    unknown exon    12613   12721   .       +       .       gene_id "DDX11L1"; transcript_id "NR_046018_1"; gene_name "DDX11L1"; tss_id "TSS14523";
chr1    unknown exon    13221   14408   .       +       .       gene_id "DDX11L1"; transcript_id "NR_046018_1"; gene_name "DDX11L1"; tss_id "TSS14523";
chr1    unknown exon    14362   14829   .       -       .       gene_id "WASH7P"; transcript_id "NR_024540"; gene_name "WASH7P"; tss_id "TSS7359";
chr1    unknown exon    14970   15038   .       -       .       gene_id "WASH7P"; transcript_id "NR_024540"; gene_name "WASH7P"; tss_id "TSS7359";
chr1    unknown exon    15796   15947   .       -       .       gene_id "WASH7P"; transcript_id "NR_024540"; gene_name "WASH7P"; tss_id "TSS7359";
chr1    unknown exon    16607   16765   .       -       .       gene_id "WASH7P"; transcript_id "NR_024540"; gene_name "WASH7P"; tss_id "TSS7359";
chr1    unknown exon    16858   17055   .       -       .       gene_id "WASH7P"; transcript_id "NR_024540"; gene_name "WASH7P"; tss_id "TSS7359";

"""

import sys

class GtfFeature(object):
    """ """
    def __init__():
       self.seqname= str
       self.source= str
       self.feature= str        
    
    def parse_strig(self):
        line= self.split('\t')
        self.seqname= line[0]
    
def gtfAttr(gtf_attribute):
    """Given the attribute string from a gtf file,
    return a dictionary of {'feature':'value'}
    
    E.g:
    gtfAttr('''gene_id "DDX11L1"; transcript_id "NR_046018_1"; gene_name "DDX11L1"; tss_id "TSS14523";''')
    >>> {'gene_id':'DDX11L1', 'gene_name':'DDX11L1', ...}
    
    Note:
        It will probably dislike attributes containing double-quotes and semicolons in the feature or value.
    
    See also:
        gtf format specs at http://www.sanger.ac.uk/resources/software/gff/spec.html
    """
    gtf_attr= gtf_attribute.strip().strip(';').split(';')
    gtf_attr= [x.strip() for x in gtf_attr]
    gtf_dict= {}
    for x in gtf_attr:
        kval= x.split(' ')
        kval= [z.strip('"') for z in kval]
        feature, value= kval
        if feature in gtf_dict:
            sys.exit('Invalid gtf attribute: %s\nSome features appera more than once' %(gtf_attribute))
        else:
            gtf_dict[feature]= value
    return(gtf_dict)
    
gtf= open(sys.argv[1])
gene_dict= {} ## Dictionary of genes and exon positions
for line in gtf:
    if line.startswith('#'):
        continue
    x= GtfFeature(line)
    
    line= line.strip().split('\t')
    fstart= int(line[3])
    fend= int(line[4])
    attr= gtfAttr(line[-1])
    if attr['gene_id'] not in gene_dict:
        gene_dict
    break

sys.exit()



