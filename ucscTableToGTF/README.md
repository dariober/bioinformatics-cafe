<!-- MarkdownTOC -->

- [Converter from UCSC genome browser table to GTF](#converter-from-ucsc-genome-browser-table-to-gtf)
    - [Example usage](#example-usage)

<!-- /MarkdownTOC -->

Converter from UCSC genome browser table to GTF
===============================================

Convert a genePred table from UCSC genome browser (e.g. refGene, vegaGene) to GTF. 
Features extracted from table and annotated into the GTF file are:
exon, intron, CDS, start_codon, stop_codon, TSS, TES, 5utr, 3utr.

Example of input: [refGene.txt.gz](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz)

This script is quite a bit redundant as [genePredToGtf](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToGtf) from
UCSC is very similar. `ucscTableToGTF.py` also adds introns, TSS and TES.

Example usage
-------------

Download and convert on the fly the VEGA genes table:

```
curl -o - -O http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/vegaGene.txt.gz | gunzip -c \
| ./ucscTableToGTF.py \
| sort -k1,1 -k4,4n -k5,5n > vega.gtf
```