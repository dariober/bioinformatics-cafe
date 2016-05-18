<!-- MarkdownTOC -->

- [Converter from UCSC genome browser table to GTF](#converter-from-ucsc-genome-browser-table-to-gtf)

<!-- /MarkdownTOC -->

Converter from UCSC genome browser table to GTF
===============================================

Convert a *gene table from UCSC genome browser (e.g. refGene, vegaGene) to GTF. 
Features extracted from table and annotated into the GTF file are:
exon, intron, CDS, start_codon, stop_codon, TSS, TES, 5utr, 3utr.

Example of input: [refGene.txt.gz](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz)

This script is quite a bit redundant as [genePredToGtf](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToGtf) from
UCSC is very similar. 