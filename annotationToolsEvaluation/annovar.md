TOC
===

<!-- MarkdownTOC -->

- [Annovar](#annovar)

<!-- /MarkdownTOC -->


Annovar
=======

First we need to collect the annotation data using `annotate_variation.pl -downdb ...`. Test files are in reference=GRCh37 
(see vcf header) and http://www.internationalgenome.org/faq/which-reference-assembly-do-you-use/

```
cd /scratch/dberaldi/projects/20170303_annotationToolsEvaluation/annovar

export PATH=$HOME/applications/annovar/annovar:$PATH

## Annotation can be obtained from UCSC tables, e.g. download table ncbiRefSeqCurated:
annotate_variation.pl -buildver hg38 -downdb ncbiRefSeqCurated humandb/

annotate_variation.pl -buildver hg19 -downdb refGene humandb/
annotate_variation.pl --buildver hg19 --downdb seq humandb/hg19_seq
retrieve_seq_from_fasta.pl humandb/hg19_refGene.txt -seqdir humandb/hg19_seq -format refGene -outfile humandb/hg19_refGeneMrna.fa

annotate_variation.pl -buildver hg19 -downdb wgEncodeGencodeBasicV19 humandb
annotate_variation.pl --buildver hg19 --downdb seq humandb/hg19_seq
retrieve_seq_from_fasta.pl humandb/hg19_wgEncodeGencodeBasicV19.txt -seqdir humandb/hg19_seq -format refGene -outfile humandb/hg19_wgEncodeGencodeBasicV19Mrna.fa

## See what annovar databases are available:
annotate_variation.pl -webfrom annovar -downdb avdblist -buildver hg38 .
cat hg38_avdblist.txt
hg38_avsnp142.txt.gz                20150106  1282852569
hg38_avsnp142.txt.idx.gz            20150106  212764014
hg38_avsnp144.txt.gz                20151102  1671400214
hg38_avsnp144.txt.idx.gz            20151102  215204030
hg38_avsnp147.txt.gz                20160601  1775686247
hg38_avsnp147.txt.idx.gz            20160601  222202148
...

## E.g. get also dbSNP
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar avsnp147 humandb/
```

Do the actual annotation

```
cd /scratch/dberaldi/projects/20170303_annotationToolsEvaluation/annovar/

## Only gencode
sbatch -J annovar --mem=4000 --wrap "table_annovar.pl ../data/ALL.wgs.integrated_phase1_release_v3_coding_annotation.20101123.snps_indels.sites.vcf.gz \
    humandb/ \
    -buildver hg19 \
    -out codingAnno \
    -remove \
    -protocol wgEncodeGencodeBasicV19 \
    -operation g \
    -nastring . \
    -vcfinput"

sacct --format="CPUTime,Elapsed,MaxRSS,JobName" -j 15471 ## All coding
#   CPUTime    Elapsed     MaxRSS    JobName 
#---------- ---------- ---------- ---------- 
#  00:37:20   00:02:20               annovar 
#  00:02:20   00:02:20   1869252K      batch 

bgzip -f codingAnno.hg19_multianno.vcf && tabix -f codingAnno.hg19_multianno.vcf.gz
```
