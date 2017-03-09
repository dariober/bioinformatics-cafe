TOC
===

<!-- MarkdownTOC -->

- [Annovar](#annovar)
    - [Comparing output from annovar with the one from 1kg](#comparing-output-from-annovar-with-the-one-from-1kg)

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

sacct --format="CPUTime,Elapsed,MaxRSS,JobName" -j 15471
#   CPUTime    Elapsed     MaxRSS    JobName 
#---------- ---------- ---------- ---------- 
#  00:37:20   00:02:20               annovar 
#  00:02:20   00:02:20   1869252K      batch 


sbatch --mem=12000 --wrap "table_annovar.pl ../data/ALL.wgs.integrated_phase1_release_v3_noncoding_annotation_20120330.20101123.snps_indels_sv.sites.vcf.gz \
    humandb/ \
    -buildver hg19 \
    -out noncodingAnno \
    -remove \
    -protocol wgEncodeGencodeBasicV19 \
    -operation g \
    -nastring . \
    -vcfinput"

sacct --format="CPUTime,MaxRSS,elapsed" -j 15388
   CPUTime     MaxRSS    Elapsed 
---------- ---------- ---------- 
  05:30:24              00:20:39 
  00:20:39   6275364K   00:20:39 
  
bgzip noncodingAnno.hg19_multianno.vcf && tabix -p vcf noncodingAnno.hg19_multianno.vcf.gz
bgzip -f codingAnno.hg19_multianno.vcf && tabix -f codingAnno.hg19_multianno.vcf.gz
```

Comparing output from annovar with the one from 1kg
---------------------------------------------------

VCF from 1kg has annotation in VA field. For SNP at 1:13302, rs180734498. It appears this filed is delimited by colon:

```
VA= 1
    AL627309.2
    ENSG00000249291.2
    +
    synonymous
    1/1
    AL627309.2-201
    ENST00000515242.2
    384_225_75_H->H;
```

This is the annotation from annovar for the same locus. This is appended to each record in the vcf:

```
ANNOVAR_DATE=2016-02-01
Func.wgEncodeGencodeBasicV7=exonic
Gene.wgEncodeGencodeBasicV7=AL627309.2
GeneDetail.wgEncodeGencodeBasicV7=.
ExonicFunc.wgEncodeGencodeBasicV7=synonymous_SNV
AAChange.wgEncodeGencodeBasicV7=AL627309.2
    ENST00000515242.2
    exon3
    c.C225T
    p.H75H
ALLELE_END
```

Let's have a look at an indel. 

Note useful page explaining how mutations are reported http://www.hgmd.cf.ac.uk/docs/mut_nom.html

```
bioawk -c vcf 'length($ref) > 1 || length($alt) > 1 {print $0}' codingAnno.hg19_multianno.vcf | less

1       1577084 .       TTCCTCCTCC      T       505     PASS    THETA=0.0021;AVGPOST=0.6147;AN=2184;VT=INDEL;LDAF=0.3179;ERATE=0.0502;AC=362;RSQ=0.2402;AF=0.17;ASN_AF=0.16;AMR_AF=0.19;AFR_AF=0.09;EUR_AF=0.21;VA=1:CDK11A:ENSG00000008128.14:-:deletionNFS:4/14:CDK11A-201:ENST00000317673.7:2352_951_317_EEEE->E:CDK11A-202:ENST00000340677.5:2319_918_306_EEEE->E:CDK11A-208:ENST00000407249.3:2358_957_319_EEEE->E:CDK11A-204:ENST00000341832.6:2217_816_272_EEEE->E;ANNOVAR_DATE=2016-02-01;Func.wgEncodeGencodeBasicV7=exonic;Gene.wgEncodeGencodeBasicV7=CDK11A;GeneDetail.wgEncodeGencodeBasicV7=.;ExonicFunc.wgEncodeGencodeBasicV7=nonframeshift_deletion;AAChange.wgEncodeGencodeBasicV7=CDK11A:ENST00000340677.5:exon9:c.909_917del:p.303_306del,CDK11A:ENST00000341832.6:exon9:c.807_815del:p.269_272del,CDK11A:ENST00000317673.7:exon10:c.942_950del:p.314_317del,CDK11A:ENST00000407249.3:exon10:c.948_956del:p.316_319del;ALLELE_END
```

From 1kg:

```
VA=1
    CDK11A
    ENSG00000008128.14
    -
    deletionNFS
    4/14
    CDK11A-201
    ENST00000317673.7
    2352_951_317_EEEE->E
    CDK11A-202
    ENST00000340677.5
    2319_918_306_EEEE->E
    CDK11A-208
    ENST00000407249.3
    2358_957_319_EEEE->E
    CDK11A-204
    ENST00000341832.6
    2217_816_272_EEEE->E
```

From Annovar

```
ANNOVAR_DATE=2016-02-01
Func.wgEncodeGencodeBasicV7=exonic
Gene.wgEncodeGencodeBasicV7=CDK11A
GeneDetail.wgEncodeGencodeBasicV7=.
ExonicFunc.wgEncodeGencodeBasicV7=nonframeshift_deletion
AAChange.wgEncodeGencodeBasicV7=CDK11A  <- gene name
    ENST00000340677.5                   <- transcript affected
    exon9                               <- exon number
    c.909_917del                        <- 'c.' = cDNA
    p.303_306del,CDK11A                 <- 'p.' = protein
    ENST00000341832.6
    exon9
    c.807_815del
    p.269_272del,CDK11A
    ENST00000317673.7
    exon10
    c.942_950del
    p.314_317del,CDK11A
    ENST00000407249.3
    exon10
    c.948_956del
    p.316_319del
ALLELE_END
```

A splicing variant:

```
1       886620  .       T       A       100     PASS    AA=T;AN=2184;RSQ=0.5826;VT=SNP;THETA=0.0006;SNPSOURCE=EXOME;AC=1;ERATE=0.0003;AVGPOST=0.9992;LDAF=0.0008;AF=0.0005;EUR_AF=0.0013;VA=1:NOC2L:ENSG00000188976.6:-:spliceOverlap:1/1:NOC2L-001:ENST00000327044.6:2250;ANNOVAR_DATE=2016-02-01;Func.wgEncodeGencodeBasicV7=splicing;Gene.wgEncodeGencodeBasicV7=NOC2L;GeneDetail.wgEncodeGencodeBasicV7=ENST00000327044.6:exon12:c.1332-2A>T;ExonicFunc.wgEncodeGencodeBasicV7=.;AAChange.wgEncodeGencodeBasicV7=.;ALLELE_END

VA=1
    NOC2L
    ENSG00000188976.6
    -
    spliceOverlap
    1/1
    NOC2L-001
    ENST00000327044.6
    2250

ANNOVAR_DATE=2016-02-01
Func.wgEncodeGencodeBasicV7=splicing
Gene.wgEncodeGencodeBasicV7=NOC2L
GeneDetail.wgEncodeGencodeBasicV7=ENST00000327044.6:exon12:c.1332-2A>T
ExonicFunc.wgEncodeGencodeBasicV7=.
AAChange.wgEncodeGencodeBasicV7=.
ALLELE_END
```

Frameshift 

```
1       874816  .       C       CT      95      PASS    AVGPOST=0.9790;AN=2184;VT=INDEL;THETA=0.0033;ERATE=0.0043;AC=33;RSQ=0.6282;LDAF=0.0236;AF=0.02;ASN_AF=0.0035;AMR_AF=0.02;AFR_AF=0.03;EUR_AF=0.01;VA=1:SAMD11:ENSG00000187634.6:+:insertionFS:1/1:SAMD11-010:ENST00000342066.3:2046_682;ANNOVAR_DATE=2016-02-01;Func.wgEncodeGencodeBasicV7=exonic;Gene.wgEncodeGencodeBasicV7=SAMD11;GeneDetail.wgEncodeGencodeBasicV7=.;ExonicFunc.wgEncodeGencodeBasicV7=frameshift_insertion;AAChange.wgEncodeGencodeBasicV7=SAMD11:ENST00000342066.3:exon7:c.682_683insT:p.P228fs;ALLELE_END

VA=1
SAMD11
ENSG00000187634.6
+
insertionFS
1/1
SAMD11-010
ENST00000342066.3
2046_682

ANNOVAR_DATE=2016-02-01
Func.wgEncodeGencodeBasicV7=exonic
Gene.wgEncodeGencodeBasicV7=SAMD11
GeneDetail.wgEncodeGencodeBasicV7=.
ExonicFunc.wgEncodeGencodeBasicV7=frameshift_insertion
AAChange.wgEncodeGencodeBasicV7=SAMD11:ENST00000342066.3:exon7:c.682_683insT:p.P228fs
ALLELE_END
```
