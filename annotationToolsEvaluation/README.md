TOC
===

<!-- MarkdownTOC -->

- [Get test data](#get-test-data)
- [Run tools](#run-tools)
    - [Annovar](#annovar)
    - [VEP](#vep)
    - [VAGReNT](#vagrent)
    - [CAVA](#cava)
- [Features](#features)

<!-- /MarkdownTOC -->


We want to evaluate some popular tools for annotatng variants.

Get test data
==============

Collect some vcf files from 1000 genomes project. Get variants from various types of mutations SNV,
Indels, etc.

These are vcf files annotated by 1000 genome project. We don't need so many variants for testing but
let's start with the whole lot  and later compare some selected sites across tools. As the names
say, one file is for coding variants, the other non-coding. See also README files in their
directory.

```
ssh t1l1.tcrc.gla.ac.uk

cd /scratch/dberaldi/projects/20170303_annotationToolsEvaluation/data

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/functional_annotation/annotated_vcfs/ALL.wgs.integrated_phase1_release_v3_coding_annotation.20101123.snps_indels.sites.vcf.gz
## wc -l 514269

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/functional_annotation/annotated_vcfs/ALL.wgs.integrated_phase1_release_v3_noncoding_annotation_20120330.20101123.snps_indels_sv.sites.vcf.gz
```

It appears the noncoding file has the header sperated by spaces instead of TAB (!?) 
which causes problems later. So change it to TABS.

```
zcat ALL.wgs.integrated_phase1_release_v3_noncoding_annotation_20120330.20101123.snps_indels_sv.sites.vcf.gz \
| awk '{if($0 ~ "#CHROM ") {gsub(" +", "\t", $0)} print $0 }' > ALL.wgs.integrated_phase1_release_v3_noncoding_annotation_20120330.20101123.snps_indels_sv.sites.vcf
bgzip -f ALL.wgs.integrated_phase1_release_v3_noncoding_annotation_20120330.20101123.snps_indels_sv.sites.vcf
tabix -f -p vcf ALL.wgs.integrated_phase1_release_v3_noncoding_annotation_20120330.20101123.snps_indels_sv.sites.vcf.gz
```

Split by chrom

```
cd /scratch/dberaldi/projects/20170303_annotationToolsEvaluation/data/
chroms=`zcat ALL.wgs.integrated_phase1_release_v3_noncoding_annotation_20120330.20101123.snps_indels_sv.sites.vcf.gz | grep -v '^#' | cut -f 1 | uniq`
for x in $chroms
do
    echo $x
    tabix -h ALL.wgs.integrated_phase1_release_v3_noncoding_annotation_20120330.20101123.snps_indels_sv.sites.vcf.gz $x \
    | bgzip -f > ALL.wgs.integrated_phase1_release_v3_noncoding_annotation_20120330.20101123.snps_indels_sv.sites.${x}.vcf.gz &&
    tabix -f -p vcf ALL.wgs.integrated_phase1_release_v3_noncoding_annotation_20120330.20101123.snps_indels_sv.sites.${x}.vcf.gz
done


cd /scratch/dberaldi/projects/20170303_annotationToolsEvaluation/data/
chroms=`zcat ALL.wgs.integrated_phase1_release_v3_coding_annotation.20101123.snps_indels.sites.vcf.gz | grep -v '^#' | cut -f 1 | uniq`
for x in $chroms
do
    echo $x
    tabix -h ALL.wgs.integrated_phase1_release_v3_coding_annotation.20101123.snps_indels.sites.vcf.gz $x \
    | bgzip -f > ALL.wgs.integrated_phase1_release_v3_coding_annotation.20101123.snps_indels.sites.${x}.vcf.gz &&
    tabix -f -p vcf ALL.wgs.integrated_phase1_release_v3_coding_annotation.20101123.snps_indels.sites.${x}.vcf.gz
done
```


Run tools
=========

_NB_ Should allele frequency be considered in the annotation? 

Annovar
-------

* [ANNOVAR](http://annovar.openbioinformatics.org/en/latest/)

[annovar here](./annovar.md)

VEP
---

* [VEP](http://www.ensembl.org/info/docs/tools/vep/index.html)

[vep here](./vep.md)


VAGReNT
-------

[vagrent](./vagrent.md)

* [VAGReNT](https://github.com/cancerit/VAGrENT). Need to install git and perl since git is not on 
  the cluster and perl is outdated (required: 5.14+). Install also Bio::Perl and Bio::DB:HTS.

Protocols is [here](http://onlinelibrary.wiley.com/doi/10.1002/0471250953.bi1508s52/full)

[vagrent here](./vagrent.md)

CAVA
----

* [CAVA](http://www.well.ox.ac.uk/cava) _done_. Source in `~/applications/CAVA-1.2.0` binaries in `~/applications/CAVA/bin/` installed as:

```
cd ~/applications/CAVA-1.2.0/
./setup.sh ../CAVA
```

Features
========

* Informative log if something goes wrong

* Can select one representative transcript among isoforms (CAVA and vagrent yes)

* 
