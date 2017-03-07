TOC
===

<!-- MarkdownTOC -->

- [Get test data](#get-test-data)
- [Tools](#tools)
- [Run on VCF files](#run-on-vcf-files)
    - [Annovar](#annovar)
    - [VEP](#vep)

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

Tools
=====

_NB_ Should allele frequency be considered in the annotation? 

* [CAVA](http://www.well.ox.ac.uk/cava) _done_. Source in `~/applications/CAVA-1.2.0` binaries in `~/applications/CAVA/bin/` installed as:

```
cd ~/applications/CAVA-1.2.0/
./setup.sh ../CAVA
```

* [ANNOVAR](http://annovar.openbioinformatics.org/en/latest/)

* [VEP](http://www.ensembl.org/info/docs/tools/vep/index.html)

Install appropriate annotation files (see [vep_cache](http://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html)). These files are massive (~700GB)!

```
cd /scratch/dberaldi/projects/20170303_annotationToolsEvaluation/vep/vep_cache

wget ftp://ftp.ensembl.org/pub/release-87/variation/VEP/homo_sapiens_vep_87_GRCh37.tar.gz

tar xfz homo_sapiens_vep_87.tar.gz
```

* [VAGReNT](https://github.com/cancerit/VAGrENT). Need to install git and perl since git is not on 
  the cluster and perl is outdated (required: 5.14+). Install also Bio::Perl and Bio::DB:HTS.

  There is no docs!!


Run on VCF files
================

Annovar
-------

[annovar here](./annovar.md)

VEP
---

[vep here](./vep.md)

