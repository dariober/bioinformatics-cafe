TOC
===

<!-- MarkdownTOC -->

- [VEP](#vep)

<!-- /MarkdownTOC -->


VEP
===

Install appropriate annotation files (see [vep_cache](http://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html)).

```
cd /scratch/dberaldi/projects/20170303_annotationToolsEvaluation/vep/vep_cache

rel=75
wget ftp://ftp.ensembl.org/pub/release-${rel}/variation/VEP/homo_sapiens_vep_${rel}.tar.gz

tar xfz homo_sapiens_vep_${rel}.tar.gz
```
    

```
cd /scratch/dberaldi/projects/20170303_annotationToolsEvaluation/vep/

export PATH=$HOME/applications/vep/ensembl-tools-release-${rel}/scripts/variant_effect_predictor:$PATH

sbatch -J vep --mem=4000 --wrap "variant_effect_predictor.pl \
    --force \
    --dir ./vep_cache \
    --offline \
    --species homo_sapiens \
    --assembly GRCh37 \
    --cache \
    -cache_version ${rel} \
    --hgvs \
    --fasta /scratch/dberaldi/ref/human_g1k_v37.fasta \
    --vcf \
    -i ../data/ALL.wgs.integrated_phase1_release_v3_coding_annotation.20101123.snps_indels.sites.vcf.gz \
    -o codingVep.vcf"

bgzip -f codingVep.vcf && tabix -f codingVep.vcf.gz

sacct --format="CPUTime,Elapsed,MaxRSS,JobName" -j 15475 ## All coding
   CPUTime    Elapsed     MaxRSS    JobName 
---------- ---------- ---------- ---------- 
  12:10:56   00:45:41                   vep 
  00:45:41   00:45:41   1406912K      batch 
```