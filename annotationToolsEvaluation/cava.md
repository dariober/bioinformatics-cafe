TOC
===

<!-- MarkdownTOC -->

- [CAVA](#cava)

<!-- /MarkdownTOC -->


CAVA
====

* Prepare database

```
cd /scratch/dberaldi/projects/20170303_annotationToolsEvaluation/cava

~/applications/CAVA-1.2.0/ensembl_prep.py -e 75 -o db/ensembl75
```

* Configuration file: Prepared from template in source dir

* Run

```
cd /scratch/dberaldi/projects/20170303_annotationToolsEvaluation/cava

sbatch -J cava --mem=4000 --wrap "~/applications/CAVA-1.2.0/cava.py -t 16 -c config.txt \
    -i ../data/ALL.wgs.integrated_phase1_release_v3_coding_annotation.20101123.snps_indels.sites.vcf.gz \
    -o codingCava"

bgzip -f codingCava.vcf && tabix -f codingCava.vcf.gz
zcat codingCava.vcf.gz | wc -l ## 514288 vs 514239

sacct --format="CPUTime,Elapsed,MaxRSS,jobname" -j 15472 ## All coding
#   CPUTime    Elapsed     MaxRSS    JobName 
#---------- ---------- ---------- ---------- 
#  04:39:28   00:17:28                  cava 
#  00:17:28   00:17:28    593800K      batch 
```

```
bcftools view -H codingCava.vcf.gz | wc -l ## 514239
```