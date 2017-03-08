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
sbatch --mem=4000 --wrap "~/applications/CAVA-1.2.0/cava.py -t 16 -c config.txt \
    -i ../data/ALL.wgs.integrated_phase1_release_v3_coding_annotation.20101123.snps_indels.sites.vcf.gz \
    -o codingCava"


chroms=`zcat ../data/ALL.wgs.integrated_phase1_release_v3_noncoding_annotation_20120330.20101123.snps_indels_sv.sites.vcf.gz | grep -v '^#' | cut -f 1 | uniq`

for x in $chroms
do
sbatch --mem=4000 --wrap "~/applications/CAVA-1.2.0/cava.py -t 16 -c config.txt \
    -i ../data/ALL.wgs.integrated_phase1_release_v3_noncoding_annotation_20120330.20101123.snps_indels_sv.sites.${x}.vcf.gz \
    -o noncodingCava.${x}"
done

sacct --format="CPUTime,elapsed,MaxRSS" -j 15438 ## chr2 non-coding
#   CPUTime    Elapsed     MaxRSS 
#---------- ---------- ---------- 
#  00:58:08   00:03:38            
#  00:03:38   00:03:38    683532K 

sacct --format="CPUTime,MaxRSS" -j 15435 ## All coding
#   CPUTime     MaxRSS 
#---------- ---------- 
#  04:39:28            
#  00:17:28    593520K 

bgzip codingCava.vcf && tabix codingCava.vcf.gz
zcat codingCava.vcf.gz | wc -l ## 514288 vs 514239

## concat vcf
for vcf in noncoding*.vcf
do
    echo $vcf
    bgzip $vcf && tabix ${vcf}.gz
done

## Check log files for rekected variants.
bcftools concat noncodingCava.*.vcf.gz | bgzip -f > noncodingCava.vcf.gz && tabix -f noncodingCava.vcf.gz
zcat noncodingCava.vcf.gz | grep -v '^#' | wc -l ## 9445132 vs 9450275!! 
zcat ../data/ALL.wgs.integrated_phase1_release_v3_noncoding_annotation_20120330.20101123.snps_indels_sv.sites.vcf.gz | grep -v '^#' | wc -l ## 9450275
rm noncodingCava.*.vcf.gz*

```

