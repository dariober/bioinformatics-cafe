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

