TOC
===

<!-- MarkdownTOC -->

- [Vagrent](#vagrent)

<!-- /MarkdownTOC -->


Vagrent
=======

* Generate cache files

It would be good to use the ensembl release corresponding to gencode v7 (same as
1kg), which is ensembl release 62 (see below). However, it appears versions < 65
do not have enough annotation in the fasta header and they can't be used. So
let's go for the latest version for GRCh37 which we have matched CDS sequence and ucsc
gencode.

See [release history](https://www.gencodegenes.org/releases/)

The `CCDS2Sequence.current.txt` file for release 75 is here:

```
cd /scratch/dberaldi/projects/20170303_annotationToolsEvaluation/vagrent/cache

curl -s ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/archive/16/BuildInfo.current.txt | column -t
#tax_id  ncbi_release_number  ensembl_release_number  assembly_name  assembly_id       ccds_release_number  date_made_public
10090    104                  75                      GRCm38.p2      GCF_000001635.22  16                   20140407

rel=75
curl ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/archive/16/CCDS2Sequence.current.txt > CCDS2Sequence.${rel}.txt
```

Make sure version numbers are consistent with the `-ccds`, `-databse` and `-ftp`.

```
cd /scratch/dberaldi/projects/20170303_annotationToolsEvaluation/vagrent/

sbatch --mem=4000 --wrap "Admin_EnsemblReferenceFileGenerator.pl \
    -species human \
    -assembly GRCh37 \
    -database homo_sapiens_core_${rel}_37 \
    -ccds cache/CCDS2Sequence.${rel}.txt \
    -ftp ftp://ftp.ensembl.org/pub/release-${rel}/fasta/homo_sapiens/cdna/ \
    -output cache/"

sacct --format="CPUTime,MaxRSS" -j 15307
   CPUTime     MaxRSS 
---------- ---------- 
  10:53:36            
  00:40:51    545480K 
```

* Annotate coding and non-coding

```
cd /scratch/dberaldi/projects/20170303_annotationToolsEvaluation/vagrent/

chroms=`zcat ../data/ALL.wgs.integrated_phase1_release_v3_noncoding_annotation_20120330.20101123.snps_indels_sv.sites.vcf.gz | grep -v '^#' | cut -f 1 | uniq`

for x in $chroms
do
sbatch --mem=4000 --wrap "AnnotateVcf.pl -i ../data/ALL.wgs.integrated_phase1_release_v3_noncoding_annotation_20120330.20101123.snps_indels_sv.sites.${x}.vcf.gz \
    -o noncodingVagrent.$x.vcf \
    -cache cache/Homo_sapiens.GRCh37.${rel}.vagrent.cache.gz \
    -sp homo_sapiens \
    -as GRCh37 \
    -tabix"

sbatch --mem=4000 --wrap "AnnotateVcf.pl -i ../data/ALL.wgs.integrated_phase1_release_v3_coding_annotation.20101123.snps_indels.sites.${x}.vcf.gz \
    -o codingVagrent.$x.vcf \
    -cache cache/Homo_sapiens.GRCh37.${rel}.vagrent.cache.gz \
    -sp homo_sapiens \
    -as GRCh37 \
    -tabix"
done

sacct --format="CPUTime,MaxRSS" -j 15339 ## This is for chr2 non-coding (~800k variants)
   CPUTime     MaxRSS 
---------- ---------- 
  18:24:00            
  01:09:00     52716K 
```

Once done concat files

```
bcftools concat codingVagrent.*.vcf.gz | bgzip -f > codingVagrent.vcf.gz && tabix -f -p vcf codingVagrent.vcf.gz
zcat codingVagrent.vcf.gz | grep -v '^#' | wc -l ## 514239
zcat ../data/ALL.wgs.integrated_phase1_release_v3_coding_annotation.20101123.snps_indels.sites.vcf.gz | grep -v '^#' | wc -l ## 514239
rm codingVagrent.*.vcf.gz*

bcftools concat noncodingVagrent.*.vcf.gz | bgzip -f > noncodingVagrent.vcf.gz && tabix -f -p vcf noncodingVagrent.vcf.gz
zcat noncodingVagrent.vcf.gz | grep -v '^#' | wc -l ## 9450275
zcat ../data/ALL.wgs.integrated_phase1_release_v3_noncoding_annotation_20120330.20101123.snps_indels_sv.sites.vcf.gz | grep -v '^#' | wc -l ## 9450275
rm noncodingVagrent.*.vcf.gz*
```