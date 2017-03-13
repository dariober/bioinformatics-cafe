TOC
===

<!-- MarkdownTOC -->

- [Vagrent](#vagrent)

<!-- /MarkdownTOC -->


Vagrent
=======

Generate cache files

It would be good to use the ensembl release corresponding to gencode v7 (same as
1kg), which is ensembl release 62 (see below). However, it appears versions < 65
do not have enough annotation in the fasta header and they can't be used. So
let's go for the latest version for GRCh37 which we have matched CDS sequence and ucsc
gencode.

See [release history](https://www.gencodegenes.org/releases/)

The `CCDS2Sequence.current.txt` file for release 75 is here:

```
cd /scratch/dberaldi/projects/20170303_annotationToolsEvaluation/vagrent/cache

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
```

Annotate 

```
cd /scratch/dberaldi/projects/20170303_annotationToolsEvaluation/vagrent/

sbatch -J vagrent --mem=4000 --wrap "AnnotateVcf.pl -i ../data/ALL.wgs.integrated_phase1_release_v3_coding_annotation.20101123.snps_indels.sites.vcf.gz \
    -o codingVagrent.vcf \
    -cache cache/Homo_sapiens.GRCh37.${rel}.vagrent.cache.gz \
    -sp homo_sapiens \
    -as GRCh37 \
    -tabix"

bgzip -f codingVagrent.vcf && tabix -f codingVagrent.vcf.gz

sacct --format="CPUTime,Elapsed,MaxRSS,JobName" -j 15473 ## All coding
#   CPUTime    Elapsed     MaxRSS    JobName 
#---------- ---------- ---------- ---------- 
#1-00:37:20   01:32:20               vagrent 
#  01:32:20   01:32:20     41068K      batch 
```