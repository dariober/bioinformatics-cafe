TOC
===

<!-- MarkdownTOC -->

- [VEP](#vep)

<!-- /MarkdownTOC -->


VEP
===

```
cd /scratch/dberaldi/projects/20170303_annotationToolsEvaluation/vep/

export PATH=$HOME/applications/vep/ensembl-tools-release-87/scripts/variant_effect_predictor:$PATH

sbatch --mem=4000 --wrap "variant_effect_predictor.pl \
    --force \
    --dir ./vep_cache \
    --offline \
    --species homo_sapiens \
    --assembly GRCh37 \
    --cache \
    -cache_version 87 \
    --vcf \
    -i ../data/ALL.wgs.integrated_phase1_release_v3_coding_annotation.20101123.snps_indels.sites.vcf.gz \
    -o codingVep.vcf"

sacct --format="CPUTime,MaxRSS" -j 15307
   CPUTime     MaxRSS 
---------- ---------- 
  10:53:36            
  00:40:51    545480K 
```

NB: Connection to DB fails due to missing perl module `BDB::mysql` perl 

There is a nice summary output in [codingVep.vcf_summary.html](codingVep.vcf_summary.html)


```
1   1577084 .   TTCCTCCTCC  T   505 PASS    THETA=0.0021;AVGPOST=0.6147;AN=2184;VT=INDEL;LDAF=0.3179;ERATE=0.0502;AC=362;RSQ=0.2402;AF=0.17;ASN_AF=0.16;AMR_AF=0.19;AFR_AF=0.09;EUR_AF=0.21;VA=1:CDK11A:ENSG00000008128.14:-:deletionNFS:4/14:CDK11A-201:ENST00000317673.7:2352_951_317_EEEE->E:CDK11A-202:ENST00000340677.5:2319_918_306_EEEE->E:CDK11A-208:ENST00000407249.3:2358_957_319_EEEE->E:CDK11A-204:ENST00000341832.6:2217_816_272_EEEE->E;
CSQ= -|inframe_deletion|MODERATE|CDK11B|ENSG00000248333|Transcript|ENST00000317673|protein_coding|10/21||||942-950|942-950|314-317|EEEE/E|gaGGAGGAGGAa/gaa|||-1||HGNC|1729,
     -|inframe_deletion|MODERATE|CDK11B|ENSG00000248333|Transcript|ENST00000340677|protein_coding|9/20||||909-917|909-917|303-306|EEEE/E|gaGGAGGAGGAa/gaa|||-1||HGNC|1729,
     -|inframe_deletion|MODERATE|CDK11B|ENSG00000248333|Transcript|ENST00000341832|protein_coding|9/20||||807-815|807-815|269-272|EEEE/E|gaGGAGGAGGAa/gaa|||-1||HGNC|1729,
     -|inframe_deletion|MODERATE|CDK11B|ENSG00000248333|Transcript|ENST00000407249|protein_coding|10/21||||948-956|948-956|316-319|EEEE/E|gaGGAGGAGGAa/gaa|||-1||HGNC|1729,
     -|inframe_deletion|MODERATE|CDK11B|ENSG00000248333|Transcript|ENST00000513088|protein_coding|4/15||||445-453|447-455|149-152|EEEE/E|gaGGAGGAGGAa/gaa|||-1|cds_start_NF|HGNC|1729
```

