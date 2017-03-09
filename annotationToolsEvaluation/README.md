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
- [Comparison](#comparison)

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

* Run and memory footprint

* Use standard sequence ontology?

* Documentation & support community

  * VEP: Well documented and supported by ENSEMBL.

  * Annovar: Well documented ([here](http://annovar.openbioinformatics.org/en/latest/)), 
     actively developed. BUT there seem to be only one developer.

  * CAVA: There is a 35 pages pdf in the source directory but no online docs. 

  * VAGReNT: The docs are in published article as [protocol](https://www.ncbi.nlm.nih.gov/pubmed/26678383).
    It's enough to make it work and understand the output anyway.

More or less in descending order of quality


* Runtime and memory usage

Tools    |    Elapsed (h:m:s) | Memory
-------- | -------------- | --------
annovar |          |
cava    |          |  
vagrent |          |  
vep     |          | 


* First release and latest version

Tools    |    First      | Current
-------- | -------------- | --------
annovar | 15 Feb 2010   | 01 Feb 2016
cava    | 2015         | 21 Dec 2016 
vagrent | 17 Jun 2015   | 31 Jan 2017  
vep     |     2010     | Dec 2016

* Popularity

Number of hits by Googling `<toolname> site:seqanswers.com` or `site:biostars.org`:

Tool    | seqanswers | biostars
--------|------------| --------
annovar | 560        | 1660
cava    | 0          | 1 
vagrent | 0          | 0
vep     | 146         | 1030

* Developing team

** [Annovar](http://annovar.openbioinformatics.org/en/latest/misc/credit/): 1 person 

** [CAVA](http://www.well.ox.ac.uk/cava): 4 authors listed (?)

** [Vagrent](https://github.com/cancerit/VAGrENT/): 2 

** [VEP]()

All: Actively developed.

* License

[Annovar](http://www.openbioinformatics.org/annovar/annovar_download_form.php)

> ANNOVAR is freely available to personal, academic and non-profit use only. You cannot redistribute ANNOVAR to other users including lab members. No liability for software usage is assumed.

> If you are a commercial user/service provider, you are required to purchase a license to ANNOVAR from BIOBASE. For more information on how to obtain a commercial license please click here.


Comparison
==========

At this point we have each annotation tools run on the same VCF files and using the same 
reference annotation (GCRh37/hg19 ensembl release 75/gencode v19).

We need to select some sites of interest

* Select a bunch of SNVs and indels from various sites, like exons, introns, UTRs
  Use the output from Annovar for this?

Let's have a look at the various categories

There is a vcf file which is pretty handy since in thos way we can deal with a standard format. 

```
cd /scratch/dberaldi/projects/20170303_annotationToolsEvaluation/comparison/

bcftools view -H ../annovar/codingAnno.hg19_multianno.vcf.gz | sed 's/.*;Func.wgEncodeGencodeBasicV19=//' | sed 's/;.*//' | sort | uniq -c | sort -k1,1n
      8 upstream\x3bdownstream
     11 UTR5\x3bUTR3
     24 ncRNA_splicing
     32 ncRNA_exonic\x3bsplicing
    165 upstream
    168 downstream
    231 exonic\x3bsplicing
    521 UTR5
    627 ncRNA_intronic
   1251 UTR3
   1617 intergenic
   2551 ncRNA_exonic
   2624 splicing
   3035 intronic
 501374 exonic

bcftools view -H ../annovar/codingAnno.hg19_multianno.vcf.gz | sed 's/.*;ExonicFunc.wgEncodeGencodeBasicV19=//' | sed 's/;.*//' | sort | uniq -c | sort -k1,1n
    173 nonframeshift_insertion
    324 frameshift_insertion
    364 stoploss
    493 frameshift_deletion
    532 nonframeshift_deletion
    600 unknown
   5674 stopgain
  12634 .
 201377 synonymous_SNV
 292068 nonsynonymous_SNV
```

```
cd /scratch/dberaldi/projects/20170303_annotationToolsEvaluation/comparison/

for reg in 17:7569542 17:7604852 17:15406542 18:19116002 18:21977636 1:156212583 1:63788983
do
    for vcf in ../annovar/codingAnno.hg19_multianno.vcf.gz ../cava/codingCava.vcf.gz ../vagrent/codingVagrent.vcf.gz ../vep/codingVep.vcf.gz
    do
        rec=`bcftools view -H $vcf $reg 2> /dev/null`
        if [[ $vcf == *"Anno"* ]]; then
            echo -e "\n---- ANNOVAR ----"
            echo $rec | sed 's/.*;ANNOVAR_DATE=/ANNOVAR_DATE=/' | sed 's/;/\n/g' | sed 's/,/\n                                 /g' | sed 's/:/ : /g'
        elif [[ $vcf == *"Cava"* ]]; then 
            echo -e "\n---- CAVA ------ "
            echo $rec | sed 's/.*;TYPE=/TYPE=/' | sed 's/;/\n/g' | sed 's/:/ : /g'
        elif [[ $vcf == *"Vagrent"* ]]; then
            echo -e "\n---- VAGReNT ------ "
            echo $rec | sed 's/.*;VT=/VT=/' | sed 's/;/\n/g'  
        elif [[ $vcf == *"Vep"* ]]; then
            echo -e "\n---- VEP ------ "
            hdr="Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|SYMBOL_SOURCE|HGNC_ID"
            dat=`echo $rec | sed 's/.*;CSQ=//'`
            echo $hdr","$dat | sed 's/,/\n/g' | sed 's/|/ | /g' | column -s '|' -t
        else
            echo "Invalid file name: $vcf"
            break
        fi
    done > cmp_out.${reg}.txt
    echo -e "\n---- LOCUS ------------------" >> cmp_out.${reg}.txt
    echo $rec | cut -d' ' -f1-7 >> cmp_out.${reg}.txt
done 

rsync --remove-source-files c1l1.tcrc.gla.ac.uk:/scratch/dberaldi/projects/20170303_annotationToolsEvaluation/comparison/cmp_out.*.txt ./
```


```
../annovar/codingAnno.hg19_multianno.vcf.gz
../annovar/noncodingAnno.hg19_multianno.vcf.gz
../cava/codingCava.vcf.gz
../cava/noncodingCava.vcf.gz
../vagrent/codingVagrent.vcf.gz
../vagrent/noncodingVagrent.vcf.gz
../vep/codingVep.vcf.gz
../vep/noncodingVep.vcf.gz

../annovar/codingAnno.hg19_multianno.vcf.gz ../cava/codingCava.vcf.gz ../vagrent/codingVagrent.vcf.gz ../vep/codingVep.vcf.gz


ASCIIGenome -fa /scratch/dberaldi/ref/human_g1k_v37.fasta /scratch/dberaldi/ref/Homo_sapiens.GRCh37.75.gtf.gz ../annovar/codingAnno.hg19_multianno.vcf.gz ../cava/codingCava.vcf.gz ../vagrent/codingVagrent.vcf.gz ../vep/codingVep.vcf.gz
```