TOC
===

<!-- MarkdownTOC -->

- [Get test data](#get-test-data)
- [Run tools](#run-tools)
- [Comparison](#comparison)
  - [Comparing output](#comparing-output)
- [Conclusions](#conclusions)

<!-- /MarkdownTOC -->


We want to evaluate some candidate tools for annotatng variants.

Get test data
==============

We use the coding variants from the 1000 genome project. See [README.coding_functional_annotations_20120326](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/functional_annotation/annotated_vcfs/README.coding_functional_annotations_20120326) 
for more information:

```
ssh t1l1.tcrc.gla.ac.uk

cd /scratch/dberaldi/projects/20170303_annotationToolsEvaluation/data

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/functional_annotation/annotated_vcfs/ALL.wgs.integrated_phase1_release_v3_coding_annotation.20101123.snps_indels.sites.vcf.gz
## wc -l 514269

## Not used:
# wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/functional_annotation/annotated_vcfs/ALL.wgs.integrated_phase1_release_v3_noncoding_annotation_20120330.20101123.snps_indels_sv.sites.vcf.gz
```

Run tools
=========

### Annovar

[website](http://annovar.openbioinformatics.org/en/latest/) 

[ref](https://academic.oup.com/nar/article/38/16/e164/1749458/ANNOVAR-functional-annotation-of-genetic-variants)

[run](./annovar.md)

### CAVA

[website](http://www.well.ox.ac.uk/cava).

[ref](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-015-0195-6)

[run](./cava.md)

### VAGReNT

[website](https://github.com/cancerit/VAGrENT)

[ref](https://www.ncbi.nlm.nih.gov/pubmed/26678383)

[run](./vagrent.md)

### VEP

[website](http://www.ensembl.org/info/docs/tools/vep/index.html) 

[ref](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0974-4)

[run](./vep.md)


At the end: Check we didn't drop any variant

```
cd /scratch/dberaldi/projects/20170303_annotationToolsEvaluation

for vcf in `find ./ -name 'coding*.vcf.gz'`
do
    bcftools view -H $vcf | wc -l ## Must return the same count as in input (514239)
    echo $vcf
done
```

Comparison
==========

### Documentation & support community

Here listed more or less in descending order of quality:

1. VEP: Well documented and supported by [Ensembl](http://www.ensembl.org/info/docs/tools/vep/script/index.html). 

2. Annovar: Well documented ([here](http://annovar.openbioinformatics.org/en/latest/)), 
  actively developed. BUT there seem to be only one developer.

3. CAVA: There is a 35 pages pdf in the source directory but no online docs.

4. VAGReNT: The docs are in published article as [protocol](https://www.ncbi.nlm.nih.gov/pubmed/26678383).
  It's enough to make it work and understand the output anyway.

**All**: Actively developed.

### First release and latest version

Tools    |    First  | Current
-------- | --------- | --------
annovar | 2010       | 01 Feb 2016
cava    | 2015       | 21 Dec 2016 
vagrent | 2014       | 31 Jan 2017  
vep     | 2010       | Dec 2016

First release may not be accurate. Note that VEP is moving to a [new version](https://github.com/Ensembl/ensembl-vep) now in beta.

### Popularity

Number of hits by Googling `<toolname> site:seqanswers.com` or `site:biostars.org`:

Tool    | seqanswers | biostars
--------|------------| --------
annovar | 560        | 1660
vep     | 146         | 1030
cava    | 0          | 1 
vagrent | 0          | 0

### Runtime and memory usage

Extracted from slurm with `sacct --format="CPUTime,Elapsed,MaxRSS,JobName" -j <ids>`

Tools    |    Elapsed (h:m:s) | Memory | Notes
-------- | -------------- | -------- | --------
annovar | 00:02:20 | 1,869,252K | 
cava    | 00:17:28 |  593,800K | Run on 16 threads
vep     | 00:45:41 | 1,406,912K | Can be forked
vagrent | 01:32:20 | 41,068K 

This is annotate 514,239 variants.

### License

* [Annovar](http://www.openbioinformatics.org/annovar/annovar_download_form.php)

> ANNOVAR is freely available to personal, academic and non-profit use only. You cannot redistribute ANNOVAR to other users including lab members. No liability for software usage is assumed.

> If you are a commercial user/service provider, you are required to purchase a license to ANNOVAR from BIOBASE. For more information on how to obtain a commercial license please click here.

* [VEP](http://www.ensembl.org/info/about/legal/code_licence.html) Apache 2.0 licence

* [VAGReNT](https://github.com/cancerit/VAGrENT/blob/dev/LICENSE.TXT) GNU v3

* CAVA MIT license (see source dir).

Comparing output
----------------

At this point we have each annotation tools run on the same VCF file and using the same 
reference annotation (GCRh37/hg19 ensembl release 75/gencode v19).

All tools can produce annotated vcf output. This bash script iterates through a list of positions and for each vcf file produced by 
the tools  extract the annotation in some sort of human readable format.

```
cd /scratch/dberaldi/projects/20170303_annotationToolsEvaluation/comparison/

POS="17:7569542 17:7604852 17:15406542 18:19116002 18:21977636 1:156212583 1:63788983 17:41201196 17:47302389 21:33744932"
VCF="../annovar/codingAnno.hg19_multianno.vcf.gz ../cava/codingCava.vcf.gz ../vagrent/codingVagrent.vcf.gz ../vep/codingVep.vcf.gz"

for pos in $POS
do
    for vcf in $VCF
    do
        rec=`bcftools view -H $vcf $pos 2> /dev/null`
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
        ## Each position goes to a new file with outputs concatenated
    done > cmp_out.${pos}.txt
    echo -e "\n---- LOCUS ------------------" >> cmp_out.${pos}.txt
    echo $rec | cut -d' ' -f1-7 >> cmp_out.${pos}.txt
done 
```

Copy files locally

```
cd /Users/db291g/svn_git/bioinformatics-cafe/trunk/annotationToolsEvaluation/out
rsync --remove-source-files c1l1.tcrc.gla.ac.uk:/scratch/dberaldi/projects/20170303_annotationToolsEvaluation/comparison/cmp_out.*.txt ./
```

See directory [out](out/) for example output.

Example output 

```
---- ANNOVAR ----
ANNOVAR_DATE=2016-02-01
Func.wgEncodeGencodeBasicV19=exonic
Gene.wgEncodeGencodeBasicV19=TVP23C
GeneDetail.wgEncodeGencodeBasicV19=.
ExonicFunc.wgEncodeGencodeBasicV19=frameshift_insertion
AAChange.wgEncodeGencodeBasicV19=TVP23C : ENST00000225576.3 : exon6 : c.466_467insAG : p.R156fs
                                 TVP23C : ENST00000519970.1 : exon6 : c.337_338insAG : p.R113fs
ALLELE_END

---- CAVA ------ 
TYPE=Insertion
TRANSCRIPT=ENST00000225576 : ENST00000519970 : ENST00000522212
GENE=TVP23C : TVP23C : TVP23C-CDRT4
GENEID=ENSG00000175106 : ENSG00000175106 : ENSG00000259024
TRINFO=-/61.3kb/6/1.5kb/276 : -/125.7kb/7/1.0kb/185 : -/127.5kb/7/2.8kb/164
LOC=Ex6 : Ex6 : In5/6
CSN=c.466_467insAG : c.337_338insAG : c.462+42556_462+42557insAG
PROTPOS=156 : 113 : .
PROTREF=R : R : .
PROTALT=. : . : .
CLASS=FS : FS : INT
SO=frameshift_variant : frameshift_variant : intron_variant
IMPACT=1 : 1 : 3
ALTANN=. : . : .
ALTCLASS=. : . : .
ALTSO=. : . : .

---- VAGReNT ------ 
VT=Ins
VD=TVP23C|CCDS11170|r.562_563insag|c.466_467insAG|p.R156fs*11|protein_coding:CDS:exon:insertion:frameshift_variant|SO:0000010:SO:0000316:SO:0000147:SO:0000667:SO:0001589
VC=frameshift
VW=TVP23C|CCDS11170|r.562_563insag|c.466_467insAG|p.R156fs*11|protein_coding:CDS:exon:insertion:frameshift_variant|SO:0000010:SO:0000316:SO:0000147:SO:0000667:SO:0001589

---- VEP ------ 
Allele    Consequence                                   IMPACT      SYMBOL          Gene               Feature_type    Feature            BIOTYPE                    EXON    INTRON    HGVSc                                             HGVSp                                   cDNA_position    CDS_position    Protein_position    Amino_acids    Codons       Existing_variation    DISTANCE    STRAND    SYMBOL_SOURCE            HGNC_ID
CT        intron_variant                                MODIFIER    TVP23C-CDRT4    ENSG00000259024    Transcript      ENST00000522212    protein_coding                     5/6       ENST00000522212.2:c.462+42556_462+42557insAG                                                                                                                                                                 -1        HGNC                               
CT        3_prime_UTR_variant&NMD_transcript_variant    MODIFIER    TVP23C-CDRT4    ENSG00000259024    Transcript      ENST00000557349    nonsense_mediated_decay    6/7               ENST00000557349.1:c.*336_*337insAG                                                        561-562                                                                                                            -1        HGNC                               
CT        frameshift_variant                            HIGH        TVP23C          ENSG00000175106    Transcript      ENST00000225576    protein_coding             6/6               ENST00000225576.3:c.466_467insAG                  ENSP00000225576.3:p.Arg156GlnfsTer11    562-563          466-467         156                 R/QX           cgg/cAGgg                                      -1        HGNC                               
CT        frameshift_variant                            HIGH        TVP23C          ENSG00000175106    Transcript      ENST00000519970    protein_coding             6/7               ENST00000519970.1:c.337_338insAG                  ENSP00000428961.1:p.Arg113GlnfsTer11    561-562          337-338         113                 R/QX           cgg/cAGgg                                      -1        HGNC                               
CT        upstream_gene_variant                         MODIFIER    AC005517.3      ENSG00000223878    Transcript      ENST00000452091    processed_pseudogene                                                                                                                                                                                                                          3637        1         Clone_based_vega_gene              
CT        intron_variant&NMD_transcript_variant         MODIFIER    TVP23C-CDRT4    ENSG00000259024    Transcript      ENST00000481756    nonsense_mediated_decay            4/5       ENST00000481756.3:c.*203+42556_*203+42557insAG                                                                                                                                                               -1        HGNC                               
CT        downstream_gene_variant                       MODIFIER    TVP23C          ENSG00000175106    Transcript      ENST00000581273    nonsense_mediated_decay                                                                                                                                                                                                                       3643        -1        HGNC                               
CT        5_prime_UTR_variant                           MODIFIER    CDRT4           ENSG00000239704    Transcript      ENST00000524205    protein_coding             1/3               ENST00000524205.2:c.-179_-178insAG                                                        380-381                                                                                                            -1        HGNC                               
CT        intron_variant&NMD_transcript_variant         MODIFIER    TVP23C-CDRT4    ENSG00000259024    Transcript      ENST00000518506    nonsense_mediated_decay            5/6       ENST00000518506.2:c.*220+42556_*220+42557insAG                                                                                                                                                               -1        HGNC                               

---- LOCUS ------------------
17 15406542 . C CCT 358 PASS
```

```
cd /scratch/dberaldi/projects/20170303_annotationToolsEvaluation/comparison/
ASCIIGenome -fa /scratch/dberaldi/ref/human_g1k_v37.fasta /scratch/dberaldi/ref/Homo_sapiens.GRCh37.75.gtf.gz ../annovar/codingAnno.hg19_multianno.vcf.gz ../cava/codingCava.vcf.gz ../vagrent/codingVagrent.vcf.gz ../vep/codingVep.vcf.gz
```

Conclusions
===========

All four tools seem to be robust enough to be good candidates for the pipeline. 
However, here they haven't been scrutinized in particularly challenging cases (see discussion [here](
http://blog.goldenhelix.com/ajesaitis/the-sate-of-variant-annotation-a-comparison-of-annovar-snpeff-and-vep/))

DB's (subjective) ranking:

1. **VEP** Well documented & popular. If we use ENSEMBL annotation, we can be fairly sure that tool and annotation are developed side by side. 
  VEP describes variants using the standard [sequence ontology nomenclature](http://www.sequenceontology.org/) and [HGVS](http://varnomen.hgvs.org/) notation.
  By default the output is quite verbose since for each variant it reports all the transcripts and features affected. 
  However, I think this is the right thing to do for an automated pipeline since any automatic prioritization may be questionable. 
  If necessary, option `--pick` *is the best method to use if you are interested only in one consequence per variant* (from docs). As a bonus, VEP
  produces an html summary of the annotation results.

1. **CAVA** CAVA implements a *Clinical Sequence Nomenclature* (CNS) to standardize annotation (but take care of this [pitfall](https://xkcd.com/927/)). 
  This is a failry new tool so it's not clear whether is going to be further developed and accepted by the community. I think for a pipeline it is
  preferable to choose a well established tool. 

1. **Annovar** Like VEP, annovar is well documented & popular. Two possible disadvantages: it uses its own nomenclaure and it is developed by virtually 
  just one person.

1. **VAGReNT** Nothing wrong with this tool but like CAVA it is not widely used and it is not clear to me what it offers in addiation to, say, VEP.