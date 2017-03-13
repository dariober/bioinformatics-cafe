TOC
===

<!-- MarkdownTOC -->

- [Make synthetic variants](#make-synthetic-variants)
- [Run tools](#run-tools)
- [Collect results](#collect-results)

<!-- /MarkdownTOC -->


Make synthetic variants
=======================


```
cd /scratch/dberaldi/projects/20170303_annotationToolsEvaluation/data

bcftools view -h ALL.wgs.integrated_phase1_release_v3_coding_annotation.20101123.snps_indels.sites.vcf.gz > synth.vcf
echo "1   39330342  .  CTA  CGCAA  100   PASS THETA=.;AN=.;AC=.;VT=.;AA=.;RSQ=.;LDAF=.;SNPSOURCE=.;AVGPOST=.;ERATE=.;AF=.;ASN_AF=.;AMR_AF=.;AFR_AF=.;EUR_AF=.;VA=." \
| sed 's/ \+/\t/g' >> synth.vcf

## Without anchor
echo "1   39339358  .  CA  GCAA  100   PASS THETA=.;AN=.;AC=.;VT=.;AA=.;RSQ=.;LDAF=.;SNPSOURCE=.;AVGPOST=.;ERATE=.;AF=.;ASN_AF=.;AMR_AF=.;AFR_AF=.;EUR_AF=.;VA=." \
| sed 's/ \+/\t/g' >> synth.vcf

bgzip -f synth.vcf && tabix -f synth.vcf.gz && rm synth.vcf
```

Run tools
=========

```
## ANNOVAR
cd /scratch/dberaldi/projects/20170303_annotationToolsEvaluation/annovar
table_annovar.pl ../data/synth.vcf.gz \
    humandb/ \
    -buildver hg19 \
    -out synthAnno \
    -remove \
    -protocol wgEncodeGencodeBasicV19 \
    -operation g \
    -nastring . \
    -vcfinput

## VEP
cd /scratch/dberaldi/projects/20170303_annotationToolsEvaluation/vep
rel=75
variant_effect_predictor.pl \
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
    -i ../data/synth.vcf.gz \
    -o synthVep.vcf

## CAVA
cd /scratch/dberaldi/projects/20170303_annotationToolsEvaluation/cava
~/applications/CAVA-1.2.0/cava.py -t 16 -c config.txt \
    -i ../data/synth.vcf.gz \
    -o synthCava

## VAGReNT
cd /scratch/dberaldi/projects/20170303_annotationToolsEvaluation/vagrent
AnnotateVcf.pl -i ../data/synth.vcf.gz \
    -o synthVagrent.vcf \
    -cache cache/Homo_sapiens.GRCh37.${rel}.vagrent.cache.gz \
    -sp homo_sapiens \
    -as GRCh37

for vcf in `find ../ -name 'synth*.vcf'`
do
    bgzip -f $vcf && tabix -f ${vcf}.gz && rm $vcf
done
```

Collect results
===============

```
cd /scratch/dberaldi/projects/20170303_annotationToolsEvaluation/comparison/

POS="1:39330342 1:39339358"
VCF="../annovar/synthAnno.hg19_multianno.vcf.gz ../cava/synthCava.vcf.gz ../vagrent/synthVagrent.vcf.gz ../vep/synthVep.vcf.gz"

for pos in $POS
do
    for vcf in $VCF
    do
        out=cmp_synth.${pos}.txt
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
    done > $out
    echo -e "\n---- LOCUS ------------------" >> $out
    echo $rec | cut -d' ' -f1-7 >> $out
done 
```

Copy files locally

```
cd /Users/db291g/svn_git/bioinformatics-cafe/trunk/annotationToolsEvaluation/out
rsync --remove-source-files c1l1.tcrc.gla.ac.uk:/scratch/dberaldi/projects/20170303_annotationToolsEvaluation/comparison/cmp_synth.*.txt ./
```