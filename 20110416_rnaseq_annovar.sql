/*
  Import output of ANNOVAR
  see /exports/work/vet_roslin_nextgen/dario/annovar/output/20110416_rnaseq_am_bmdm
  and 20110416_annotate_snp.sh
  Annotation relates to SNP calling 20110414_rnaseq_mpileup done on concatenated CTRL and LPS libs
*/

-------------------------------------------------------------------------------
-- Import *.variant_function
-------------------------------------------------------------------------------

select read_table($$ file:'D:/Tritume/var.flt.bmdm_cat.annvar.variant_function', table: 'vcf_variants.variant_function', allow_rugged:True, overwrite:False,
    header: ['file', 'region', 'gene_annot', 'chrom', 'start', 'end', 'allele_1', 'allele_2', 'genotype', 'snp_qual', 'depth', 'rms_mq'],
    apply:
      """def apply(line):
          line= line.replace('\t\t', '\t') ## This is to remove occasional double \t after the second allele base(s) (bug in annovar?)
          line= 'var.flt.bmdm_cat.vcf' + '\t' + line
          return(line)
"""$$);

select read_table($$ file:'D:/Tritume/var.flt.am_cat.annvar.variant_function', table: 'vcf_variants.variant_function', allow_rugged:True, append:True,
    header: ['file', 'region', 'gene_annot', 'chrom', 'start', 'end', 'allele_1', 'allele_2', 'genotype', 'snp_qual', 'depth', 'rms_mq'],
    apply:
      """def apply(line):
          line= line.replace('\t\t', '\t') ## This is to remove occasional double \t after the second allele base(s) (bug in annovar?)
          line= 'var.flt.am_cat.vcf' + '\t' + line
          return(line)
"""$$);
comment on table variant_function is 'Annotation for SNPs in table vcf see labbook 16/04/2011. Output of ANNOVAR *.variant_function: "annotate_variation.pl -geneanno -buildver susScr2 -dbtype ensgene var.flt.bmdm_cat.annvar $annodb" see 20110416_annotate_snp.sh.'

-------------------------------------------------------------------------------
-- Import *.exonic_variant_function
-------------------------------------------------------------------------------

-- drop table if exists exonic_variant_function;
select read_table($$ file: 'D:/Tritume/var.flt.am_cat.annvar.exonic_variant_function', table:'vcf_variants.exonic_variant_function',
apply:
"""def apply(line):
    line= 'var.flt.am_cat.vcf' + '\t' + line  ## This should be the original .vcf file and should match column variant_function.file
    return(line)
""", header: ['file', 'line_no', 'consequence', 'annotation', 'chrom', 'start', 'end', 'allele_1', 'allele_2', 'genotype', 'snp_qual', 'depth', 'rms_mq']  $$)

select read_table($$ file: 'D:/Tritume/var.flt.bmdm_cat.annvar.exonic_variant_function', table:'vcf_variants.exonic_variant_function', append:True, allow_rugged:True,
apply:
"""def apply(line):
    line= 'var.flt.bmdm_cat.vcf' + '\t' + line  ## This should be the original .vcf file and should match column variant_function.file
    return(line)
""", header: ['file', 'line_no', 'consequence', 'annotation', 'chrom', 'start', 'end', 'allele_1', 'allele_2', 'genotype', 'snp_qual', 'depth', 'rms_mq']  $$)

-- 

select * from exonic_variant_function where rms_mq is null;

-------------------------------------------------------------------------------
-- Summary stats
-------------------------------------------------------------------------------

-- No. variants by file and consequence
select file, consequence, count(*) as no_snp
from exonic_variant_function 
group by file, consequence
order by file, consequence;

-- Same as above, by genotype also:
select file, consequence, genotype, count(*) as no_snp
from exonic_variant_function 
group by file, consequence, genotype
order by file, consequence, genotype;

select * from cufflinks_transcript_gtf
where transcript_id in 
    (select distinct substring(annotation from 20 for 18 ) from exonic_variant_function where consequence like 'stopgain SNV' and file = 'var.flt.am_cat.vcf')
and feature = 'transcript' and source like '20110202_am_%'
order by gene_id, source;

------------------------------------------------------------------------------
-- Genes carrying variants
-------------------------------------------------------------------------------

drop table vcf_genes;
create temp table vcf_genes AS(
  select distinct ensembl_transcript_id, external_transcript_id, exonic_variant_function.*
  from exonic_variant_function_genes, exonic_variant_function
  where exonic_variant_function.annotation like '%' || exonic_variant_function_genes.ensembl_transcript_id || '%'
);

select distinct * 
from vcf_genes
where consequence = 'stopgain SNV' 

select distinct * 
from vcf_genes, cufflinks_transcript_gtf
where consequence = 'stopgain SNV' and
      vcf_genes.ensembl_transcript_id = cufflinks_transcript_gtf.transcript_id and
  feature = 'transcript' and 
  file = 'var.flt.am_cat.vcf' and source in('20110202_am_ctrl', '20110202_am_lps') and
  fpkm > 10


UNION select distinct * 
from vcf_genes, cufflinks_transcript_gtf
where 
    consequence = 'stopgain SNV' and
    vcf_genes.ensembl_transcript_id = cufflinks_transcript_gtf.transcript_id and
    feature = 'transcript' and 
file = 'var.flt.bmdm_cat.vcf' and source in('20110207_bmdm_ctrl', '20110207_bmdm_lps')
order by transcript_id, file, source;

select external_transcript_id, 
       count(external_transcript_id) as count_all, 
       count(distinct external_transcript_id) AS count_distinct 
from vcf_genes 
group by external_transcript_id
order by count_all desc;


------------------------------------------------------------------------------
-- TRITUME
-------------------------------------------------------------------------------

select * from cufflinks_transcript_gtf where transcript_id like 'ENSSSCT00000001533' and feature like 'transcript'
union select * from cufflinks_transcript_gtf_old where transcript_id like 'ENSSSCT00000001533' and feature like 'transcript'
order by source;

select * from deseq_nbinom where id like 'ENSSSCG00000001404';


select distinct chrom, start, "end" from exonic_variant_function where file ='var.flt.am_cat.vcf';

select distinct(region), count(*)
from exonic_variant_function inner join variant_function on 
    exonic_variant_function.file = variant_function.file and
    exonic_variant_function.chrom = variant_function.chrom and
    exonic_variant_function.start = variant_function.start and
    exonic_variant_function.end = variant_function.end
where exonic_variant_function.file = 'var.flt.am_cat.vcf'
group by region;    

select file, region, count(distinct chrom||start) AS no_snp from variant_function group by file, region order by file, no_snp desc;

select count(*), file, genotype from variant_function group by file, genotype;

select * from variant_function where region like 'exonic;splicing';
select * from variant_function order by snp_qual desc limit 10;
select * from variant_function order by rms_mq desc limit 10;
select * from vcf where pos = 37797689;
select distinct genotype from variant_function limit 10;

select * from cufflinks_transcript_gtf where gene_id like 'ENSSSCG00000007585' and feature like 'transcript';
select * from cufflinks_transcript_gtf_old where gene_id like 'ENSSSCG00000007585' and feature like 'transcript';
