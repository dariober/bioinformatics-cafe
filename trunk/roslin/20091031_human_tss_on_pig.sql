-- Which human CAGE reads (library hg95ctrls, single mappers and located near 
-- 5' ends of refseq human genes) map on the pig genome. 
-- Which of these reads is near known pig transcripts from ensembl.

------------------------------[ pig_ensembl_tr ]-------------------------------
-- Pig transcripts from ensembl

-- drop table if exists pig_ensembl_tr;
create temp table pig_ensembl_tr AS(
  select 
    ensembl_gene_id,
    ensembl_transcript_id,
    chr,
    strand,
    transcript_start,
    transcript_end
    from mart_tome
  where genome like 'ss9'
  );

---------------------------------[ hs_ss_tss ]---------------------------------
-- tss (reads) near human refseq genes mapping on the pig genome

create temp table hs_ss_tss AS (
select 
  out_hg95ctrls_promoters.read_name,
  out_hg95ctrls_promoters.chr AS human_chr,
  out_hg95ctrls_promoters.strand AS human_strand,
  out_hg95ctrls_promoters.tss AS human_tss,
  out_hg95ctrls_promoters.tag AS human_tag,
  out_hg95ctrls_promoters.seq_name AS human_refseq,
  bowtie.strand AS pig_tss_strand,
  bowtie.refseq_name AS pig_chr,
  bowtie.position AS pig_position
from out_hg95ctrls_promoters inner join bowtie on
  out_hg95ctrls_promoters.read_name = bowtie.read_name
where bowtie.job_id = 1
  );

----------------------------------[ Tritume ]----------------------------------

select * from pig_ensembl_tr limit 10;

select length(ensembl_gene_id) from mart_tome group by length(ensembl_gene_id);

select * from hs_ss_tss limit 10;
select count(*) from hs_ss_tss;