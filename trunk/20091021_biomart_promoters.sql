drop table if exists mart_promoters;
create temp table mart_promoters AS (
  select *, 
    CASE WHEN strand like '1' THEN transcript_start - 300
         WHEN strand like '-1' THEN transcript_end - 50 END AS promoter_start,
    CASE WHEN strand like '1' THEN transcript_start + 50 
         WHEN strand like '-1' THEN transcript_end - 300 END AS promoter_end,
    transcript_end - transcript_start AS transcript_length
  from mart_tome
  );

-- No. of genes and transcripts in pig and human transcriptomes
select 
  genome,
  count(ensembl_transcript_id) AS no_transcripts, 
  count(distinct ensembl_gene_id) AS no_genes   
from mart_tome  
group by genome;

-- Tags common to both H.sapiens and S.scrofa genomes
drop table if exists c_tss ;
create temp table c_tss AS (
  select distinct tmp_hg95_ss9.read_name,
    tmp_hg95_ncbi37.strand AS hs_strand,
    tmp_hg95_ncbi37.refseq_name AS hs_chr,
    tmp_hg95_ncbi37.position AS hs_position,
    tmp_hg95_ss9.strand AS pig_strand,
    tmp_hg95_ss9.refseq_name AS pig_chr, 
    tmp_hg95_ss9.position AS pig_position
  from tmp_hg95_ss9 inner join tmp_hg95_ncbi37
    on tmp_hg95_ss9.read_name = tmp_hg95_ncbi37.read_name
  );
-- drop index Ind_c_tss_hs_chr ;
create index Ind_c_tss_hs_chr on c_tss (hs_chr, hs_position);
create index Ind_c_tss_pig_chr on c_tss (pig_chr, pig_position);


-- drop table if exists mart_prom_hs_1;
-- create temp table mart_prom_hs_1 AS (
--  select * from mart_promoters where genome = 'h37' -- and strand = '1'
--  ); 


drop table if exists tags_near_promoters ;
create temp table tags_near_promoters AS(
  -- Tags in both pig and H.sap. near 5' end of PIG transcripts
  SELECT *, 'ss9' AS ref_tome
  FROM  c_tss INNER JOIN mart_promoters ON pig_chr = chr
  WHERE 
      pig_position > promoter_start AND 
      pig_position < promoter_end AND
      mart_promoters.genome = 'ss9'
      
  -- Tags in both pig and H.sap. near 5' end of HUMAN transcripts
  UNION SELECT *, 'h37' AS ref_tome
  FROM  c_tss INNER JOIN mart_promoters ON hs_chr = chr
  WHERE 
      hs_position > promoter_start AND 
      hs_position < promoter_end AND
      mart_promoters.genome = 'h37'
  );

-- No. of tags mapped against pig and human transcriptome
select 
  ref_tome AS ref_transcriptome,
  count(distinct read_name) AS no_unique_tags,
  count(distinct read_name) / (select count(*) from c_tss)::numeric AS proportion_mapped,
  count(read_name) AS tot_hits  
from tags_near_promoters
group by ref_tome;

-- Agreement between transcript and tag sense/antisense
select ref_tome, 

------------------------------[ Tritume ]--------------------------------------



select count(distinct read_name) from hs_promoters;

select count(*), strand from hs_promoters
group by strand;

select * from mart_promoters limit 10;



limit 10000;



select * from mart_promoters limit 10;
select * from c_tss limit 10;


select * from c_tss limit 10;
select refseq_name, count(*) 
from c_tss
group by refseq_name;

select count(*) from c_tss 
where hs_strand = '+';
limit 10;

select chr, genome, count(*) 
from mart_promoters
group by chr, genome
order by genome, count(*) desc;

