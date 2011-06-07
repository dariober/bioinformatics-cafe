/*
  Import CAGE tags 050810 mapped to pig genome and converted to BED.
  See labook 01/10/2011 about CAGE050810
*/

select read_table($$ file:'F:/data/20101218_CAGE050810/cage_050810_bwa_20101218.nordna.single.bed.gz', table: 'cage_bed', temp:True $$); -- 2528998 rows
select * from cage_bed limit 10;


-- drop table cage050810_tss;
create table cage.cage050810_tss AS (
  -- Identify and quantify the TSS by grouping tags mapping to the same position and strand in the BED file
  -- REMINDER: BED format has coordinates 0-based, add 1 to start and end to have it 1-based.
  select
  -- Forward TSS
    min('TSS_' || v1 || '_' || v2+1 || '_0') AS tss_id,
    v1 AS rname,
    (v2+1)::int AS tss_pos,
    min(v6) AS strand,
    count(v4)::int AS tss_count,
    (count(v4)::numeric/(select count(*) from cage_bed))*10^6 AS tss_tpm
  from cage_bed 
  where v6 like '+' 
  group by v1, v2+1
  union select     
  -- Reverse TSS
    min('TSS_' || v1 || '_' || v3+1 || '_16') AS tss_id,
    v1 AS rname,
    (v3+1)::int AS tss_pos,
    min(v6) AS strand,
    count(v4)::int AS tss_count,
    (count(v4)::numeric/(select count(*) from cage_bed))*10^6 AS tss_tpm
  from cage_bed where v6 like '-' group by v1, v3+1
  );
COMMENT ON TABLE cage050810_tss IS 'See 20110110_cage050810_tss.sql. Pig TSS obtained by mapping CAGE050810 to Sscrofa9, SAM converted to BED and then tags mapping to same position and strand counted. This table is in 1-based genome coordinates.'
select * from cage050810_tss limit 10;


-- Annotate TSS with Sus GTF: Which TSS map within (100) bp from the beginning of any exon?
-- drop table cage050810_gtf;
create table cage.cage050810_gtf AS(
  -- Collect all TSSs that map within (100) bases from the beginning or end of any annotated pig exon.
  select cage050810_tss.tss_id,
       cage050810_tss.rname, 
       cage050810_tss.tss_pos,
       cage050810_tss.strand,
       cage050810_tss.tss_tpm,
       attributes,
       feature,
       f_start,
       f_end,
       (CASE WHEN pig_exons.strand = '+' THEN (tss_pos - f_start)
             WHEN pig_exons.strand = '-' THEN (f_end - tss_pos) END)::int AS tss_exon_dist
  from cage050810_tss inner join (select * from sus_scrofa_sscrofa9_59_gtf where feature like 'exon') AS pig_exons ON
    cage050810_tss.rname = pig_exons.rname AND cage050810_tss.strand = pig_exons.strand
  where 
    (pig_exons.strand = '+' and tss_pos between f_start - 1000 and f_start + 1000) OR
    (pig_exons.strand = '-' and tss_pos between f_end - 1000 and f_end + 1000)
  ); -- 1925744 ms (32 min)

select 
    (select count(distinct tss_id) from cage050810_tss) AS no_tss,
    (select count(distinct gtf_attribute(attributes, 'gene_id')) from sus_scrofa_sscrofa9_59_gtf) AS no_ensembl_59_genes,
    (select count(distinct tss_id) from cage050810_gtf) AS no_annotated_tss, 
    (select count(distinct tss_id) from cage050810_gtf where abs(tss_exon_dist) <= 100) AS no_annotated_tss_100bp,
    (select count(distinct gtf_attribute(attributes, 'gene_id')) from cage050810_gtf) AS no_tagged_genes,
    (select count(distinct gtf_attribute(attributes, 'gene_id')) from cage050810_gtf where abs(tss_exon_dist) <= 100) AS no_tagged_genes_100bp;

-- Amount of expression mapped to exons 5'ends.
select sum(tss_tpm) from
    (select distinct tss_id, tss_tpm from cage050810_gtf) AS expr; -- 809129 <- Total expression (tag count) mapped to pig exons
select sum(tss_tpm) from cage050810_tss; -- 1000000 Total expression.

------------------------ Produce Transcription Start Clusters (CTSS) ----------------------

select cluster_reads($$ select rname, strand, tss_pos AS tss_start, tss_pos AS tss_end, tss_id from cage050810_tss $$, 'cage050810_ctss', -50);

-- bind column of trascription clusters to mphage_mcyte_hg19_on_ss9_ctss:

alter table cage050810_tss add column ctss_id text;
update cage050810_tss set ctss_id = tag_cluster from cage050810_ctss
  where cage050810_tss.tss_id = cage050810_ctss.tss_id
comment on column cage050810_tss.ctss_id is 'ID for transcription start clusters produced by grouping TSSs separated by no more than 50 bp. See 20110110_cage050810.sql';
drop table cage050810_ctss;

select * from cage050810_tss limit 10;

select * from public.cage050810_ctss;
select * from cage050810_tss where strand = '+' order by rname, tss_pos limit 100;
select * from cage050810_tss where tss_pos between 3868537 - 100 and 3868537 + 100 and rname like '3' order by tss_pos;
select * from mphage_mcyte_hg19_on_ss9_ctss where ctss_pos between 3868537 - 100 and 3868537 + 100 and rname like '3';

select * from macrophage_monocyte_derived_hg19_ctss_bed where rname like 'chr5' and f_start > 149432854 and f_end < 149492935 order f_start;
