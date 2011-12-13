set work_mem = '512MB';
show work_mem;

drop table if exists cage_tags;
create temp table cage_tags AS(
  -- Extract the mapped CAGE tags
  select 
    rname || '_' || CASE WHEN flag = 0 THEN pos + length(seq)
                         WHEN flag = 16 THEN pos END 
          || '_' || CASE WHEN flag = 0 THEN '+' 
                         WHEN flag = 16 THEN '-' END AS tss_id,
    qname,
    rname,
    pos,
    CASE WHEN flag = 0 THEN pos + length(seq)
         WHEN flag = 16 THEN pos END AS tss,
    (CASE WHEN flag = 0 THEN '+' 
         WHEN flag = 16 THEN '-' END)::varchar(1) AS strand
  from sam_bowtie_human_vs_pig
  where flag in (0, 16)
--  limit 1000000
  ); -- 8030 ms for limit 1000000; 311750 ms entire table.
-- Some or all of these indexes are not necessary. Especially 2) and 3)
-- 1)
create index index_cage_tags_tss_id ON cage_tags (tss_id); -- 635922 ms entire table
-- 2)
create index index_cage_tags ON cage_tags (tss, rname, strand); -- 2844 ms when 'limit 1000000'; For both 'create table' and 'create index': 337320 ms (5.62 min) on entire table
-- 3)
create index index_cage_tags_tss ON cage_tags (tss);

vacuum analyze cage_tags; -- 141954 ms. entire table

--------------------------------[ Prepare GTF ]---------------------------------

drop table if exists gtf;
create temp table gtf AS(
  -- Pull out all the exons from the gtf file
  select distinct 
      gtf_attribute(attributes, 'transcript_id') AS transcript_id,
      gtf_attribute(attributes, 'exon_number') AS exon_number,
      rname,
      strand::varchar(1),
      f_start,
      f_end,
      CASE WHEN strand = '+' THEN f_start WHEN strand = '-' THEN f_end END AS tss -- Start of the exon (start is the "end" if strand is '-')
  FROM sus_scrofa_sscrofa9_56_gtf
  WHERE feature like 'exon'
  ); -- 6.3 sec for 168515 rows
create index index_gtf on gtf (tss, rname, strand);

vacuum analyze gtf;

drop table if exists gtf_exp;
create temp table gtf_exp AS(
explain analyze
    select gtf.*,
    tss + range AS tss_range,
    rname || '_' || (tss + range) || '_' || strand AS tss_id
    from gtf, (select generate_series(-100, +100) AS range) AS t1 -- Change the range in generate_series() to extend the range surrounding each tss
    ); -- 95875 ms (1.5 min)
create index index_gtf_exp_tss_id on gtf_exp (tss_id); -- 895906 ms. when generate_series(-100, +100)
create index index_gtf_exp on gtf_exp (tss_range, rname, strand);
vacuum analyze gtf_exp; -- 341179 ms. (5.7 min) for generate_series(-100, +100)

--------------------------[ Fetch exons with cage tags ]-----------------------

drop table if exists gtf_to_cage;
create temp table gtf_to_cage AS (
  -- This table contains all the CAGE tags mapping near the tss of ensembl exons
  -- 'Near' is defined by generate_series() above
  explain analyze
  select gtf_exp.transcript_id, 
         exon_number, 
         gtf_exp.rname, 
         gtf_exp.strand, 
         f_start, 
         f_end, 
         gtf_exp.tss_id, 
         cage_tags.tss AS cage_tag_tss, 
         cage_tags.qname AS cage_qname
  from gtf_exp inner join cage_tags on gtf_exp.tss_id = cage_tags.tss_id
  ); -- 1641369 ms (27.35 min) for generate_series(-100, +100); 1804722 ms
  
/* EXPLAIN ANALYZE for the above query -------------------------------------------------------------------------------------------------------+

Merge Join  (cost=10071866.28..10921194.60 rows=29363192 width=71) (actual time=1336536.674..1802728.333 rows=870209 loops=1)
  Merge Cond: (gtf_exp.tss_id = cage_tags.tss_id)
  ->  Sort  (cost=5997651.21..6082330.00 rows=33871516 width=46) (actual time=764598.798..923303.778 rows=33871515 loops=1)
        Sort Key: gtf_exp.tss_id
        Sort Method:  external merge  Disk: 1968048kB
        ->  Seq Scan on gtf_exp  (cost=0.00..719467.16 rows=33871516 width=46) (actual time=37.244..83653.887 rows=33871515 loops=1)
  ->  Materialize  (cost=4074198.14..4373644.04 rows=23955672 width=38) (actual time=571914.491..746072.422 rows=24021816 loops=1)
        ->  Sort  (cost=4074198.14..4134087.32 rows=23955672 width=38) (actual time=571914.481..662178.282 rows=23954997 loops=1)
              Sort Key: cage_tags.tss_id
              Sort Method:  external merge  Disk: 1155376kB
              ->  Seq Scan on cage_tags  (cost=0.00..482928.72 rows=23955672 width=38) (actual time=21.940..78089.354 rows=23955671 loops=1)
Total runtime: 1804722.576 ms
*/

create table ensembl_exon_caged AS ( select * from gtf_to_cage );
comment on table ensembl_exon_caged is 'See script 20100828_map_cage2exons-2.sql. Human CAGE tags mapping in +/- 100 bp from the TSS of each emsembl exon';
select count(*) from gtf_to_cage limit 100;


--------------------------------[ Tritume ]-------------------------------------

create index ind_gtf_tss_rname_strand on gtf(tss, rname, strand);
create index ind_gtf_tss on gtf(tss);
create index ind_gtf_rname on gtf(rname);
create index ind_gtf_strand on gtf(strand);
vacuum analyze gtf;

explain analyze
select gtf.transcript_id, 
       exon_number, 
       gtf.rname, 
       gtf.strand, 
       f_start, 
       f_end,
       gtf.tss, 
       cage_tags.tss AS cage_tag_tss, 
       cage_tags.qname AS cage_qname
from gtf, cage_tags
where gtf.rname = cage_tags.rname and
      gtf.strand = cage_tags.strand and 
      cage_tags.tss between gtf.tss - 1000 and gtf.tss + 1000;

select * from gtf limit 100;

select * from sam_bowtie_human_vs_pig limit 100;


"Merge Join  (cost=2867151.68..3200195.16 rows=15067569 width=58)" -- No indexes


select * from gtf_to_cage order by transcript_id, exon_number, tss limit 100;
    gtf_exp.tss_range = cage_tags.tss AND
    cage_tags.rname = gtf_exp.rname AND
    cage_tags.strand = gtf_exp.strand


  57289765 - 57290574


    "Merge Join  (cost=3722351.01..3 428 086 438.28 rows=13 696 737 679 width=41)"

