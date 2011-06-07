show work_mem;
set work_mem= '512MB'
show shared_buffers;
show effective_cache_size;
show fsync;
show random_page_cost;
set random_page_cost = 100;

---- CAGE tags
drop table if exists tritume.cage_tags;
create table tritume.cage_tags AS(
  -- Extract the mapped CAGE tags
  select 
    rname || '_' || CASE WHEN flag = 0 THEN pos + length(seq)
                         WHEN flag = 16 THEN pos END AS tss_id,
    rname,
    (CASE WHEN flag = 0 THEN true 
         WHEN flag = 16 THEN false END)::boolean AS strand,
    CASE WHEN flag = 0 THEN pos + length(seq)
         WHEN flag = 16 THEN pos END AS tss
  from sam_bowtie_human_vs_pig
  where flag in (0, 16)
--  LIMIT 5000000
  ); -- 380984 ms. whole
create index ind_cage_tss_rname on tritume.cage_tags(tss, rname);
create index ind_cage_tss on tritume.cage_tags(tss); -- 43281 ms. for LIMIT 10,000,000; 266734 ms. whole
create index ind_cage_rname on tritume.cage_tags(rname);
vacuum analyze tritume.cage_tags; -- 47125 ms.

select count(*) from tritume.cage_tags limit 100;

----- GTF file with the tss of each exon
drop table if exists tritume.gtf;
create table tritume.gtf AS(
  -- Pull out all the exons from the gtf file
  select distinct 
      gtf_attribute(attributes, 'transcript_id') AS transcript_id,
      rname,
      (CASE WHEN strand = '+' THEN true WHEN strand = '-' THEN false END)::boolean AS strand,
      CASE WHEN strand = '+' THEN f_start WHEN strand = '-' THEN f_end END AS tss -- Start of the exon (start is the "end" if strand is '-')
  FROM sus_scrofa_sscrofa9_56_gtf
  WHERE feature like 'exon'
  ); -- 6.3 sec for 168515 rows
create index ind_gtf_tss_rname on tritume.gtf(tss, rname);
create index ind_gtf_tss on tritume.gtf(tss);
create index ind_gtf_rname on tritume.gtf(rname);
vacuum analyze tritume.gtf;

select * from tritume.gtf limit 100;

/*
METHOD 1: USE BETWEEN JOIN
*/
show all;
set random_page_cost = 4;
set enable_bitmapscan = true;
set enable_indexscan= off;
set enable_hashjoin= ;

set seq_page_cost = 0.1;
set random_page_cost = 0.1;


explain analyze
SELECT tritume.cage_tags.tss_id 
FROM tritume.cage_tags, tritume.gtf
WHERE cage_tags.rname LIKE gtf.rname AND
      cage_tags.strand = gtf.strand AND
      cage_tags.tss BETWEEN gtf.tss - 1000 AND gtf.tss + 1000;

"Nested Loop  (cost=27581.60..22955575714.38 rows=1121381426 width=11) (actual time=22.318..1511560.272 rows=870209 loops=1)"
"  Join Filter: (((cage_tags.rname)::text ~~ (gtf.rname)::text) AND (cage_tags.strand = gtf.strand))"
"  ->  Seq Scan on gtf  (cost=0.00..1809.15 rows=168515 width=7) (actual time=0.051..509.806 rows=168515 loops=1)"
"  ->  Bitmap Heap Scan on cage_tags  (cost=27581.60..69679.21 rows=2661741 width=18) (actual time=3.551..8.845 rows=42 loops=168515)"
"        Recheck Cond: ((cage_tags.tss >= (gtf.tss - 100)) AND (cage_tags.tss <= (gtf.tss + 100)))"
"        ->  Bitmap Index Scan on ind_cage_tss  (cost=0.00..26916.16 rows=2661741 width=0) (actual time=3.363..3.363 rows=42 loops=168515)"
"              Index Cond: ((cage_tags.tss >= (gtf.tss - 100)) AND (cage_tags.tss <= (gtf.tss + 100)))"
"Total runtime: 1513183.031 ms"


"23955671 rows bitmap"
"Nested Loop  (cost=27581.60..22,955,575,714.38 rows=1121381426 width=11) (actual time=4.831..120745.946 rows=870209 loops=1)"
"  Join Filter: (((cage_tags.rname)::text ~~ (gtf.rname)::text) AND (cage_tags.strand = gtf.strand))"
"  ->  Seq Scan on gtf  (cost=0.00..1809.15 rows=168515 width=7) (actual time=4.586..423.088 rows=168515 loops=1)"
"  ->  Bitmap Heap Scan on cage_tags  (cost=27581.60..69679.21 rows=2661741 width=18) (actual time=0.064..0.596 rows=42 loops=168515)"
"        Recheck Cond: ((cage_tags.tss >= (gtf.tss - 100)) AND (cage_tags.tss <= (gtf.tss + 100)))"
"        ->  Bitmap Index Scan on ind_cage_tss  (cost=0.00..26916.16 rows=2661741 width=0) (actual time=0.042..0.042 rows=42 loops=168515)"
"              Index Cond: ((cage_tags.tss >= (gtf.tss - 100)) AND (cage_tags.tss <= (gtf.tss + 100)))"
"Total runtime: 122359.256 ms"

"20 million rows bitmap"
"Nested Loop  (cost=22955.97..19149379592.97 rows=936264164 width=11) (actual time=24.589..92613.293 rows=739814 loops=1)"
"  Join Filter: (((cage_tags.rname)::text ~~ (gtf.rname)::text) AND (cage_tags.strand = gtf.strand))"
"  ->  Seq Scan on gtf  (cost=0.00..1809.15 rows=168515 width=7) (actual time=24.168..378.305 rows=168515 loops=1)"
"  ->  Bitmap Heap Scan on cage_tags  (cost=22955.97..58080.49 rows=2222222 width=18) (actual time=0.059..0.448 rows=35 loops=168515)"
"        Recheck Cond: ((cage_tags.tss >= (gtf.tss - 100)) AND (cage_tags.tss <= (gtf.tss + 100)))"
"        ->  Bitmap Index Scan on ind_cage_tss  (cost=0.00..22400.42 rows=2222222 width=0) (actual time=0.039..0.039 rows=35 loops=168515)"
"              Index Cond: ((cage_tags.tss >= (gtf.tss - 100)) AND (cage_tags.tss <= (gtf.tss + 100)))"
"Total runtime: 93968.328 ms"

"cage_tags = 15 million rows, index scan"
"Nested Loop  (cost=0.01..14,313,681,427.72 rows=702172221 width=11) (actual time=0.396..51861.672 rows=515717 loops=1)"
"  Join Filter: (((cage_tags.rname)::text ~~ (gtf.rname)::text) AND (cage_tags.strand = gtf.strand))"
"  ->  Seq Scan on gtf  (cost=0.00..1809.15 rows=168515 width=7) (actual time=0.025..350.929 rows=168515 loops=1)"
"  ->  Index Scan using ind_cage_tss on cage_tags  (cost=0.01..43273.42 rows=1666667 width=18) (actual time=0.033..0.237 rows=26 loops=168515)"
"        Index Cond: ((cage_tags.tss >= (gtf.tss - 100)) AND (cage_tags.tss <= (gtf.tss + 100)))"
"Total runtime: 52791.029 ms"

"cage_tags = 15,000,000 rows bitmap scan"
"Nested Loop  (cost=18173.37..29,236,971,166.19 rows=702172221 width=11) (actual time=0.410..55211.746 rows=515717 loops=1)"
"  Join Filter: (((cage_tags.rname)::text ~~ (gtf.rname)::text) AND (cage_tags.strand = gtf.strand))"
"  ->  Seq Scan on gtf  (cost=0.00..2925.15 rows=168515 width=7) (actual time=0.023..348.694 rows=168515 loops=1)"
"  ->  Bitmap Heap Scan on cage_tags  (cost=18173.37..131831.04 rows=1666667 width=18) (actual time=0.047..0.251 rows=26 loops=168515)"
"        Recheck Cond: ((cage_tags.tss >= (gtf.tss - 100)) AND (cage_tags.tss <= (gtf.tss + 100)))"
"        ->  Bitmap Index Scan on ind_cage_tss  (cost=0.00..17756.71 rows=1666667 width=0) (actual time=0.032..0.032 rows=26 loops=168515)"
"              Index Cond: ((cage_tags.tss >= (gtf.tss - 100)) AND (cage_tags.tss <= (gtf.tss + 100)))"
"Total runtime: 56145.662 ms"


"Nested Loop  (cost=0.01..8426347686.63 rows=468121433 width=11) (actual time=0.230..24344.143 rows=391738 loops=1)"
"  Join Filter: (((cage_tags.rname)::text ~~ (gtf.rname)::text) AND (cage_tags.strand = gtf.strand))"
"  ->  Seq Scan on gtf  (cost=0.00..2925.15 rows=168515 width=7) (actual time=0.008..341.685 rows=168515 loops=1)"
"  ->  Index Scan using ind_cage_tss on cage_tags  (cost=0.01..22225.75 rows=1111111 width=18) (actual time=0.022..0.095 rows=18 loops=168515)"
"        Index Cond: ((cage_tags.tss >= (gtf.tss - 100)) AND (cage_tags.tss <= (gtf.tss + 100)))"
"Total runtime: 25031.346 ms"

"Nested Loop  (cost=0.01..7958247063.22 rows=936194444 width=11) (actual time=0.222..25215.047 rows=709108 loops=1)"
"  Join Filter: ((cage_tags.rname)::text ~~ (gtf.rname)::text)"
"  ->  Seq Scan on gtf  (cost=0.00..2925.15 rows=168515 width=6) (actual time=0.006..384.630 rows=168515 loops=1)"
"  ->  Index Scan using ind_cage_tss on cage_tags  (cost=0.01..22225.73 rows=1111111 width=17) (actual time=0.023..0.096 rows=18 loops=168515)"
"        Index Cond: ((cage_tags.tss >= (gtf.tss - 100)) AND (cage_tags.tss <= (gtf.tss + 100)))"
"Total runtime: 26446.319 ms"


/*
METHOD 2: EXPAND TABLE gtf
*/
--- Expand gtf table
drop table if exists tritume.gtf_exp;
create table tritume.gtf_exp AS(
    select gtf.*,
    tss + range AS tss_range,
    rname || '_' || (tss + range) AS tss_id
    from tritume.gtf, (select generate_series(-100, +100) AS range) AS t1 -- Change the range in generate_series() to extend the range surrounding each tss
    ); -- 90328 ms.
create index ind_gtf_exp_tss_id on tritume.gtf_exp(tss_id); 
vacuum analyze tritume.gtf_exp;

explain analyze
SELECT tritume.cage_tags.tss_id 
FROM tritume.cage_tags inner join tritume.gtf_exp on
   cage_tags.tss_id = gtf_exp.tss_id;