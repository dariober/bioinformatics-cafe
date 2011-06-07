/*
  Analysis fo pileup file produced by aligning CAGE tags 050810
  against Sscrofa 9 using bwa. see also labbook 20/12/2010
  /exports/work/vet_roslin_nextgen/dario/bwa/output/20101218_cage_050810_increment
*/

-- This stuff will go to schema "cage". Make sure you have it first on the path or transfer 
-- tables as appropriate
show search_path;

-- Import data

-- drop table cage050810
select read_table_alpha($$ file:'F:/data/20101218_CAGE050810/cage_050810_bwa_20101218.clean.single.tss.pileup.gz', dest_table: 'cage050810', header:True, sep:'\t' $$);
alter table cage050810 set schema pileup;
create index ind_cage_pos_rname_strand on cage050810(pos, rname, strand);
create index ind_cage_tss_cluster on cage050810(tss_cluster_id);
vacuum analyze cage050810;
select * from cage050810 limit 100;


-- Max expression of each cluster
create temp table peak_expr AS(
    select tss_cluster_id, max(read_count) AS peak_expr from cage050810 group by tss_cluster_id
    );
create index ind_peak on peak_expr(tss_cluster_id, peak_expr);
vacuum analyze peak_expr;
select * from peak_expr limit 10;

create temp table tss AS(
/* For each trascription cluster, the position of the peak (i.e. most expressed base). 
   When a peak is 'flat', i.e. different bases have the same max(read_count), take the 
   most 5' position. These could be considered as TSS?
*/
    select 
       cage050810.tss_cluster_id, cage050810.rname, cage050810.strand,
       CASE WHEN strand = '+' THEN min(pos) 
            WHEN strand = '-' THEN max(pos) END AS peak_pos,  -- Consider as position of the peak the most 5' base.
       peak_expr AS peak_expr
    from cage050810 inner join peak_expr on
     cage050810.tss_cluster_id = peak_expr.tss_cluster_id and
     cage050810.read_count = peak_expr.peak_expr     
    group by cage050810.tss_cluster_id, cage050810.rname, cage050810.strand, peak_expr
    ); 
-- alter table tss drop column tss_id;
ALTER TABLE tss ADD COLUMN tss_id text;
UPDATE tss SET tss_id = 
    CASE WHEN strand = '+' THEN 'TSSPEAK_'||rname||'_F_'||lpad(peak_pos::text, 9, '0'::text)
         WHEN strand = '-' THEN 'TSSPEAK_'||rname||'_R_'||lpad(peak_pos::text, 9, '0'::text) END;
alter table tss add constraint pk_tss primary key (tss_id);
create unique index ind_tss on tss(tss_cluster_id);
vacuum analyze tss;
select * from tss order by peak_pos limit 100;


-- Stats of each TSS cluster
-- drop table tss_cluster_stats;
create table tss_cluster_stats AS(
  select 
    cage050810.tss_cluster_id,
    min(cage050810.rname)AS rname,
    min(cage050810.strand) AS strand,
    min(cage050810.pos)::int AS cluster_start,
    max(cage050810.pos)::int AS cluster_end,
    max(cage050810.pos)::int-min(cage050810.pos)::int AS cluster_width,
    max(peak_expr)::int as peak_height,
    max(peak_pos)::int as peak_pos,
    sum(read_count)::int AS cluster_expr,
    kurtosis(array_prepend(0::double precision,                         -- |
             array_append(array_agg(read_count),                        -- | Cluster kurtosis. Note a zero has been appended and prepended to the array.
                           0::double precision))) AS cluster_kurtosis   -- |
  from cage050810 inner join tss on cage050810.tss_cluster_id = tss.tss_cluster_id
  group by cage050810.tss_cluster_id
); -- 153698 ms.
comment on table tss_cluster_stats is 'Statistics about cage clusters generated from library CAGE050810; see 20101219_cage_gtf.sql; labbook 20/12/2010'
alter table tss_cluster_stats add constraint pk_tss_cluster primary key (tss_cluster_id);
create index ind_cage_tss on tss_cluster_stats(peak_pos, rname, strand);
create index ind_cage_tc_start on tss_cluster_stats(cluster_start, rname, strand);
create index ind_cage_tc_end on tss_cluster_stats(cluster_end, rname, strand);
vacuum analyze tss_cluster_stats;
select count(*) from tss_cluster_stats limit 100;
-- alter table tss_cluster_stats set schema cage;

/*
 Identify clusters near the 5'ends of ensembl genes:
*/
-- Extract 5' ends from sscrofa gtf:
-- create index ind_transcript on sus_scrofa_sscrofa9_59_gtf (gtf_attribute(attributes, 'transcript_id'));
-- vacuum analyze sus_scrofa_sscrofa9_59_gtf;
create temp table ensembl_tss AS(
  select gtf_attribute(attributes, 'transcript_id') AS transcript_id, min(rname) AS rname, min(strand) AS strand, min(f_start) AS ensembl_tss
  from sus_scrofa_sscrofa9_59_gtf
  where strand = '+'
  group by gtf_attribute(attributes, 'transcript_id')

  union select gtf_attribute(attributes, 'transcript_id') AS transcript_id, min(rname) AS rname, min(strand) AS strand, max(f_end) AS ensembl_tss
  from sus_scrofa_sscrofa9_59_gtf
  where strand = '-'
  group by gtf_attribute(attributes, 'transcript_id')
  );
create index ind_ensembl_tss on ensembl_tss(ensembl_tss, rname, strand);
vacuum analyze ensembl_tss;
select count(*) from ensembl_tss limit 100;

-- Attach tss-clusters to ensembl transcripts:
/*
  Collect all tag clusters that map near the 5' end of ensembl transctipts (table ensembl_tss).
  Transcripts are associated to clusters if their 5'-end (col. ensembl_tss) is no more than 1000 bp from either 
  the start or the end of the cluster.
*/
-- drop table tss_to_transcripts;
create table tss_to_transcripts AS(
  select 
       tss_cluster_stats.*,       
       ensembl_tss.transcript_id,
       ensembl_tss.ensembl_tss,
       cluster_start - ensembl_tss AS enstss2cluster_start, -- Distance of the ensembl tss to the start (left-most pos) of the cluster. When strand= '+', a negative dist. means the cluster is upstream
       cluster_end - ensembl_tss AS enstss2cluster_end      -- Distance of the ensembl tss to the end (right-most pos) of the cluster.
  from tss_cluster_stats, ensembl_tss
  where tss_cluster_stats.rname = ensembl_tss.rname and
      tss_cluster_stats.strand = ensembl_tss.strand and
      -- Set stringency for annotating transcripts here:
      ((tss_cluster_stats.cluster_start between ensembl_tss.ensembl_tss - 1000 and ensembl_tss.ensembl_tss + 1000) or
       (tss_cluster_stats.cluster_end   between ensembl_tss.ensembl_tss - 1000 and ensembl_tss.ensembl_tss + 1000))
  );
comment on table tss_to_transcripts is $$Transcripts from S.scrofa 9.59 (gtf) annotated with transcription clusters from CAGE050810, see script 20101219_cage_gtf.sql and labbook 20/12/2010.
Cage tags have been grouped into clusters and transcripts are associated to clusters if their 5'-end (col. ensembl_tss) is no more than 1000 bp from either the start or the end of the cluster.$$ --'

-- Summary stats for the relation clusters-transcripts:
select 0 AS col_pos, 'No. clusters' AS feature, count(*) AS stats from tss_cluster_stats
union select 1, 'Tot. no. ensembl transcripts', count(distinct gtf_attribute(attributes, 'transcript_id')) from sus_scrofa_sscrofa9_59_gtf
union select 2, 'No. annotated transcripts', count(distinct transcript_id) from tss_to_transcripts
union select 3, 'No. clusters assigned to transcripts', count(distinct tss_cluster_id) from tss_to_transcripts
union select 4, 'Total expression in clusters', sum(cluster_expr) from tss_cluster_stats
union select 5, 'Total expression of cluster assigned to transcr.', sum(cluster_expr) from tss_to_transcripts
union select 6, 'Median expression of cluster peak (height)', (r_median(array_agg(peak_height)))::int from tss_to_transcripts
union select 7, 'Mean distance [cluster (peak) - transcript (5-end)]\n(Positive: cluster is downstream)', avg(CASE WHEN strand = '+' THEN peak_pos - ensembl_tss WHEN strand = '-' THEN ensembl_tss - peak_pos END) from tss_to_transcripts
union select 8, 'St.Dev. distance [cluster (peak) - transcript (5-end)]', stddev(CASE WHEN strand = '+' THEN peak_pos - ensembl_tss WHEN strand = '-' THEN ensembl_tss - peak_pos END) from tss_to_transcripts
order by col_pos


/*
  Relation between CAGE clusters and RNAseq: ENSEMBL based analysis
*/

-- Intersect RNAseq transcript (modelled on ensembl GTF) with tss clusters mapping near ensembl transcripts

-- RPKM of core promoter region (promoter= transcription clusters mapped to the same transcript)
-- drop table promoters;
create table promoters AS(
  select 'SSPROM_' || min(rname) || CASE WHEN min(strand) = '+' THEN '_F_' WHEN min(strand) = '-' THEN '_R_' END  || lpad(min(cluster_start)::text, 9, 0::text) AS promoter_id,
       transcript_id,
       min(rname) AS rname, min(strand) AS strand, 
       min(cluster_start) AS promoter_start, 
       max(cluster_end) AS promoter_end,
  --           |------ C ------|      |---------------------- N ----------------------|   |------- L ------|  See labbook 26/04/2010 for how RPKM has been calculated. 
       (10^9 * sum(cluster_expr)) / ( (select sum(cluster_expr) from tss_cluster_stats) * sum(cluster_width) ) AS promoter_rpkm
  from tss_to_transcripts 
  group by transcript_id
  );
select * from promoters;

drop table rnaseq_ctrl_gtf;
create temp table rnaseq_ctrl_gtf AS (
    select * from cufflinks_transcript_gtf where source like '20100602_CTRL2_gtf'
    );

select distinct rnaseq_ctrl_gtf.transcript_id, 
                promoter_id, 
                fpkm AS transcript_fpkm, 
                promoters.promoter_rpkm 
from rnaseq_ctrl_gtf inner join promoters on 
    rnaseq_ctrl_gtf.transcript_id = promoters.transcript_id
where feature like 'transcript' and fpkm > 5;

"20100602_CTRL2_denovo"
"20100602_CTRL2_gtf"
"20100602_LPS2_denovo"
"20100602_LPS2_gtf"


--------------------------------[ Tritume ]------------------------------------

select * from cage050810 where tss_cluster_id = 'SSTSC_1_R_069858393';

select count(distinct gtf_attribute(attributes, 'transcript_id')) from sus_scrofa_sscrofa9_59_gtf;

show search_path;

alter table tss_cluster_stats rename column peak_hight to peak_height;

select * from cage050810 where tss_cluster_id = 'SSTSC_1_F_003145759';

select * from cage050810 where tss_cluster_id = 'SSTSC_1_F_005798809';

create table cage050810_gtf AS(
    select min(rname)::text AS rname, 
          'transcript_start'::text AS source,
          'trancription_start_cluster'::text AS feature,
          min(pos) AS f_start,
          max(pos) AS f_end,
          min(strand) AS strand,
          ' tss_cluster "'|| tss_cluster_id || '";' AS attributes
    from cage050810
    group by tss_cluster_id
    );
alter table cage050810_gtf set schema pileup;
select * from cage050810_gtf limit 100;
select count(*) from cage050810_gtf;
select * from cage050810 where tss_cluster_id like 'SSTSC_1_F_000006656' order by pos;
"1";"transcript_start";"trancription_start_cluster";6656;6677;"+";"SSTSC_1_F_000006656"

select distinct source from sus_scrofa_sscrofa9_59_gtf;

select lpad('123', 9, '0')

alter table cage050810 set tablespace hdd_free_agent;

select tss_cluster_id, kurtosis(array_prepend(0::double precision, array_append(array_agg(read_count), 0::double precision))) from cage050810 group by tss_cluster_id limit 10;