-- Maps cage reads from RIKEN library hg95ctrls against the vicinity of 5' end 
-- of human RefSeq genes (hg19 build from Genome Browser).
-- CAGE reads filtered to have only single mappers (163750)


-------------------[ Promoter boundaries from refseq genes ]-------------------

-- drop table if exists refseq_promoters ;
create temp table refseq_promoters AS(
  select 
    seq_name, 
    ltrim(chrom, 'chr') AS chr,
    CASE WHEN strand like '+' THEN true WHEN strand like '-' THEN false END AS strand,
    chromstart,
    chromend,
    -- Promoter window size
    CASE WHEN strand like '+' THEN chromstart - 1000 
         WHEN strand like '-' THEN chromend + 1000 END AS promoter_start,
    CASE WHEN strand like '+' THEN chromstart + 1000 
         WHEN strand like '-' THEN chromend - 1000 END AS promoter_end
  from refseq_hg19
  );


------------[ bwt_single: Human CAGE tags hg95ctrls: single mappers ]----------

-- drop table if exists single_mappers ;
create temp table single_mappers AS ( -- List of tags with one hit only
  select read_name
  from bowtie 
  where job_id = 2
  group by read_name 
  having count(*) = 1
  );

-- drop table if exists bwt_single ;
create temp table bwt_single AS (   -- Annotate tags from previous query 
                                    -- with genome position and strand info
  select 
    bowtie.read_name, 
    bowtie.refseq_name AS chr,
    bowtie.strand, 
    CASE WHEN strand = true THEN bowtie."position" 
         WHEN strand = false THEN bowtie."position" + length(read_sequence) 
           END AS tss
  from bowtie inner join single_mappers on single_mappers.read_name = bowtie.read_name
  where job_id = 2
  );

-----------------------[ Single mappers near promoters ]-----------------------

-- drop table if exists out_hg95ctrls_promoters ;
create table out_hg95ctrls_promoters AS (
  select 
    bwt_single.read_name, 
    bwt_single.chr, 
    bwt_single.strand,
    bwt_single.chr || '_' || tss || '_' || bwt_single.strand AS tag,
    seq_name,
    refseq_promoters.strand AS refseq_strand,
    chromstart,
    chromend,
    promoter_start,
    promoter_end,
    -- Distance between tag (TSS) and transcript 5' end. 
    -- Negative values: The TSS is upstream to the 5' end.
    CASE WHEN refseq_promoters.strand = true THEN tss - chromstart
         WHEN refseq_promoters.strand = false THEN chromend - tss END AS tss_distance
  from bwt_single inner join refseq_promoters 
    on bwt_single.chr = refseq_promoters.chr
  where 
    (bwt_single.tss < promoter_end AND -- For sense transcripts
    bwt_single.tss > promoter_start)
    OR
    (bwt_single.tss > promoter_end AND -- For anti-sense transcripts
    bwt_single.tss < promoter_start)
  );

comment on table out_hg95ctrls_promoters is 
'Tags from RIKEN library hg95ctrls mapping near the 5'' end of human refseq genes (UCSC). 
Single mappers only. 
Table generated with 20091028_hg95ctrls_refseq_promoters.sql';

-------------[ Group reads mapping to the same position into tags ]------------

-- drop view if exists qry_hg95ctrls_promoters ;
create view qry_hg95ctrls_promoters AS (
  select distinct 
    tag, 
    chr,
    tss,
    strand,
    refseq_strand,
    seq_name,
    chromstart,
    chromend,
    promoter_start,
    promoter_end,
    tss_distance
  from out_hg95ctrls_promoters
  );
  
comment on view qry_hg95ctrls_promoters is 
'CAGE tags from library hg95ctrls mapping near 5'' end of human refseq genes. 
Single mappers only. 
Reads mapping to the same position grouped into tags.
Query executed with 20091028_hg95ctrls_refseq_promoters.sql';

-----------------------------[ Summary stats ]---------------------------------

-- Reads in the vicinity of a known 5' end
select 
  (select count(*) from single_mappers) AS total_single_mappers,
  count(*) AS no_hits,
  count(distinct read_name) AS mapped_reads,
  (select count(distinct read_name) from out_hg95ctrls_promoters) / (select count(*) from single_mappers)::numeric AS percentage_mapped,
  count(distinct seq_name) AS no_unique_genes_hit,
  count(distinct seq_name)::numeric / (select count(distinct seq_name) from refseq_hg19) AS percentage_genes_hit
from out_hg95ctrls_promoters;

-- tags in the vicinity of a known 5' end
select 
  count(*) AS no_mapped,
  (select count(*) from (select distinct chr, strand, tss from bwt_single) as t1)
    AS total_single_tags, 
  count(*)::numeric / (select count(*) from (select distinct chr, strand, tss from bwt_single) as t1)
    AS percentage_mapped
from qry_hg95ctrls_promoters;

-- Disagreement bwt transcript and tag orientation. I.e. tags mapped in one 
-- orientation map near a transcript with opposite orientation
select count(*), count(*)::numeric /  (select count(*) from single_mappers) AS percentage_opposite
from out_hg95ctrls_promoters
where strand != refseq_strand;

-- Summary table of cage and gene orientation
select strand AS cage_strand, refseq_strand, 
  CASE WHEN tss_distance <= 0 THEN 'upstream' 
       WHEN tss_distance > 0 THEN 'downstream' END AS tss_rel_pos,
  count(*), 
  100*count(*)::numeric / (select count(*) from out_hg95ctrls_promoters) AS percentage,
  avg(tss_distance) AS avg_tss_distance,
  stddev(tss_distance) AS stdev_tss_distance
from out_hg95ctrls_promoters
group by strand, refseq_strand, tss_rel_pos
order by percentage, strand, refseq_strand, tss_rel_pos;

-- RefSeq genes with associated a tag (tagged refseq genes)
select 
  (select count(distinct seq_name) from refseq_hg19) AS no_refseq_genes,
  (select count(distinct seq_name) from out_hg95ctrls_promoters ) 
    AS no_tagged_refseq_genes,
  (select count(distinct seq_name) from out_hg95ctrls_promoters ) / 
    (select count(distinct seq_name) from refseq_hg19)::numeric 
    AS percentage_tagged
  ;

-------------------------------[ Tritume ]-------------------------------------

select * from out_hg95ctrls_promoters limit 100;

select * from bwt_single limit 10;
select * from refseq_promoters limit 10;
select * from single_mappers limit 10;


select count(*) 
from bowtie inner join single_mappers on
  bowtie.read_name = single_mappers.read_name
where job_id = 2
group by refseq_name, "position", strand;


select * from tags limit 10;
select count(*), sum(no_unique_tags) from tags;
select * from tags order by no_unique_reads desc limit 20;



--tags: Group together tags mapping in the same locations
-- drop table if exists tags;
create temp table tags AS (
  select
    -- tag: identifier of each locus  
    refseq_name || '_' || "position" || '_' || strand AS tag, 
    refseq_name, 
    "position", 
    strand, 
    count(distinct read_sequence) AS no_unique_reads
  from bowtie
  where job_id = 2
  group by refseq_name, "position", strand
  );
alter table tags add constraint "PK_tags__tag" primary key (tag);
