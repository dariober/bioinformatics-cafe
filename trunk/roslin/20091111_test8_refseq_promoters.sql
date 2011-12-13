-- Maps cage reads from pig AM Test8 against the vicinity of 5' end 
-- of human RefSeq genes (hg19 build from Genome Browser).
-- CAGE reads filtered to have only single mappers


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


-----------------[ bwt_single: CAGE tags Test8: single mappers ]---------------

-- drop table if exists single_mappers ;
create temp table single_mappers AS ( -- List of tags with one hit only
  select read_name
  from bowtie 
  where job_id = 3 -- bowtie job: Test8 reads with quality score always >= 20 against
                   -- human genome hg37
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
  where job_id = 3
  );

-----------------------[ Single mappers near promoters ]-----------------------

-- Make sure the table with tags has chromosome names and positions consistent 
-- with the table of refseq genes!

-- drop table if exists out_test8_promoters ;
create table out_test8_promoters AS (
  select 
    bwt_single.read_name, 
    bwt_single.chr, 
    bwt_single.strand,
    tss,
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

comment on table out_test8_promoters is 
'CAGE tags from pig library test8 mapping near the 5'' end of human refseq genes (UCSC). 
Single mappers only. 
Table generated with 20091111_test8_refseq_promoters.sql';

-------------[ Group reads mapping to the same position into tags ]------------

-- drop view if exists qry_test8_promoters ;
create view qry_test8_promoters AS (
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
  from out_test8_promoters
  );
  
comment on view qry_test8_promoters is 
'CAGE tags from pig Test8 mapping near 5'' end of human refseq genes. 
Single mappers only. 
Reads mapping to the same position grouped into tags.
Query executed with 20091111_test8_refseq_promoters.sql';

-----------------------------[ Summary stats ]---------------------------------

-- Reads in the vicinity of a known 5' end
select 
  (select count(*) from single_mappers) AS total_single_mappers,
  count(*) AS no_hits,
  count(distinct read_name) AS mapped_reads,
  (select count(distinct read_name) from out_test8_promoters) / (select count(*) from single_mappers)::numeric AS percentage_mapped,
  count(distinct seq_name) AS no_unique_genes_hit,
  count(distinct seq_name)::numeric / (select count(distinct seq_name) from refseq_hg19) AS percentage_genes_hit
from out_test8_promoters;

-- tags in the vicinity of a known 5' end
select 
  count(*) AS no_mapped,
  (select count(*) from (select distinct chr, strand, tss from bwt_single) as t1)
    AS total_single_tags, 
  count(*)::numeric / (select count(*) from (select distinct chr, strand, tss from bwt_single) as t1)
    AS percentage_mapped
from qry_test8_promoters;

-- Disagreement bwt transcript and tag orientation. I.e. tags mapped in one 
-- orientation map near a transcript with opposite orientation
select count(*), count(*)::numeric /  (select count(*) from single_mappers) AS percentage_opposite
from out_test8_promoters
where strand != refseq_strand;

-- Summary table of cage and gene orientation
select strand AS cage_strand, refseq_strand, 
  CASE WHEN tss_distance <= 0 THEN 'upstream' 
       WHEN tss_distance > 0 THEN 'downstream' END AS tss_rel_pos,
  count(*), 
  100*count(*)::numeric / (select count(*) from out_test8_promoters) AS percentage,
  avg(tss_distance) AS avg_tss_distance,
  stddev(tss_distance) AS stdev_tss_distance
from out_test8_promoters
group by strand, refseq_strand, tss_rel_pos
order by percentage, strand, refseq_strand, tss_rel_pos;

-- RefSeq genes with associated a tag (tagged refseq genes)
select 
  (select count(distinct seq_name) from refseq_hg19) AS no_refseq_genes,
  (select count(distinct seq_name) from out_test8_promoters ) 
    AS no_tagged_refseq_genes,
  (select count(distinct seq_name) from out_test8_promoters ) / 
    (select count(distinct seq_name) from refseq_hg19)::numeric 
    AS percentage_tagged
  ;

-- How many of the reads mapping near genes are mapped against Sscorfa9.53?
select count(distinct bowtie.read_name) 
from bowtie inner join out_test8_promoters 
  on out_test8_promoters.read_name = bowtie.read_name 
where job_id = 1;

-------------------------------[ Tritume ]-------------------------------------

select count(*) from bwt_single;
select * from bwt_single limit 10;
select * from refseq_promoters limit 10;
select count(*) from refseq_promoters limit 10;
