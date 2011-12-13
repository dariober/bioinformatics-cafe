/*
  Identify CAGE tags mapping in the vicinity of exons.
  CAGE tags have been mapped to the pig genome and they are in SAM format.
  Exons are pulled out from a GTF file (version 9.56, S. scrofa).
  CAGE tags are searched within n bases from the start of each exon (where
  start is the "end" in the gtf file if the transcript is on - strand).
*/

---------------------------[ Performance parameters ]---------------------------

show work_mem;
set work_mem = '256MB';

--------------------------------[ Prepare GTF ]---------------------------------
/*
   Extract from the GTF file the features for which you want to see which CAGE
   tags are mapped in the nearby. E.g. extract exons.
   For each feature (exon) define an interval surrounding its starting position (i.e.
   surrounding the "end" if sitting on the '-' strand, otherwise surrounding the
   "start").
   Take the unique starting positions in a separate table and assign it a primary 
   key (this will be joined to the cage tss)
*/

drop table gtf_1;
create temp table gtf_1 AS(
  -- Pull out all the exons from the gtf file
  select distinct 
      gtf_attribute(attributes, 'transcript_id') AS transcript_id,
      gtf_attribute(attributes, 'exon_number') AS exon_number,
      rname,
      strand,
      f_start,
      f_end,
      CASE WHEN strand = '+' THEN f_start WHEN strand = '-' THEN f_end END AS tss_start -- Start of the exon (start is the "end" if strand is '-')
  FROM sus_scrofa_sscrofa9_56_gtf
  WHERE feature like 'exon'
  ); -- 6.3 sec for 168515 rows
analyze gtf_1;
select count(*) from gtf_1;

drop table if exists gtf;
create temp table gtf AS(
    -- Expand each starting position so that each row is a position surrounding each starting point
    select 
       rname || '_' || (tss_start + range) || '_' || strand AS tss_id,
       gtf_1.*, 
       tss_start + range AS tss_range 
    from gtf_1, (select generate_series(-100, +100) AS range) AS t1 -- Change the range in generate_series() to extend the range surrounding each tss
--    order by tss_id
    ); -- 10 secs for generate_series(-10, +10); 90 sec. for generate_series(-100, +100)
create index index_tss_id ON gtf (tss_id); -- 75 sec. when generate_series(-10, +10) and 3538815 returned. 755172 ms. for generate_series(-100, +100)
analyze gtf; -- 39278 ms. for generate_series(-100, +100)

drop table if exists tss_id;
-- This table contains all the unique tss surrounding each exon
create temp table tss_id AS(
    select distinct tss_id from gtf
    ); -- 38 sec
alter table tss_id add primary key (tss_id); -- 66 sec. -- create + alter for generate_series(-100, +100) = 1,275,568 ms (21.2 min.)
analyze tss_id; -- 77165 ms. for generate_series(-100, +100)

select count(*) from gtf;

-----------------------------[ Prepare CAGE tags ]-----------------------------

drop table cage_tags;
create temp table cage_tags AS(
  -- Extract the mapped CAGE tags and create an identifier for each tag (tss_id,
  -- which will be joined to the gtf tss_id)
  select 
    rname || '_' || CASE WHEN flag = 0 THEN pos + length(seq)
                         WHEN flag = 16 THEN pos END || '_' || CASE 
                         WHEN flag = 0 THEN '+' 
                         WHEN flag = 16 THEN '-' END AS tss_id,
    sam_bowtie_human_vs_pig.*
  from sam_bowtie_human_vs_pig
  where flag in (0, 16)
  ); -- c.ca 7.83 min.; 6.75 after changing log.
create index index_cage_tags on cage_tags (tss_id); -- 651469 ms; 10.85 min. c.ca 27 million rows.
analyze cage_tags; -- 64140 ms.

select * from cage_tags limit 1000;

-- Option 1: Join the range of tss produced in gtf.tss_id directly to the  
--           tag table

select * from cage_tags limit 10; 
create temp table cage_to_gtf AS(
  select cage_tags.tss_id,
         qname,
         flag,
         rname,
         pos,
         mapq,
         cigar,
         tags
  from cage_tags inner join tss_id ON cage_tags.tss_id = tss_id.tss_id 
  ); -- c.ca 5678141 ms. (c.ca 94 min for 19 million tags and 3538815 tss_id from gtf
     -- c.ca 2 min. (!) after having indexed cage_tags, analyzed, and set work_mem = '512MB'
     -- 545342 ms (9.08 min) with generate_series(-100, +100)

select * from cage_to_gtf limit 10;

     
------------------------------[ exons tagged by CAGE ]-------------------------

-- analyze gtf;
select cage_id_gtf_id.tss_id, gtf.transcript_id, exon_number, rname, strand, f_start, f_end
from gtf inner join cage_id_gtf_id on gtf.tss_id = cage_id_gtf_id.tss_id; -- few seconds (3-17 sec).

select cage_tags.* from cage_tags inner join cage_id_gtf_id on cage_id_gtf_id.tss_id = cage_tags.tss_id;

--___________________________________________________________________________--
--
--                                 Tritume
--___________________________________________________________________________--

"Hash Join  (cost=1012380.11..3 841 764.93 rows=23 958 740 width=79)"

-- Option 2: Produce table of unique CAGE tags ids, primary key this table, 
--           map primary key against gtf.tss_id.

drop table cage_tags_ids;
create temp table cage_tags_ids AS(
    -- Unique cage tags: add primary key.
    select distinct tss_id from cage_tags
    ); -- 198875 ms (3.31 min)
alter table cage_tags_ids add constraint pk_cage_tags_ids primary key (tss_id); -- 369671 ms. (6.16 min) for 19223754 rows output'd
analyze cage_tags_ids;
select count(*) from cage_tags_ids;


-- analyze cage_tags;

drop table cage_id_gtf_id;
create temp table cage_id_gtf_id AS(
  select cage_tags_ids.tss_id 
  from cage_tags_ids inner join tss_id on 
     cage_tags_ids.tss_id = tss_id.tss_id
  ); -- 103469 ms.; 1.72 min (36475 rows output'd)
alter table cage_id_gtf_id add constraint pk_cage_id_gtf_id primary key (tss_id);
analyze cage_id_gtf_id;



select * from gtf limit 10;
select names('sam_bowtie_human_vs_pig')
"dataset_id", "qname", "flag", "rname", "pos", "mapq", "cigar", "mrnm", "mpos", "isize", "seq", "qual", "tags"

select count(*) from tss_id;



limit 100;



CREATE TABLE "Pigs".sus_scrofa_sscrofa9_56_gtf
(
  rname character varying,
  source character varying,
  feature character varying,
  f_start integer,
  f_end integer,
  score character varying,
  strand character varying,
  frame character varying,
  attributes character varying
)
WITH (
  OIDS=FALSE
);
ALTER TABLE "Pigs".sus_scrofa_sscrofa9_56_gtf OWNER TO dberaldi;

set work_mem = '512MB';