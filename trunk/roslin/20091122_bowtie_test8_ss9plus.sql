-- Select reads with a *single* best alignment

-- drop table if exists bowtie_plus;
create temp table bowtie_plus AS (
  select 
    read_name, 
    strand, 
    refseq_name, 
    "position", 
    read_sequence, 
    list_mism, 
    count_pattern(list_mism, '>') AS no_mism
  from bowtie
  where job_id = 4
  );


-- drop table if exists bowtie4_mism_stratum;
create temp table bowtie4_mism_stratum AS (
  select 
    read_name, 
    no_mism,
    count(count_pattern(list_mism, '>')) AS no_hits_stratum -- Number of hits with "no_mism" no of mismatches
  from bowtie_plus
  group by read_name, no_mism
  );

-- drop table if exists bowtie4;
create temp table bowtie4 AS (
  select bowtie_plus.*, bowtie4_mism_stratum.no_hits_stratum
  from bowtie_plus inner join bowtie4_mism_stratum on
    bowtie4_mism_stratum.read_name = bowtie_plus.read_name AND
    bowtie4_mism_stratum.no_mism = bowtie_plus.no_mism
  order by bowtie_plus.read_name, bowtie_plus.no_mism
  );

-- drop table if exists best_stratum;
create temp table best_stratum AS (
  select read_name, min(no_mism) AS best_stratum
  from bowtie4
  group by read_name
  );

-- drop table if exists best_single;
create table best_single AS (
  select bowtie4.*
  from bowtie4 inner join best_stratum on
    bowtie4.read_name = best_stratum.read_name AND
    bowtie4.no_mism = best_stratum.best_stratum
  where no_hits_stratum = 1
  );

--------------------------------[ Stats ]--------------------------------------

-- Total reads mapped and single, best mappers
select 
  (select count(distinct read_name) from bowtie4) AS no_aligned_reads,
  (select count(*) from best_single) AS no_single_mappers,
  (select count(*) from best_single)::numeric / (select count(distinct read_name) from bowtie4) AS perc_single_mappers
  ;

-- Overwiev of alignment quality: Mismatches in single, best mappers
select 
  no_mism, 
  count(*) AS no_reads, 
  count(*)::numeric / (select count(*) from best_single) AS mism_frequency
from best_single 
group by no_mism;

-- Where single, best mappers align, chromosome by chromosome
select 
  refseq_name AS reference_chr, 
  count(*) AS no_reads, 
  count(*)::numeric / (select count(*) from best_single) AS proportion
from best_single
group by refseq_name
order by no_reads;

--------------------------------[ Tritume ]------------------------------------



select no_mism, count(distinct read_name)
from bowtie4
where refseq_name = 'gi|555853|gb|U13369.1|HSU13369'
group by no_mism;

select no_mism, count(distinct read_name), count(distinct read_name)::numeric / (select count(distinct read_name) from bowtie4)
from bowtie4
group by no_mism;

select refseq_name, count(*) AS no_mapped_reads
from best_single  
group by refseq_name
order by no_mapped_reads;

select * from bowtie4 limit 200;


-- Select best strata
-- drop table if exists bowtie4_best_strata;
create temp table bowtie4_best_strata AS(
  select 
    read_name,
    min(no_mism) AS best_stratum
  from bowtie4_mism_stratum 
  group by read_name
  );


-- Which reads have a single hit in the best startum?
-- drop table if exists bowtie4_best_strata2;
create temp table bowtie4_best_strata2 AS(
  select
    bowtie4_best_strata.read_name,
    bowtie4_best_strata.best_stratum,
    bowtie4_mism_stratum.no_hits_stratum
  from bowtie4_best_strata inner join bowtie4_mism_stratum on
    bowtie4_best_strata.read_name = bowtie4_mism_stratum.read_name AND
    bowtie4_best_strata.best_stratum = bowtie4_mism_stratum.no_mism
  where bowtie4_mism_stratum.no_hits_stratum = 1
  );

select * from bowtie4_mism_stratum order by read_name, no_mism limit 100;
select * from bowtie4_best_strata order by read_name limit 100; --"EBRI093151:6:1:1133:1728#0/1"
select * from bowtie4_best_strata2 order by read_name limit 100; --"EBRI093151:6:1:1133:1728#0/1"


select count(read_name) from bowtie4_best_strata2; -- 233527

select count(*) from test8_ss9_single;

select * from bowtie where job_id = 4 and read_name = 'EBRI093151:6:1:1129:649#0/1';

select * from bowtie where job_id = 4 and read_name = 'EBRI093151:6:1:1130:457#0/1';

alter table test8_ss9_single rename to bowtie4_mism_stratum;