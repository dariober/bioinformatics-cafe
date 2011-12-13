/*
  Which reads in hg95ctrls (single mappers) can be found mapped on Sscrofa 9.53 (and give a single hit)
*/

-------------------------[ Single mappers in hg19 and ss9.53]------------------

-- drop table if exists hg_single_mappers ;
create temp table hg_single_mappers AS ( -- List of hg95ctrls tags with one hit only against hg19
  select read_name
  from bowtie 
  where job_id = 2
  group by read_name 
  having count(*) = 1
  );

-- drop table if exists ss_single_mappers ;
create temp table ss_single_mappers AS ( -- List of hg95ctrls tags with one hit only against Sscrofa 9.53
  select read_name
  from bowtie 
  where job_id = 1
  group by read_name 
  having count(*) = 1
  );

-- Intersect single mappers in hg19 and ss9 (i.e. single mappers common to human and Sus)

-- drop table if exists common_single_mappers ;
create temp table common_single_mappers AS(
  select ss_single_mappers.* 
  from ss_single_mappers inner join hg_single_mappers 
    on ss_single_mappers.read_name = hg_single_mappers.read_name
  );

-- drop table if exists tmp_ss_bwt_single ;
create table tmp_ss_bwt_single AS (   -- Annotate tags from previous query 
                                      -- with genome position and strand info from Ss9 alignment
  select 
    bowtie.* 
  from bowtie inner join common_single_mappers on 
    common_single_mappers.read_name = bowtie.read_name
  where job_id = 1
  );

-- Cluster overlapping tags in tmp_ss_bwt_single
-- drop table if exists tmp_cluster_hg95_ss;
select cluster_reads(
  'select 
    refseq_name, 
    strand, 
    position AS read_start, 
    position + length(read_sequence) AS read_end 
  from tmp_ss_bwt_single', 'tmp_cluster_hg95_ss');

-- drop table if exists tmp_ss_bwt_single ;

-- How expressed are the reads in the common set?
create temp table ss_bwt_single_expr AS(
  select tmp_ss_bwt_single.*, hg95ctrls.total
  from tmp_ss_bwt_single inner join hg95ctrls on
    tmp_ss_bwt_single.read_name = hg95ctrls.id
  );

----------------------------------[ Summary stats ]----------------------------

select 
  (select count(*) from hg_single_mappers) AS "single mappers in hg19" ,
  (select count(*) from ss_single_mappers) AS "single mappers in ss9.53",
  (select count(*) from tmp_ss_bwt_single) AS "common reads",
  (select count(distinct tag_cluster) from tmp_cluster_hg95_ss) 
    AS "No. of clusters in common set"
  ;
  
select 
  (select avg(total) from ss_bwt_single_expr),
  (select avg(total) from hg95ctrls)
  ;

------------------------------[ Tritume ]--------------------------------------
select count(*) from ss_single_mappers;
select count(*) from hg_single_mappers;
select count(*) from common_single_mappers;
select count(*) from ss_bwt_single;
select tmp_cluster_hg95_ss.* from tmp_cluster_hg95_ss limit 1000;
select count(distinct tag_cluster) from tmp_cluster_hg95_ss;

select * from ss_bwt_single_expr limit 100;

drop table tmp_cluster_hg95_ss;