---------------------------------[ hg95ctrls vs Pig ]--------------------------
drop table if exists mapped;
create temp table mapped AS(
  select hg95ctrls.*
  from tmp_hg95ctrls_ss9 inner join hg95ctrls on
  read_name = id
  );

select 
  (select sum(cd14) from hg95ctrls) AS total_reads_cd14,
  (select sum(cd14) from mapped)    AS mapped_reads_cd14,
  (select sum(cd14) from mapped)/(select sum(cd14) from hg95ctrls)::numeric AS prop_mapped
  ;

select 
  (select sum(control) from hg95ctrls) AS total_reads_control,
  (select sum(control) from mapped)    AS mapped_reads_control,
  (select sum(control) from mapped)/(select sum(control) from hg95ctrls)::numeric AS prop_mapped
  ;

select 
  (select sum(total) from hg95ctrls) AS total_reads_total,
  (select sum(total) from mapped)    AS mapped_reads_total,
  (select sum(total) from mapped)/(select sum(total) from hg95ctrls)::numeric AS prop_mapped
  ;

---------------------------------[ hg95ctrls vs Human]-------------------------

drop table if exists mapped;
create temp table mapped AS(
  select hg95ctrls.*
  from tmp_hg95ctrls_hg19 inner join hg95ctrls on
  read_name = id
  );

-- Library cd14 only
select 
  (select sum(cd14) from hg95ctrls) AS total_reads_cd14,
  (select sum(cd14) from mapped)    AS mapped_reads_cd14,
  (select sum(cd14) from mapped)/(select sum(cd14) from hg95ctrls)::numeric AS prop_mapped
  ;

-- library control only
select 
  (select count(distinct id) from hg95ctrls where control > 0) AS distinct_reads_in_control,
  (select sum(control) from hg95ctrls) AS total_reads_control,
  (select sum(control) from mapped)    AS mapped_reads_control,
  (select sum(control) from mapped)/(select sum(control) from hg95ctrls)::numeric AS prop_mapped
  ;



-- Both libraries
select 
  (select sum(total) from hg95ctrls) AS total_reads_total,
  (select sum(total) from mapped)    AS mapped_reads_total,
  (select sum(total) from mapped)/(select sum(total) from hg95ctrls)::numeric AS prop_mapped
  ;