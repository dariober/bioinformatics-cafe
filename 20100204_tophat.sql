---------------[ Extract junction coordinates from tophat output ]-------------

-- drop table if exists tophat_juncs;
create temp table tophat_juncs AS ( 
  select 
    tophat_job, 
    junction_name, 
    chr, 
    (block_start + (string_to_array(block_sizes, ','))[1]::int) - 1 AS j_left,
    (block_end - (string_to_array(block_sizes, ','))[2]::int) AS j_right,
    strand
  from tophat
  );
alter table tophat_juncs  add constraint PK_tophat_juncs primary key (tophat_job, chr, j_left, j_right, strand);
comment on column tophat_juncs.j_left is 'Coordinate of the last character of the left sequence to be spliced, inclusive. This column corresponds to the junction file produced by tophat when --GFF is specified';
comment on column tophat_juncs.j_left is 'Coordinate of the first character of the right sequence, inclusive. This column corresponds to the junction file produced by tophat when --GFF is specified';

create temp table mapped_juncs AS (
  select tophat_juncs.* 
  from tophat_juncs inner join sscrofa_juncs on
    tophat_juncs.chr = sscrofa_juncs.chr
    and tophat_juncs.j_left = sscrofa_juncs.j_left
    and  tophat_juncs.j_right = sscrofa_juncs.j_right
    and  tophat_juncs.strand = sscrofa_juncs.strand
-- where tophat_job like '20100115_LPS'
  order by junction_name
  );

select distinct t1.tophat_job, total_mapped, found_in_ensembl,  (total_mapped - found_in_ensembl) AS not_in_ensembl, (total_mapped - found_in_ensembl)/total_mapped::numeric AS perc_new
  from 
  (select tophat_job, count(*) AS total_mapped from tophat group by tophat_job) as t2 
  inner join 
  (select mapped_juncs.tophat_job, count(*) AS found_in_ensembl from mapped_juncs group by tophat_job) as t1
  on
  t1.tophat_job = t2.tophat_job;

select count(*)
  from
  (select * from tophat_juncs where tophat_job like '20100122_LPS_gff') as lps
  inner join
  (select * from tophat_juncs where tophat_job like '20100122_CTRL_gff') as ctrl
  on
  lps.chr = ctrl.chr and
  lps.j_left = ctrl.j_left and
  lps.j_right = ctrl.j_right and
  lps.strand = ctrl.strand
--  where ctrl.chr is null;

3445258

select tophat_job, chr, j_left, j_right from tophat_juncs

select * from tophat_juncs limit 10;
  
select get_column_names('tophat_juncs');

select tophat_job, count(*) from tophat 
group by tophat_job;

select (string_to_array(block_sizes, ','))[1] from tophat limit 10;