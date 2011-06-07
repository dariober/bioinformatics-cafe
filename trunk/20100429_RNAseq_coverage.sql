/*
   Explore read coverage from RNAseq across known transcripts extracted from 
   GTF file.

   Pileup dataset from 9/4/2010 (RNAseq CTRL and LPS). Tophat
   run from 22/1/2010 (selected reads aligned with the aid of GTF file)

   The files imported by read_table() have been produced by 
   20100429_RNAseq_CTRL_annotate_gtf-1.2.py and 20100429_RNAseq_LPS_annotate_gtf-1.2.py
   
*/

-- Import annotated GTF
select read_table($$ file:'C:/Tritume/20100409_RNAseq_CTRL_sscrofa9.56.cov', dest_table: 'tmp_ss9_56_gtf_coverage_ctrl', 
    data_type: 'varchar, varchar, int, int, int, varchar, int, text', header:True, skip:1, limit:-1, overwrite:True  $$);
select names($$ t:'tmp_ss9_56_gtf_coverage_ctrl', which: [8], rename: ['ctrl_nreads'], echo:False $$);
select names($$ t:'tmp_ss9_56_gtf_coverage_ctrl', which: [7], rename: ['ctrl_sum'], echo:False $$);
select * from tmp_ss9_56_gtf_coverage_ctrl order by transcript_id, exon_number limit 1000;

select read_table($$ file:'C:/Tritume/20100409_RNAseq_LPS_sscrofa9.56.cov', dest_table: 'tmp_ss9_56_gtf_coverage_lps', 
    data_type: 'varchar, varchar, int, int, int, varchar, int, text', header:True, skip:1, limit:-1, overwrite:True  $$);
select names($$ t:'tmp_ss9_56_gtf_coverage_lps', which: [8], rename: ['lps_nreads'], echo:False $$);
select names($$ t:'tmp_ss9_56_gtf_coverage_lps', which: [7], rename: ['lps_sum'], echo:False $$);


-- Cross tab cufflink RPKM
-- drop table tmp_rpkm ;
select cross_tab($$select transcript_id, fpkm from cufflinks_transcript_gtf where feature like 'transcript' and source in('20100317_RNAseq_CTRL','20100317_RNAseq_LPS') $$, 'tmp_rpkm');
select * from tmp_rpkm where transcript_id like 'ENSSSCT00000001533' limit 10;
select names($$ t:'tmp_rpkm' $$)  -- ""transcript_id", "20100317_RNAseq_CTRL", "20100317_RNAseq_LPS""

select count(*) from tmp_rpkm; -- 9948 transcripts


-- drop table tmp_ss9_56_gtf_coverage;
create table tmp_ss9_56_gtf_coverage AS(
  -- Combine the coverage tables
  select distinct 
         tmp_ss9_56_gtf_coverage_ctrl.transcript_id, 
         tmp_ss9_56_gtf_coverage_ctrl.rname, 
         tmp_ss9_56_gtf_coverage_ctrl.exon_number,
         tmp_ss9_56_gtf_coverage_ctrl.strand,
         tmp_ss9_56_gtf_coverage_ctrl.f_start, 
         tmp_ss9_56_gtf_coverage_ctrl.f_end,
         "20100317_RNAseq_CTRL" AS ctrl_rpkm,
         "20100317_RNAseq_LPS" AS lps_rpkm,
         tmp_ss9_56_gtf_coverage_ctrl.ctrl_sum,
         lps_sum, 
         array_to_string((string_to_array(ctrl_nreads, ',')),',') AS ctrl_cov_5p3p, 
         array_to_string((string_to_array(lps_nreads, ',')),',') AS lps_cov_5p3p
  from tmp_ss9_56_gtf_coverage_ctrl inner join 
    tmp_ss9_56_gtf_coverage_lps on 
      tmp_ss9_56_gtf_coverage_ctrl.transcript_id = tmp_ss9_56_gtf_coverage_lps.transcript_id and
      tmp_ss9_56_gtf_coverage_ctrl.exon_number = tmp_ss9_56_gtf_coverage_lps.exon_number
    -- Add RPMKs
    left join tmp_rpkm on tmp_ss9_56_gtf_coverage_ctrl.transcript_id = tmp_rpkm.transcript_id
  order by tmp_ss9_56_gtf_coverage_ctrl.transcript_id, tmp_ss9_56_gtf_coverage_ctrl.exon_number
  );
select names($$ t:'tmp_ss9_56_gtf_coverage' $$)   

-- Some summary check
select 
  count(*) AS tot_num_exons, 
  count(distinct transcript_id) AS no_transcripts,
  (select count(distinct transcript_id) from tmp_ss9_56_gtf_coverage where ctrl_rpkm > 0) AS no_ctrl_detected,
  (select count(distinct transcript_id) from tmp_ss9_56_gtf_coverage where lps_rpkm > 0) AS no_lps_detected,
  (select count(distinct transcript_id) from tmp_ss9_56_gtf_coverage where lps_rpkm > 0 and ctrl_rpkm > 0) AS no_int_detected,
  (select count(distinct transcript_id) from tmp_ss9_56_gtf_coverage where ctrl_rpkm > 2) AS no_ctrl_rpkm2,
  (select count(distinct transcript_id) from tmp_ss9_56_gtf_coverage where lps_rpkm > 2) AS no_lps_rpkm2,
    (select count(distinct transcript_id) from tmp_ss9_56_gtf_coverage where lps_rpkm > 2 and ctrl_rpkm > 2) AS no_int_detected,
  (select count(distinct transcript_id) from tmp_ss9_56_gtf_coverage where ctrl_rpkm > 5) AS no_ctrl_rpkm5,
  (select count(distinct transcript_id) from tmp_ss9_56_gtf_coverage where lps_rpkm > 5) AS no_lps_rpkm5,
    (select count(distinct transcript_id) from tmp_ss9_56_gtf_coverage where lps_rpkm > 5 and ctrl_rpkm > 5) AS no_int_detected,
  (select count(distinct transcript_id) from tmp_ss9_56_gtf_coverage where ctrl_rpkm > 10) AS no_ctrl_rpkm10,
  (select count(distinct transcript_id) from tmp_ss9_56_gtf_coverage where lps_rpkm > 10) AS no_lps_rpkm10,
  (select count(distinct transcript_id) from tmp_ss9_56_gtf_coverage where lps_rpkm > 10 and ctrl_rpkm > 10) AS no_int_detected
from tmp_ss9_56_gtf_coverage;

---------------------------[ Stitch exons togheter into transcripts ]----------

create temp table ensembl AS(
  -- Make a table of transcripts and their gene names
  select distinct array_to_string(regexp_matches(attributes, '.*transcript_id (.*?);.*'), '')::varchar AS transcript_id,
         array_to_string(regexp_matches(attributes, '.*gene_name (.*?);.*'), '|')::varchar AS gene_name
  from sus_scrofa_sscrofa9_56_gtf
  );

select distinct 
         array_to_string(regexp_matches(attributes, '.*transcript_id (.*?);.*'), '')::varchar AS transcript_id,
         array_to_string(regexp_matches(attributes, '.*gene_name (.*?);.*'), '|')::varchar AS gene_name,
         f_start, f_end,
         array_to_string(regexp_matches(attributes, '.*exon_number (.*?);.*'), '|')::int AS exon_number
from sus_scrofa_sscrofa9_56_gtf 
limit 10;

  
select count(*) from ensembl;
select * from ensembl where gene_name like 'TNFA_PIG';

-- drop table tmp_transcript_cov ;
create table tmp_transcript_cov AS(
    select 
      t1.transcript_id,
      gene_name,
      rname, 
      strand, 
      min(f_start) AS t_start, 
      max(f_end) AS t_end,
      sum((1 + f_end) - f_start) AS transcript_lenght,
      ctrl_rpkm, 
      lps_rpkm,
      sum(ctrl_sum) AS ctrl_sum, 
      sum(lps_sum) AS lps_sum,      
      array_to_string(array_agg(ctrl_cov_5p3p), ',') AS ctrl_transcript_cov_5p3p,
      array_to_string(array_agg(lps_cov_5p3p), ',') AS lps_transcript_cov_5p3p
    from (
        select transcript_id, 
               rname, 
               exon_number, 
               strand, 
               f_start, f_end, 
               ctrl_rpkm, lps_rpkm,
               ctrl_sum, 
               lps_sum,
               -- Add a string of zeros where an exon is not covered at all.
               CASE WHEN ctrl_cov_5p3p is null THEN rtrim(repeat('0,', (f_end - f_start)+1), ',') ELSE  ctrl_cov_5p3p END AS ctrl_cov_5p3p,
               CASE WHEN lps_cov_5p3p is null THEN rtrim(repeat('0,', (f_end - f_start)+1), ',') ELSE  lps_cov_5p3p END AS lps_cov_5p3p
        from tmp_ss9_56_gtf_coverage 
        where ctrl_rpkm is not null or lps_rpkm is not null
        order by transcript_id, exon_number
        ) as t1 left join ensembl on ensembl.transcript_id = t1.transcript_id
    group by t1.transcript_id, gene_name, rname, strand, ctrl_rpkm, lps_rpkm
    order by t1.transcript_id
);
select count(*) from tmp_transcript_cov ;
select * from tmp_transcript_cov limit 100;
select * from tmp_transcript_cov where transcript_id like 'ENSSSCT00000001533';


select transcript_id, array_to_string(array_agg(ctrl_cov_5p3p), ',')
from (select transcript_id, exon_number, ctrl_cov_5p3p from tmp_ss9_56_gtf_coverage where ctrl_rpkm is not null order by transcript_id, exon_number limit 1000) as t1
group by transcript_id
limit 10;


------------------------------------[ Tritume ]--------------------------------

select * from cufflinks_transcript_gtf where transcript_id like 'ENSSSCT00000001533';


select transcript_id, exon_number, ctrl_cov_5p3p from tmp_ss9_56_gtf_coverage where transcript_id like 'ENSSSCT00000001125' order by transcript_id, exon_number limit 1000;

select * from tmp_ss9_56_gtf_coverage where transcript_id like 'ENSSSCT00000000004' limit 100 ;

drop table tmp_ss9_56_gtf_coverage;

select max(f_end - f_start) from tmp_ss9_56_gtf_coverage limit 10;


drop table ss9_56_gtf_coverage
select sum(ctrl_sum) AS ctrl_lib, sum(lps_sum) AS lps_lib from ss9_56_gtf_coverage;

-- drop table tmp_first_exon;
create table tmp_first_exon AS(
  " Extract  "
  select transcript_id, rname, exon_number, f_start, f_end, 
      array_to_string((string_to_array(lps_nreads, ','))[1:100],',') AS lps_nreads, 
      array_to_string((string_to_array(lps_nreads, ','))[1:100],',') AS ctrl_nreads
  from ss9_56_gtf_coverage
  where exon_number = 1 and (lps_nreads is not null OR ctrl_nreads is not null)
  );
 
drop table tmp_ss9_56_gtf_coverage_ctrl;
drop table tmp_ss9_56_gtf_coverage_lps;
drop table tmp_first_exon;

"transcript_id, rname, exon_number, f_start, f_end, ctrl_nreads, lps_nreads"
select * from tmp_ss9_56_gtf_coverage_ctrl limit 10;

select sum(fpkm) from cufflinks_transcript_gtf where source like '20100317_RNAseq_LPS' and feature like 'transcript';

drop table tmp_ss9_56_gtf_coverage;

select * from tmp_edger_transcript where transcript_id like 'ENSSSCT00000001533'

select names($$ t:'ss9_56_gtf_coverage', quote: '' $$);


select * from sus_scrofa_sscrofa9_56_gtf where attributes like '%ENSSSCT00000020953%'