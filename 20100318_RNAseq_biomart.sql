truncate table cufflinks_transcript_gtf;
-- transcript.gtf RNAseq_CTRL and RNAseq_LPS assembled against GTF file sscrofa9.56
select read_table($$ file:'C:/Tritume/20100317_RNAseq_CTRL/transcripts.gtf.split', 
  dest_table:'cufflinks_transcript_gtf', header: True, overwrite:True,
  data_type: 'varchar, varchar, varchar, int, int, double precision, varchar, varchar, varchar, varchar, varchar, double precision, double precision, double precision, double precision, double precision' 
  $$);
select read_table($$ file:'C:/Tritume/20100317_RNAseq_LPS/transcripts.gtf.split', 
  dest_table:'cufflinks_transcript_gtf', header: True, append: True $$);

select fpkm_ctrl.*, fpkm_lps.fpkm_lps from 
  (select distinct transcript_id, fpkm AS fpkm_ctrl from cufflinks_transcript_gtf where source like '20100317_RNAseq_CTRL') AS fpkm_ctrl 
  INNER JOIN 
  (select distinct transcript_id, fpkm AS fpkm_lps from cufflinks_transcript_gtf where source like '20100317_RNAseq_LPS') AS fpkm_lps ON
  fpkm_ctrl.transcript_id= fpkm_lps.transcript_id;

-- drop table fpkm;
create temp table fpkm AS(
  select fpkm_ctrl.*, fpkm_lps.fpkm_lps,  log(2, (fpkm_lps.fpkm_lps/fpkm_ctrl.fpkm_ctrl)::numeric) AS fold_change from 
    (select distinct transcript_id, fpkm AS fpkm_ctrl from cufflinks_transcript_gtf where source like '20100317_RNAseq_CTRL') AS fpkm_ctrl 
    INNER JOIN 
    (select distinct transcript_id, fpkm AS fpkm_lps from cufflinks_transcript_gtf where source like '20100317_RNAseq_LPS') AS fpkm_lps ON
    fpkm_ctrl.transcript_id= fpkm_lps.transcript_id
    order by fold_change
  );
alter table fpkm add constraint pk_fpkm primary key  (transcript_id);
copy (select transcript_id from fpkm) to 'C:/Tritume/RNAseq_transcript_id.txt' with csv delimiter E'\t';

create temp table overexpr AS(
  select * from fpkm where fold_change >= (select fold_change from fpkm where transcript_id like 'ENSSSCT00000001533')
  );

select read_table($$file:'C:/Downloads/mart_export_RNAseq.txt', dest_table:'tmp_mart_rnaseq', limit:-1, header:True, overwrite:True$$)
select set_column_names('tmp_mart_rnaseq', $$['ensembl_gene_id', 'ensembl_transcript_id', 'description', 'associated_gene_name', 'associated_transcript_name', 'go_term_name_bp', 'go_term_name_cc', 'go_term_name_mf']$$)
select * from tmp_mart_rnaseq limit 10;

select distinct overexpr.*, tmp_mart_rnaseq.description, associated_gene_name
from overexpr inner join tmp_mart_rnaseq on ensembl_transcript_id = transcript_id
order by fold_change;

select transcript_id, log(2, (fpkm_lps/fpkm_ctrl)::numeric) AS expr_ratio from fpkm order by expr_ratio desc limit 10;
select *, log(2, (fpkm_lps/fpkm_ctrl)::numeric) from fpkm where transcript_id like 'ENSSSCT00000001533'; -- TNF-alpha
select *, log(2, (fpkm_lps/fpkm_ctrl)::numeric) from fpkm where transcript_id like 'ENSSSCT00000013865'; -- HPRT
select *, log(2, (fpkm_lps/fpkm_ctrl)::numeric) from fpkm where transcript_id like 'ENSSSCT00000008324'; -- ACTB



select * from cufflinks_transcript_gtf where transcript_id like 'ENSSSCT00000008324';
select * from mart_tome where ensembl_transcript_id like 'ENSSSCT00000008324';

select * from cufflinks_transcript_gtf where transcript_id like 'ENSSSCT00000001533' and feature like 'transcript'; -- TNF-alpha
select * from cufflinks_transcript_gtf where transcript_id like 'ENSSSCT00000013865' and feature like 'transcript'; -- HPRT

(select distinct transcript_id from cufflinks_transcript_gtf where feature like 'transcript') AS transcript_id;
  
select distinct source from cufflinks_transcript_gtf limit 10;



select read_table($$ file: 'C:/Tritume/Sus_scrofa.Sscrofa9.56.gtf', dest_table: 'sscrofa9_56_gtf' $$);
-- drop table sscrofa9_56_gtf;
select * from sscrofa9_56_gtf limit 10;
    
select * from sscrofa9_56_gtf limit 10;

select read_table($$file:'C:/Tritume/20100316_cuffcompare_RNAseq/20100316_cuffcompare_RNAseq_summary.combined.gtf', dest_table: 'tmp_cuffcompare_gtf'$$);
-- drop table tmp_cuffcompare ;
select set_column_names('tmp_cuffcompare_gtf', $$['rname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']$$)
select * from tmp_cuffcompare_gtf limit 1000;

select read_table($$file:'C:/Tritume/20100316_cuffcompare_RNAseq/20100316_cuffcompare_RNAseq_summary.tracking', dest_table: 'tmp_cuffcompare_tracking'$$);
select * from tmp_cuffcompare_tracking limit 10;

--------------------------------------[ Tritume ]------------------------------
select source, count(*)from  cufflinks_transcript_gtf group by source;


select read_table($$file:'C:/Tritume/20100301_RNAseq_CTRL/transcripts.gtf', 
  dest_table:'cufflinks_transcript_expr', overwrite:True, header: True,
  dataset_id:'20100301_RNAseq_CTRL' $$)

select get_column_names('cufflinks_transcript_gtf')

select read_table($$file:'C:/Tritume/20100301_RNAseq_LPS/transcripts.expr', 
  dest_table:'cufflinks_transcript_expr', append:True, header: True,
  dataset_id:'20100301_RNAseq_LPS' $$)

update cufflinks_transcript_expr set trans_id = replace(trans_id, 'CUFF', 'RNAseq_CTRL') WHERE dataset_id like '20100301_RNAseq_CTRL';
update cufflinks_transcript_expr set trans_id = replace(trans_id, 'CUFF', 'RNAseq_LPS') WHERE dataset_id like '20100301_RNAseq_LPS';

select read_table($$file:'C:/Tritume/20100301_RNAseq_CTRL/transcripts.expr', 
  dest_table:'cufflinks_transcript_expr', overwrite:True, header: False, skip: 1000,
  dataset_id:'20100301_RNAseq_CTRL', limit: 10$$)
select * from cufflinks_transcript_expr limit 10;
select * from cufflinks_transcript_expr limit 10;
select dataset_id, count(*) from cufflinks_transcript_expr group by dataset_id;

select column_name, data_type from information_schema.columns where table_name like 'cufflinks_transcript_expr' order by ordinal_position;
(select get_column_names('cufflinks_transcript_gtf')))