drop table sam_rebase;

select read_sam($$ file:'D:/Tritume/20100706_RNAseq_CTRL_repbase_se.aligned.sam', dest_table: 'sam_bowtie_rebase', dataset_id: 'CTRL_se' $$)
select read_sam($$ file:'D:/Tritume/20100706_RNAseq_LPS_repbase_se.aligned.sam', dest_table: 'sam_bowtie_rebase', dataset_id: 'LPS_se', append:True $$)

comment on table sam_bowtie_rebase is 'SAM alignment produced by Bowtie. Reference is the pig database of repeats. See script 20100706_sam_bowtie_rebase.sql'

drop table tmp_flag_count;
select cross_tab('select flag, dataset_id, count(flag) AS flag_count from sam_bowtie_rebase group by flag, dataset_id order by dataset_id, flag', 'tmp_flag_count');
select * from tmp_flag_count;

drop table tmp_rname_count;
select cross_tab('select rname, dataset_id, count(qname) from sam_bowtie_rebase group by rname, dataset_id', 'tmp_rname_count');
copy tmp_rname_count to 'D:/Tritume/RNAseq_repbase.txt' with csv header delimiter E'\t' ;
select * from tmp_rname_count;

-- Number of repeats hit by at least 1 read
select distinct dataset_id, count(distinct rname) from sam_bowtie_rebase group by dataset_id;


-------------------------[ This repeat found diff.expr. by DESeq ]-------------
--
--LTR16A1#LTR/ERVL_RepbaseID:_LTR16A1

select * from sam_bowtie_rebase where rname like 'LTR16A1#LTR/ERVL_RepbaseID:_LTR16A1';
select qname, count(qname) from sam_bowtie_rebase where rname like 'LTR16A1#LTR/ERVL_RepbaseID:_LTR16A1' group by qname;

select * from ucsc_repeat_masker_pig where rep_name like 'LTR16A1';





-----------------------------------[ Tritume ]---------------------------------

select read_sam($$ file:'D:/Tritume/20100629_RNAseq_CTRL_repbase.aligned.sam', dest_table: 'sam_rebase', dataset_id: 'CTRL' $$)

select distinct rep_name from ucsc_repeat_masker_pig;

select * from sam_bowtie_rebase limit 10;

select distinct mrnm from sam_rebase;