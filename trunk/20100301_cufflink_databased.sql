select read_table('C:/Documents and Settings/postgres/My Documents/Tritume/20100301_LPS_transcripts.expr', 'cufflink_transcript_expr', True, E'\t', 0, '');
update cufflink_transcript_expr set trans_id = replace(trans_id, 'CUFF.', 'RNAseq_LPS.') ;
alter table cufflink_transcript_expr add constraint pk_trans_id primary key (trans_id);

copy cufflink_transcript_expr from 'C:/Documents and Settings/postgres/My Documents/Tritume/20100301_CTRL_transcripts.expr' with csv header delimiter E'\t';
update cufflink_transcript_expr set trans_id = replace(trans_id, 'CUFF.', 'RNAseq_CTRL.') ;

select * from cufflink_transcript_expr 
where trans_id like 'RNAseq_CTRL%'
limit 100;


select * from tmp_cufflink_lps 
limit 100;

select count(*) from tmp_cufflink_lps limit 10;