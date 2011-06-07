create temp table ss9_ttome AS(select * from mart_tome where genome like 'ss9');

select distinct 
  tmp_inserts_feb2010.qname, flag, 
  seq, 
  ss9_ttome.*, 
  tmp_inserts_feb2010.pos, 
  (CASE WHEN strand = 't' THEN pos - transcript_start WHEN strand = 'f' THEN transcript_end - pos END) AS tss_distance
from ss9_ttome inner join tmp_inserts_feb2010 on chr = rname
where 
  (CASE WHEN strand = 't' THEN pos - transcript_start WHEN strand = 'f' THEN transcript_end - pos END) between -1000 and +1000 and
  strand = (CASE WHEN flag = 0 THEN 't' WHEN flag = 16 THEN 'f' END)::boolean
order by tss_distance;

--------------------------------[ Tritume ]------------------------------------

select * from mart_tome limit 10;
select * from tmp_inserts_feb2010;

select get_column_names('tmp_inserts_feb2010');

select distinct chr from ss9_ttome;
select distinct rname from tmp_inserts_feb2010;

"['genome', 'ensembl_gene_id', 'ensembl_transcript_id', 'chr', 'transcript_start', 'transcript_end', 'strand']"

"['qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'mrnm', 'mpos', 'isize', 'seq', 'qual', 'opt_field_1', 'opt_field_2', 'opt_field_3']"