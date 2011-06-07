/*
  Annotation of the snowball array consisted in annotating FASTA files (non_PSR_fastas.zip)
  against refseq. Now we want to link FASTA names to probe set IDs
*/

-- This file from "Snowball_array_annotation_8-3-11.xls" from Chris Tuggle.
select read_table($$ file:'F:/data/20110504_snowball_annotation/Snowball_array_annotation_8-3-11.txt', 
    temp: True, header:True, table: 'array_ann' $$);
create unique index ind_symbolprobe on array_ann ("Symbol;Probe");
create unique index ind_snowball_id on array_ann ("SNOWBALL_ID") where "SNOWBALL_ID" is not null;
create unique index ind_seqname on array_ann ("Seqname") where "SNOWBALL_ID" is not null;
select * from array_ann limit 1000;

select count(*), count(distinct qseqid), count(distinct "SNOWBALL_ID") from array_ann inner join non_psr_fastas on non_psr_fastas.qseqid = array_ann."SNOWBALL_ID"; -- 47722 rows

select count(*) from array_ann where "SNOWBALL_ID" is null;

select * from array_ann where "SNOWBALL_ID"||'' != ;

-- This file from Kim 
select read_table($$ file:'F:/data/20110504_snowball_annotation/Pig_array_annotation.txt', 
    temp: True, header:True, table: 'pig_array' $$);
select * from pig_array limit 1000;

select "Probeset", "Symbol", qseqid, symbol, refseq_id, species_name
from pig_array inner join blastn_nonpsr_refseqrna_best on substring("Probeset" from 1 for length("Probeset")-3 ) = qseqid 
where (symbol is null or symbol like 'LOC%') and "Symbol" is not null and "Symbol" != '-' -- 122
limit 1000;

select * from pig_array where "Probeset" like 'SNOWBALL_012304%'
select * from blastn_nonpsr_refseqrna_annotated where qseqid = 'SNOWBALL_000141'

select 600 / 29000.0