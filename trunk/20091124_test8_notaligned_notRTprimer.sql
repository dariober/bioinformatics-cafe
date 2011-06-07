/* 
Sequences of high quality reads (each base quality >= 20, AGACAGCAG correct)
not mapped against sscrofa 9.53, not matching RT primers.
Grouped and counted
*/ 

-- drop table if exists test8_qs20_unaligned_sequences ;
create temp table test8_qs20_unaligned_sequences AS (
  select 
    "tmpTest8_qs20".read_sequence, 
    count(*) AS read_count
  from "tmpTest8_qs20" left join "tmpSoapTest8_qs20_r1" on
    "tmpTest8_qs20".read_name = "tmpSoapTest8_qs20_r1".read_name 
  where 
    "tmpTest8_qs20".read_sequence not like '%AAGGTCTAT%' AND -- Remove RT primers
    "tmpSoapTest8_qs20_r1".read_name is null                 -- Remove aligned reads
  group by "tmpTest8_qs20".read_sequence                
  order by read_count desc
  );

-- Write out fasta file (N.B. writefasta() can't read temporary tables)
create table tmp_notaligned AS (select * from test8_qs20_unaligned_sequences); 
select writefasta('select read_count, read_sequence from tmp_notaligned order by read_count desc limit 1000')
drop table tmp_notaligned ;

-------------------------------------[ Stats ]---------------------------------

select 
  sum(read_count) AS unidentified_reads,
  sum(read_count)::numeric / (select count(*) from "tmpTest8_qs20") AS proportion_unidentified,
  count(*) AS unidentified_sequences,
  count(*)::numeric / (select count(distinct read_sequence) from "tmpTest8_qs20") AS proportion_unid_reads
from test8_qs20_unaligned_sequences ;

----------------------------------[ Tritume ]----------------------------------

-- Sequence expression in aligned, not RT primers
create temp table aligned_not_rdna AS(
  select "tmpTest8_qs20".read_sequence, count("tmpTest8_qs20".read_sequence) AS sequence_count
  from "tmpTest8_qs20" 

    left join "tmpTest8_qs20_rDNA" on
      "tmpTest8_qs20".read_name = "tmpTest8_qs20_rDNA".read_name

    inner join "tmpSoapTest8_qs20_r1" on                 -- Retain aligned reads
      "tmpTest8_qs20".read_name = "tmpSoapTest8_qs20_r1".read_name

  where "tmpTest8_qs20_rDNA" is null AND                 -- Remove rDNA
    "tmpTest8_qs20".read_sequence not like '%AAGGTCTAT%' -- Remove RT primers (if any aligned)
  group by "tmpTest8_qs20".read_sequence
  order by sequence_count desc
  );
select * from aligned_not_rdna 


select "tmpTest8_qs20".read_sequence, count("tmpTest8_qs20".read_sequence) AS sequence_count
from "tmpTest8_qs20"
group by read_sequence
order by sequence_count desc
limit 100;

select * from aligned_not_rdna limit 100;

select * from test8_qs20_unaligned_sequences order by read_count desc limit 100;

select sum(read_count) from test8_qs20_unaligned_sequences;

select sum(read_count)
from test8_qs20_unaligned_sequences
where read_sequence like upper('%tggcctgctttgaacactctaatt%');


select sum(read_count) from test8_qs20_unaligned_sequences;

select 1324065 + 160860 + 421618 ;