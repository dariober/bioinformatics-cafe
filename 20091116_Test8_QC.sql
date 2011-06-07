/*
  Detecting primer and sequence contamination in library CAGE Test8
*/

---------------------------[ All test8 reads ]---------------------------------

select count(read_sequence), count(read_sequence)::numeric / (select count(*) from fastq where fastq_name = 'test8')
from fastq 
where read_sequence like '%AAGGTCTAT%' 
  and fastq_name = 'test8';

select count(read_sequence), count(read_sequence)::numeric / (select count(*) from fastq where fastq_name = 'test8')
from fastq 
where read_sequence like '%CCACCG%' or
      read_sequence like '%CTACAG%' or
      read_sequence like '%CTGTAG%' or
      read_sequence like '%CGGTGG%' or
      read_sequence like  '%TCGTAT%' or
      read_sequence like  '%TGCTTG%' or
      read_sequence like  '%CAAGCA%' or
      read_sequence like  '%ATACGA%'
  and fastq_name = 'test8';

------------------------[ High quality reads ]---------------------------------

select count(read_sequence), count(read_sequence)::numeric / (select count(*) from "tmpTest8_qs20")
from "tmpTest8_qs20"
where read_sequence like '%CTGTAG%'

select count(read_sequence), count(read_sequence)::numeric / (select count(*) from "tmpTest8_qs20")
from "tmpTest8_qs20"
where read_sequence like '%CCACCG%' or
      read_sequence like '%CTACAG%' or
      read_sequence like '%CTGTAG%' or
      read_sequence like '%CGGTGG%' or
      read_sequence like  '%TCGTAT%' or
      read_sequence like  '%TGCTTG%' or
      read_sequence like  '%CAAGCA%' or
      read_sequence like  '%ATACGA%';

----------------------------[ Tritume ]----------------------------------------

select * from fastq where read_sequence like '%AAGGTCTAT%' limit 100;

select count("TSS"), count("TSS")::numeric / (select count(*) from "tmpSeq")
from "tmpSeq"
where "TSS" like '%AAGGTCTAT%'

select count(read_sequence), count(read_sequence)::numeric / (select count(*) from "tmpTest8_qs20")
from "tmpTest8_qs20"
where read_sequence like '%AAGGTCTAT%'


