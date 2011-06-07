select read_sam($$ file:'F:/data/20101206_CAGE/20101207_CAGE_vs_sscrofa9_human_rDNA.sam', dest_table:'sam_cage_20101207'  $$)
alter table sam_cage_20101207 set tablespace hdd_free_agent;
comment on table sam_cage_20101207 is 'SAM file produced by pig CAGE tags template "CAGE 05/08/2010" see labbook 06/12/2010'

-- No of mismatches:
select sam_get_tagvalue(tags, 'NM') AS no_mismatches, count(*) AS count from sam_cage_20101207 group by sam_get_tagvalue(tags, 'NM');

select qname, sam_get_tagvalue(tags, 'NM'), tags from sam_cage_20101207 limit 10;

select * from sam_cage_20101207 limit 10;

select rname,
       count(distinct qname) AS n_reads, 
       (count(distinct qname)::numeric/(select count(distinct qname) from sam_cage_20101207))*100 AS freq 
from sam_cage_20101207 group by rname
order by freq desc;

select qname, count(qname) 
from sam_cage_20101207 
group by qname 
having count(qname) > 8;

select * from sam_cage_20101207 where qname like 'EBRI093151_0001:8:48:865:1743#0/1'


select count(distinct qname) from sam_cage_20101207 
where seq like '%AGACAGCAG%' or 
      seq like dna_revcomp('%AGACAGCAG%'); -- 40770 seqs

select count(distinct qname) from sam_cage_20101207 
where seq like '%AGACAG%' or 
      seq like dna_revcomp('%AGACAG%');

select 40770::numeric / 3899341;
select * from sam_cage_20101207 limit 36;

select qual, CASE WHEN flag & 16 = 16 THEN phred(reverse_string(qual), 'phred')
            ELSE phred(qual, 'phred') END
from sam_cage_20101207 
limit 36
where flag = 16

limit 10;
limit 1000; 

seq like '%AGACAGACG%' or 

select dna_revcomp('%AGACAGCAG%');

"HashAggregate  (cost=87722.06..87722.33 rows=22 width=3) (actual time=20122.801..20122.843 rows=22 loops=1)"
"  ->  Seq Scan on sam_cage_20101207  (cost=0.00..62143.27 rows=5115757 width=3) (actual time=0.008..9683.302 rows=5115757 loops=1)"
"Total runtime: 20122.942 ms"


drop index ind_sam_rname;
create index ind_sam_rname on sam_cage_20101207 (rname) tablespace hdd_free_agent; -- 83385 ms.
vacuum analyze sam_cage_20101207; -- 40919 ms.

create index ind_sam_qname on sam_cage_20101207 (qname) tablespace hdd_free_agent; -- 94235 ms.
vacuum analyze sam_cage_20101207; -- 2656 ms.