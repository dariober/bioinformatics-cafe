drop table multi_tags;
create temp table multi_tags AS(
  select "TagChrPos", "Read", sum("ReadCount") AS no_reads
  from "SolexaTags"
  group by "TagChrPos", "Read"
  order by "TagChrPos"
  );
select * from multi_tags limit 100;

drop table tag_summary;
create temp table tag_summary AS(
  select "TagChrPos", count(distinct "Read") AS no_distinct_reads, sum(no_reads) AS tot_count
  from multi_tags
  group by "TagChrPos" 
  order by "TagChrPos"
  );
select * from tag_summary limit 10;

drop table tag_freq;
create temp table tag_freq AS(
  select tag_summary.*, multi_tags."Read", multi_tags.no_reads/tot_count AS read_frequency, 
         dense_rank() OVER(PARTITION BY tag_summary."TagChrPos" ORDER BY multi_tags.no_reads DESC) AS read_rank
  from tag_summary inner join multi_tags on tag_summary."TagChrPos" = multi_tags."TagChrPos"
  order by "TagChrPos", read_frequency desc
  )
copy tag_freq to 'F:/data/20101023_blackface_DGE/tag_frequency.txt' with csv delimiter E'\t' header;
select * from tag_freq limit 100;

copy (
  select read_rank, 
       avg(read_frequency) AS avg_frequency, 
       stddev("read_frequency"), 
       count("TagChrPos") AS n
  from tag_freq
  where read_frequency < 1
  group by read_rank
  order by read_rank
  ) to 'F:/data/20101023_blackface_DGE/read_rank.txt' with csv header delimiter E'\t';

select * from tag_freq where read_rank = 63;

create temp table est_tags AS(
  select "SolexaTagsTGI"."TGI_ID", "TGIposition", sum("ReadCount") AS no_reads, count(distinct "ReadSeq") AS no_distinct_reads
    from "SolexaTagsTGI" inner join "SolexaRuns" on "SolexaRuns"."RunID" = "SolexaTagsTGI"."RunID"
    where "SolexaRuns"."InUse" = True
    group by "TGI_ID", "TGIposition"
    order by "TGI_ID", no_reads desc
    );
select * from est_tags limit 20;

copy (
  select "TGI_ID", count("TGIposition") AS no_positions, sum(no_reads) AS total_reads
  from est_tags
  group by "TGI_ID"
  order by no_positions desc
  ) TO 'F:/data/20101023_blackface_DGE/est_positions.txt' with csv header delimiter E'\t';

------------------------

copy (
  select out_solexa_cumcountall_tpm."LibraryID", sum("CumCount") AS no_tags, count(distinct "TagChrPos") AS no_distinct_positions, "TreatmentGroup"
  from out_solexa_cumcountall_tpm inner join "qrySolexaLambs" on out_solexa_cumcountall_tpm."LibraryID" = "qrySolexaLambs"."LibraryID"
  where "CumCount" > 0
  group by out_solexa_cumcountall_tpm."LibraryID", "TreatmentGroup"
  order by "TreatmentGroup", "LibraryID"
  ) TO 'F:/data/20101023_blackface_DGE/counts_per_lib.txt' with csv header delimiter E'\t';


create temp table est_reads AS(
  select "TGI_ID", count(distinct "ReadSeq") AS no_distinct_reads, sum("ReadCount") AS no_reads
  from "SolexaTagsTGI" inner join "SolexaRuns" on "SolexaRuns"."RunID" = "SolexaTagsTGI"."RunID"
  where "SolexaRuns"."InUse" = True
  group by "TGI_ID"
  order by "TGI_ID", no_reads desc
  ); -- 88319 rows

select * from est_reads limit 10;



select "SolexaTagsTGI".*, est_reads.* 
from "SolexaTagsTGI" inner join est_reads on "SolexaTagsTGI"."TGI_ID" = est_reads."TGI_ID"


