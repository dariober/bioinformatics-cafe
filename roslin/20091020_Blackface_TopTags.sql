drop table if exists sage_rs;
create temp table sage_rs AS(
  select "SAGEID", "Dataset", "TagChrPos", "E"
  from "outSAGEtop" 
  where "Dataset" like 'RvsS' and "E"<0.01
  order by "E"
  );

drop table if exists sage_rs2;
create temp table sage_rs2 AS(
  select sage_rs.*, "Read", "ReadCount", "SolexaRuns"."RunID", "SolexaRuns"."LibraryID", "TreatmentGroup"
  from sage_rs inner join "SolexaTags" on sage_rs."TagChrPos" = "SolexaTags"."TagChrPos" inner join
  "SolexaRuns" on "SolexaRuns"."RunID" = "SolexaTags"."RunID" inner join
  "qrySolexaLambs" on "qrySolexaLambs"."LibraryID" = "SolexaRuns"."LibraryID"
  where "InUse" = True
  order by "E"
  );

drop table if exists sage_rs3;
create table sage_rs3 AS(
  select "SAGEID", "Dataset", "TagChrPos", "E", max("ReadCount") AS maxcount
  from sage_rs2 
  group by "SAGEID", "Dataset", "TagChrPos", "E"
  );

select distinct sage_rs2."SAGEID", sage_rs2."Dataset", sage_rs2."TagChrPos", sage_rs2."E", sage_rs2."Read"
from sage_rs2 inner join sage_rs3 on 
  sage_rs2."SAGEID" = sage_rs3."SAGEID" and 
  sage_rs2."ReadCount" = sage_rs3.maxcount
  order by "E";

-- Expression levels by TagChrPos
select "TreatmentGroup", "TagChrPos", "E", avg("ReadCount"), stddev("ReadCount")
from sage_rs2
where "TreatmentGroup" not like 'C'
group by "TreatmentGroup", "TagChrPos", "E"
order by "E", "TreatmentGroup";

select * from sage_rs2;
