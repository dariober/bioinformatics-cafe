/*
  Output from cuffcompare, formatted and parsed
*/
select read_table($$ file:'C:/Tritume/20100505_RNAseq_cuffcompare/RNAseq.tracking2', dest_table:'cuffcompare_tracking', dataset_id:'20100317_RNAseq', header:True, overwrite:False $$)
comment on table cuffcompare_tracking is 'Output from cuffcompare, file: <basename>.tracking and parsed by 20100505_cuffcompare_tracking.py'
select index('cuffcompare_tracking', 'dataset_id, cufflinks_transcript_id', 'p');

select * from cuffcompare_tracking limit 100;
select * from cufflinks_transcript_gtf where source in ('20100317_RNAseq_LPS_noGTF', '20100317_RNAseq_CTRL_noGTF') limit 10;

-- How many denovo transcripts have 1, 2, 3... n exons?
-- Number of exons per transcript, with transcript FPKM and class_code from cuffcompare
-- This query takes a while
drop table exons;
create temp table exons AS (
  select max(exon_number::int) AS no_exons,
         source, 
         transcript_id, 
         max(fpkm) AS fpkm, 
         class_code,
         cuffcompare_transcript_id,
         cuffcompare_locus_id,
         ref_transcript_id
  from cufflinks_transcript_gtf 
  -- Add class_code
  left join cuffcompare_tracking on 
      cufflinks_transcript_gtf.transcript_id = cuffcompare_tracking.cufflinks_transcript_id
  -- 
  where cuffcompare_tracking.dataset_id like '20100317_RNAseq' or cuffcompare_tracking.dataset_id is null
  group by source, transcript_id, class_code,cuffcompare_transcript_id, cuffcompare_locus_id, ref_transcript_id
  );

-- This table sent to R 20100514_cuffcompare_denovo.R
copy exons to 'C:/Tritume/exons.txt' with csv header delimiter E'\t';
select * from exons limit 10;
select count(*) from exons;
select * from exons where transcript_id like 'ENSSSCT00000008324';
select * from exons where transcript_id like 'CTRL.10.0';

------------------------------[ Select de novo transcripts ]-------------------
-- See labbook 11/5/10 for criteria
-- These data used for edgeR
drop table tmp_denovo;
create table tmp_denovo AS(
select exons.cuffcompare_transcript_id, source, exons.fpkm from
    (
    select count(cuffcompare_transcript_id) as no_libs, sum(fpkm) as tot_fpkm, max(no_exons) as no_exons, class_code, cuffcompare_transcript_id
    from exons 
    where source in ('20100317_RNAseq_LPS_noGTF', '20100317_RNAseq_CTRL_noGTF')
    group by cuffcompare_transcript_id, class_code
--    having (sum(fpkm) > 24 or max(no_exons) > 1 or class_code not like 'u') and class_code not like 'i'
    -- Accept transcripts found in both libraries
    having count(cuffcompare_transcript_id) > 1
    ) AS t1 inner join exons on exons.cuffcompare_transcript_id = t1.cuffcompare_transcript_id
where source in ('20100317_RNAseq_LPS_noGTF', '20100317_RNAseq_CTRL_noGTF')
order by source, transcript_id
);
select * from tmp_denovo limit 100;
drop table tmp_denovo_edger;
select cross_tab('select * from tmp_denovo', 'tmp_denovo_edger');
select names($$ t:'tmp_denovo_edger', which:(2,3), rename:('fpkm_ctrl', 'fpkm_lps'), echo:False  $$)

select * from tmp_denovo_edger limit 1000;
select count(*) from tmp_denovo_edger where (fpkm_lps, fpkm_ctrl) is not null;

select * from exons limit 20;

-----------------------[ Analyze differential expression ]---------------------

-- This dataset produced by 20100514_edger_denovo.R
create temp table edger_denovo AS(
    select * from edger_toptags 
    where dataset_id like '20100514_CTRLvsLPS_denovo' and "FDR" < 0.05
    );
select * from edger_denovo;

-- 
select distinct cuffcompare_tracking.*
from cuffcompare_tracking inner join edger_denovo on
    edger_denovo.transcript_id = cuffcompare_tracking.cuffcompare_transcript_id
where class_code in ('', 'j');


-- Annotate differential transcripts
select class_code, count(distinct transcript_id)
from cuffcompare_tracking inner join edger_denovo on
    edger_denovo.transcript_id = cuffcompare_tracking.cuffcompare_transcript_id
group by class_code

select distinct ref_transcript_id from 
cuffcompare_tracking inner join edger_denovo on
edger_denovo.transcript_id = cuffcompare_tracking.cuffcompare_transcript_id
where class_code like 'j';

"TCONS_00051450"
--------------[ Table of average expression and single exons ]-----------------
-- The following two queries make the table on 12/5/10 after reformatting in Excel

select * from pg_stat_activity;
select pg_cancel_backend(2380);

-- Number of sngle-exon transcripts with their count and stdev, sterr
select count(distinct transcript_id), no_exons, source, avg(fpkm), stddev(fpkm)/sqrt(count(transcript_id)) AS sterr
from exons 
where no_exons = 1
group by no_exons, source
order by source;

select avg(fpkm) AS avg_rpkm, stddev(fpkm) AS stdev_rpkm, 
    count(transcript_id) AS n_transcripts, stddev(fpkm)/sqrt(count(transcript_id)) AS sterr_fpkm, source 
from exons 
group by source;

-------------------[ Tritume ]-------------------------------------------------
select distinct *   -- count(distinct tmp_denovo_edger.cuffcompare_transcript_id), class_code 
from tmp_denovo_edger inner join cuffcompare_tracking on 
tmp_denovo_edger.cuffcompare_transcript_id = cuffcompare_tracking.cuffcompare_transcript_id
inner join edger_toptags on tmp_denovo_edger.cuffcompare_transcript_id = edger_toptags.transcript_id
where edger_toptags.dataset_id like '20100514_CTRLvsLPS_denovo' and "FDR" < 0.01;
group by class_code;

select 

select distinct edger_toptags.*, cuffcompare_tracking.ref_gene_id, cuffcompare_tracking.ref_transcript_id, cuffcompare_tracking.class_code, abs("logFC")
from edger_toptags inner join cuffcompare_tracking on
    edger_toptags.transcript_id = cuffcompare_tracking.cuffcompare_transcript_id
where edger_toptags.dataset_id like '20100514_CTRLvsLPS_denovo' and "FDR" < 0.01 and class_code like 'j'
order by abs("logFC") desc;

select * from edger_toptags where edger_toptags.dataset_id like '20100514_CTRLvsLPS_denovo' limit 10;
select names('cuffcompare_tracking'); -- "dataset_id", "cuffcompare_transcript_id", "cuffcompare_locus_id", "ref_gene_id", "ref_transcript_id", "class_code", "cufflinks_gene_id", "cufflinks_transcript_id"


select * 
from cuffcompare_tracking inner join cufflinks_transcript_gtf 
 on cufflinks_transcript_gtf.transcript_id = cuffcompare_tracking.cufflinks_transcript_id
where ref_gene_id like 'PARP14' and class_code like 'j';

select * from cufflinks_transcript_gtf where transcript_id like 'ENSSSCT00000012988';

where class_code like 'j';

select * from cuffcompare_tracking where cuffcompare_transcript_id like 'TCONS_00000686';

select * from cufflinks_transcript_gtf where transcript_id like 'LPS.155595.0'
union select * from cufflinks_transcript_gtf where transcript_id like 'ENSSSCT00000001055';

select distinct * from cuffcompare_tracking  
inner join cufflinks_transcript_gtf on cuffcompare_tracking.cufflinks_transcript_id = cufflinks_transcript_gtf.transcript_id
where 
    -- feature like 'transcript' and 
    class_code like 'u'
order by rname, "start", "end"
limit 1000;

select count(*) from tmp_denovo_edger;

select * from cufflinks_transcript_gtf where transcript_id like 'ENSSSCT00000000020';

drop function array_median(double precision[]);
drop aggregate median(double precision);


create or replace function median(anyarray) 
returns double precision as $$
  select ($1[array_upper($1,1)/2+1]::double precision + $1[(array_upper($1,1)+1) / 2]::double precision) / 2.0; 
$$ language sql immutable strict;

    select source,
        'CTRL ref' AS dataset,
        count(distinct transcript_id) AS n,
        median(array(select fpkm 
                     from exons 
                     where fpkm is not null and 
                         source like '20100317_RNAseq_CTRL' and
                         no_exons = 1
                     order by fpkm)) AS median,
        avg(fpkm), stddev(fpkm)/sqrt(count(transcript_id)) AS sterr
    from exons where source like '20100317_RNAseq_CTRL' and no_exons = 1 group by source

    UNION select source,
        'CTRL denovo' AS dataset, 
        count(distinct transcript_id) AS n,
        median(array(select fpkm 
                     from exons 
                     where fpkm is not null and 
                         source like '20100317_RNAseq_CTRL_noGTF'
                     order by fpkm)) AS median,
        avg(fpkm), stddev(fpkm)/sqrt(count(transcript_id)) AS sterr
    from exons where source like '20100317_RNAseq_CTRL_noGTF' group by source

select * from exons limit 10;
    
 "20100317_RNAseq_CTRL"
"20100317_RNAseq_CTRL_noGTF"
"20100317_RNAseq_LPS"
"20100317_RNAseq_LPS_noGTF"


select distinct source from exons;

select count(distinct transcript_id), no_exons, source, median(array[fpkm])
from exons 
where no_exons = 1
group by no_exons, source
order by source;

select median(array(select fpkm from exons));


drop table denovo;
create temp table denovo AS (
    -- Assign FPKM to the denovo transcripts
    select cuffcompare_tracking.*, fpkm 
    from cuffcompare_tracking inner join cufflinks_transcript_gtf on
        cuffcompare_tracking.cufflinks_transcript_id = cufflinks_transcript_gtf.transcript_id
    where cufflinks_transcript_gtf.source in ('20100317_RNAseq_LPS_noGTF', '20100317_RNAseq_CTRL_noGTF') AND
        cuffcompare_tracking.dataset_id = '20100317_RNAseq' AND
        cufflinks_transcript_gtf.feature like 'transcript'
    );

-- Select genes 
select count(distinct cuffcompare_transcript_id), class_code from denovo group by class_code;
select count(*) 
select * from denovo limit 100;


limit 20;

select class_code, count(transfag_id) AS no_transfags from cuffcompare_tracking group by class_code order by no_transfags;

select * from cuffcompare_tracking where class_code like 'o' order by ref_gene_id;
select * from cuffcompare_tracking where class_code like '=' order by ref_gene_id limit 100;
select * from cuffcompare_tracking where locus_id like 'XLOC_029532';
select * from cuffcompare_tracking where ref_transcript_id like 'ENSSSCT00000000756';


select * from cufflinks_transcript_gtf where transcript_id in ('LPS.115774.0', 'CTRL.122412.0', 'ENSSSCT00000008324');
select * from cufflinks_transcript_gtf where transcript_id in ('LPS.152114.0', 'CTRL.160181.0', 'ENSSSCT00000000756');
select * from cufflinks_transcript_gtf where transcript_id in ('LPS.152798.0', 'CTRL.160868.0', 'ENSSSCT00000000829');

select * from cufflinks_transcript_gtf 
inner join cuffcompare_tracking on cufflinks_transcript_gtf.transcript_id = cuffcompare_tracking.transcript_id_1
where fpkm > 10 and feature like 'transcript' and source in ('20100317_RNAseq_LPS_noGTF', '20100317_RNAseq_CTRL_noGTF');


drop table denovo;
create temp table denovo AS(
    select distinct cuffcompare_tracking.ref_transcript_id, cufflinks_transcript_gtf.*
    from cufflinks_transcript_gtf 
    inner join cuffcompare_tracking on 
        cufflinks_transcript_gtf.transcript_id = cuffcompare_tracking.transcript_id_1
    where class_code like '='
    UNION 
    select distinct cuffcompare_tracking.ref_transcript_id, cufflinks_transcript_gtf.*
    from cufflinks_transcript_gtf 
    inner join cuffcompare_tracking on 
        cufflinks_transcript_gtf.transcript_id = cuffcompare_tracking.transcript_id_2
    where class_code like '='
    );
select * from denovo limit 10;

drop table ref;
create temp table ref AS(
    select distinct cufflinks_transcript_gtf.*
    from cufflinks_transcript_gtf
    inner join (select distinct ref_transcript_id from denovo) AS t1 on t1.ref_transcript_id = cufflinks_transcript_gtf.transcript_id
    );
select * from ref limit 10;

select distinct denovo.exon_number::numeric, denovo.transcript_id AS denovo_id, denovo.fpkm, denovo.rname, denovo.start AS denovo_start, denovo.end AS denovo_end, denovo.strand AS denovo_strand,
                ref.transcript_id AS ref_id, ref.rname, ref.start AS ref_start, ref.end AS ref_end, ref.strand AS ref_strand,
                (ref.end - ref.start) AS ref_length,
                (denovo.end - denovo.start) AS denovo_length,
                (denovo.start - ref.start) AS start_diff,  -- Negative means denovo begins earlier
                (denovo.end - ref.end) AS end_diff         -- Positive means denovo ends later
from denovo inner join ref on 
    denovo.ref_transcript_id = ref.transcript_id AND
    denovo.exon_number = ref.exon_number
where ref.strand like '+' and denovo.exon_number like '1'
order by denovo.rname, denovo.start;

select * from cufflinks_transcript_gtf where transcript_id like 'ENSSSCT00000014300'

from cufflinks_transcript_gtf
inner join (select * from cuffcompare_tracking where class_code like '=') AS denovo
      on cufflinks_transcript_gtf.transcript_id = denovo.transcript_id_1     
inner join (select * from cuffcompare_tracking where class_code like '=') AS ref
      on denovo.transcript_id = ref.ref_transcript_id
limit 20;



select cufflinks_transcript_gtf.transcript_id, cufflinks_transcript_gtf.rname, cufflinks_transcript_gtf.start, cufflinks_transcript_gtf.end, ref.ref_transcript_id
from cufflinks_transcript_gtf
    inner join cuffcompare_tracking on cufflinks_transcript_gtf.transcript_id = cuffcompare_tracking.transcript_id_1
    inner join  (select * from cuffcompare_tracking) AS ref on cufflinks_transcript_gtf.transcript_id = ref.ref_transcript_id
where exon_number like '1' and ref.class_code like '=' 
limit 10;

select * from cufflinks_transcript_gtf where transcript_id like 'LPS.144792.0';

122238


select cuffcompare_transcript_id, count(cuffcompare_transcript_id)
from cuffcompare_tracking 
where class_code like '='
group by cuffcompare_transcript_id
limit 10;