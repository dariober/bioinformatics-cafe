/*
  Efficiency at reaching 5' end
*/

-- LPS library assembled denovo
drop table denovo_lps;
create temp table denovo_lps AS (
  select cufflinks_transcript_gtf.transcript_id,
         cufflinks_transcript_gtf.start,
         cufflinks_transcript_gtf.end,
         cufflinks_transcript_gtf.rname,
         cufflinks_transcript_gtf.strand,
         cufflinks_transcript_gtf.exon_number,
         cufflinks_transcript_gtf.fpkm,
         cuffcompare_tracking.ref_transcript_id
  from cuffcompare_tracking inner join cufflinks_transcript_gtf on
      cuffcompare_tracking.cufflinks_transcript_id = cufflinks_transcript_gtf.transcript_id
  where cuffcompare_tracking.ref_transcript_id is not null and 
        cuffcompare_tracking.class_code like '=' AND
        cufflinks_transcript_gtf.strand = '+' and
        cufflinks_transcript_gtf.exon_number like '1' and
        cufflinks_transcript_gtf.source like '20100317_RNAseq_LPS_noGTF'
  );

copy (
select denovo_lps.*, 
       ref.f_start,
       ref.f_end,
       (denovo_lps.start - ref.f_start) AS denovo_rel_position -- negative means denovo starts earlier
from denovo_lps inner join 
     (
      select *, gtf_attribute(sus_scrofa_sscrofa9_56_gtf.attributes, 'transcript_id') AS ref_transcript
      from sus_scrofa_sscrofa9_56_gtf
      where gtf_attribute(sus_scrofa_sscrofa9_56_gtf.attributes, 'exon_number') like '1' and
      strand like '+' and
      feature like 'exon'
      ) AS ref on
     denovo_lps.ref_transcript_id = ref.ref_transcript
     order by denovo_rel_position
) to 'C:/Tritume/denovo_tss.txt' with csv header delimiter E'\t';

select * from denovo_lps limit 210;
select ENSSSCT00000008324
-- 

-------------------------------[ Tritume ]---------------------------

select * from cuffcompare_tracking  limit 100;

select * from cufflinks_transcript_gtf limit 10