-- Calculate transcriptome size (i.e. length of all non-overlapping exons)

-- drop table tmp_cluster_exon;
select cluster_reads($$
  -- Cluster exons so not to double count bases spanned by more than one exon
  select rname, '+'::text, f_start, f_end, substring(attributes from position('transcript_id ' in attributes)+14 for 18) AS transcript_id, feature, score, strand, frame, attributes 
  from sus_scrofa_sscrofa9_56_gtf 
  where 
    substring(attributes from position('transcript_id ' in attributes)+14 for 18) 
      -- Include only exons belonging to these transcripts
      in (select distinct transcript_id from cufflinks_transcript_gtf where source like '20100317_RNAseq_CTRL')
      and feature like 'exon'
  $$, 
'tmp_cluster_exon');

select * from tmp_cluster_exon limit 20;

select 
  -- Sum the size of all clustered exons
  count(distinct transcript_id) AS no_transcripts, count(tag_cluster) AS non_overlap_exons, sum(exon_length) AS transcriptome_length
  from ( 
    select transcript_id, tag_cluster, min(f_start), max(f_end), max(f_end) - min(f_start) AS exon_length 
    from tmp_cluster_exon 
    group by tag_cluster, transcript_id) as t1; -- 16967188


-- Quick version for calculating the transcriptome size (ignore overlapping, incorrect but minimal difference)
select count(*) AS no_features, count(distinct transcript_id) AS no_transcripts, sum(exon_length) AS transcriptome_length
from (
  select *, substring(attributes from position('transcript_id ' in attributes)+14 for 18) AS transcript_id, f_end - f_start AS exon_length from sus_scrofa_sscrofa9_56_gtf 
   where substring(attributes from position('transcript_id ' in attributes)+14 for 18) 
      in (select distinct transcript_id from cufflinks_transcript_gtf where source like '20100317_RNAseq_CTRL')
      and feature like 'exon'
  ) AS t1
 limit 20; 
 
select count(tag_cluster) as no_clusters, tag_cluster from tmp_cluster_exon group by tag_cluster having count(tag_cluster) > 1 order by no_clusters desc limit 100;



-- where tag_cluster in ('2_+_55843036', '7_+_27718556', '7_+_27715820')
select * from tmp_cluster_exon where tag_cluster in ('2_+_55843036', '7_+_27718556', '7_+_27715820')



select distinct transcript_id from cufflinks_transcript_gtf where source like '20100317_RNAseq_CTRL';

select names($$ t:'sus_scrofa_sscrofa9_56_gtf', quote:'' $$)

"rname, source, feature, f_start, f_end, score, strand, frame, attributes"

select distinct * from tmp_edger_transcript inner join cufflinks_transcript_gtf on cufflinks_transcript_gtf.transcript_id = tmp_edger_transcript.transcript_id
where "PValue" < 0.01
  and feature like 'transcript';