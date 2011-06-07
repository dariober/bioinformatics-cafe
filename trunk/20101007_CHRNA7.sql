select read_table($$ file:'F:/data/20101006_CHRNA7/tblastn_CHRNA7_human_vs_ss9.56.txt', dest_table: 'tblastn_chrna7', header:True $$);
select * from tblastn_chrna7;

select subject_id || ':' || s_start || ':' || s_end from tblastn_chrna7;

drop table blast_to_gtf;
create temp table blast_to_gtf AS (
  -- Blasted regions overlapping an annotated transcript (i.e. feature in GTF)
  select distinct *, 
      gtf_attribute(attributes, 'transcript_id') AS transcript_id, 
      gtf_attribute(attributes, 'gene_id') AS gene_id, 
      gtf_attribute(attributes, 'gene_name') AS gene_name, 
      gtf_attribute(attributes, 'exon_number') AS exon_number
  from tblastn_chrna7 inner join sus_scrofa_sscrofa9_56_gtf on subject_id::varchar = rname
  where (s_start between f_start and f_end) OR
        (s_end between f_start and f_end) OR
        (s_start <= f_start and s_end >= f_end) OR
        (s_end <= f_start and s_start >= f_end)
  order by subject_id, s_start
  );
copy (select distinct query_id, subject_id, identity, alignment_length, mismatches, gap_openings, q_start, q_end, s_start, s_end, e_value, bit_score, transcript_id,gene_id, gene_name
      from blast_to_gtf order by bit_score) 
  to 'F:/data/20101006_CHRNA7/blast_to_gtf.txt' with csv header delimiter E'\t';
/*
These are the cases to cover by the WHERE clauses

s_start |-----------------| s_end                                             (s_end between f_start and f_end) 
                                        s_start |-----------------| s_end     (s_start between f_start and f_end)
s_start |---------------------------------------------------------| s_end     (s_start <= f_start and s_end >= f_end)
s_end   |---------------------------------------------------------| s_start   (s_end <= f_start and s_start >= f_end)
           f_start |-----------------------------------| f_end
*/

-- Blasted regions without annotation
copy
  (select tblastn_chrna7.* from tblastn_chrna7 left join blast_to_gtf on
      tblastn_chrna7.subject_id = blast_to_gtf.subject_id AND
      tblastn_chrna7.s_start = blast_to_gtf.s_start AND
      tblastn_chrna7.s_end = blast_to_gtf.s_end
  where blast_to_gtf.attributes is null
  ) to 'F:/data/20101006_CHRNA7/no_annotation.txt' with csv header delimiter E'\t';



-- RNAseq expression for blasted regions/transcripts
copy (
  select distinct cufflinks_transcript_gtf.* from cufflinks_transcript_gtf, blast_to_gtf
  where cufflinks_transcript_gtf.transcript_id like blast_to_gtf.transcript_id
    and cufflinks_transcript_gtf.feature like 'transcript' and cufflinks_transcript_gtf.source like '20100602%'
  order by cufflinks_transcript_gtf.transcript_id, cufflinks_transcript_gtf.source
  ) to 'F:/data/20101006_CHRNA7/chrna7_expression.txt' with csv header delimiter E'\t';

-------------------------------[ Tritume ]-------------------------------------

select * from sus_scrofa_sscrofa9_56_gtf where attributes like '%ENSSSCT00000012894%' order by f_start;
select * from blast_to_gtf where attributes like '%ENSSSCT00000001958%' order by f_start;


select * from limma_toptable inner join porcine_na30_ensembl on limma_toptable.id = porcine_na30_ensembl.probe_set_id
where transcript_id in ('ENSSSCT00000009622', 'ENSSSCT00000009620');

s_start 53415793        s_end 53414777   subject
f_start 53408570  	f_end 53426370   feature