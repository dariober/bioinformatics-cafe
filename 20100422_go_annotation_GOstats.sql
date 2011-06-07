------------------[ Annotation datafile mapping transcripts => GO ]------------

-- C:/Downloads/mart_export.txt downloaded from Biomart.
select read_table($$ file:'C:/Downloads/mart_export.txt', dest_table:'go_annotation_ensembl_sscrofa57', header:True, limit:-1, overwrite:True $$)
select names($$ t:'go_annotation_ensembl_sscrofa57', rename:True, apply:['.lower()', '.replace(" ", "_")', 're.sub("\(|\)", "", x)', "re.sub('^go_term_accession$', 'go_term_accession_mf', x)"], echo:False  $$)
select names($$ t:'go_annotation_ensembl_sscrofa57', quote:'' $$)

-- ensembl_transcript_id, go_term_accession_bp, go_term_evidence_code_bp, go_term_accession_cc, go_term_evidence_code_cc, go_term_accession_mf, go_term_evidence_code_mf

DROP TABLE go_annotation_long;
CREATE TEMP TABLE go_annotation_long AS (
  -- GO annotation format as per GOstats
  select distinct go_term_accession_bp AS go_id, go_term_evidence_code_bp AS evidence, ensembl_transcript_id AS gene_id 
    from go_annotation_ensembl_sscrofa57 
    where go_term_accession_bp is not null and go_term_evidence_code_bp is not null
  union select distinct go_term_accession_cc AS go_id, go_term_evidence_code_cc AS evidence, ensembl_transcript_id AS gene_id 
    from go_annotation_ensembl_sscrofa57
        where go_term_accession_cc is not null and go_term_evidence_code_cc is not null
  union select distinct go_term_accession_mf AS go_id, go_term_evidence_code_mf AS evidence, ensembl_transcript_id AS gene_id 
    from go_annotation_ensembl_sscrofa57
        where go_term_accession_mf is not null and go_term_evidence_code_mf is not null
  );
select distinct * from go_annotation_long limit 100;

copy go_annotation_long to 'C:/Tritume/go_annotation_ensembl_sscrofa57_GOstats.txt' with csv header delimiter E'\t';

-- Move this to U:/Documents/GO_annotation/GOstats. Cygwin:
-- mv C:/Tritume/go_annotation_ensembl_sscrofa57_GOstats.txt U:/Documents/GO_annotation/GOstats

----------------------------------[ Gene universe ]----------------------------

select names($$ t:'tmp_edger_transcript'  $$); -- "transcript_id", "logConc", "logFC", "PValue", "FDR"

-- FPKM for LPS and CTRL libs
select cross_tab($$ 
  select transcript_id, source, fpkm
  from cufflinks_transcript_gtf 
  where source in ('20100317_RNAseq_CTRL', '20100317_RNAseq_LPS') and feature like 'transcript' $$,
  'tmp_cufflinks_transcript_gtf_ct');

-- Rename columns to better names
select names($$ t:'tmp_cufflinks_transcript_gtf_ct', which:[2,3], rename:True, apply: ["re.sub('.*CTRL.*', 'fpkm_ctrl', x)", "re.sub('.*LPS.*', 'fpkm_lps', x)"], echo:False  $$)

-- Make temp table and delete permanent one (not to litter the DB)
create temp table rnaseq_expr AS(select * from tmp_cufflinks_transcript_gtf_ct);
drop table tmp_cufflinks_transcript_gtf_ct;

-- Gene universe alterenative selection. Transcripts where both libraries combined make >n FPKM
copy (
  select distinct tmp_edger_transcript.transcript_id, ensembl_annotation_sscrofa9_57.associated_gene_name, "logFC", "FDR", "PValue", fpkm_lps, fpkm_ctrl, (fpkm_lps + fpkm_ctrl) AS combined
  from tmp_edger_transcript 
    -- Restrict to transctipts having GO annotation
    inner join go_annotation_long on tmp_edger_transcript.transcript_id = go_annotation_long.gene_id
    -- Add FPKM
    inner join rnaseq_expr on rnaseq_expr.transcript_id = tmp_edger_transcript.transcript_id
    -- Add gene annotation if available
    left join ensembl_annotation_sscrofa9_57 on ensembl_annotation_sscrofa9_57.ensembl_transcript_id = tmp_edger_transcript.transcript_id
  -- Remove transcripts where the combined expression is less than x
  where fpkm_lps + fpkm_ctrl > 5
  order by combined
  ) to 'C:/Tritume/rnaseq_gene_universe_fpkm_5.txt' with csv header delimiter E'\t';

-- Move to GO_annotation. Cygwin shell:
-- mv C:/Tritume/rnaseq_gene_universe_fpkm_5.txt U:/Documents/GO_annotation/
