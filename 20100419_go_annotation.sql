/*
*  Generate a list of go terms for each ensembl transcript from S. scrofa
*  File C:/Downloads/mart_export.txt produced by biomart setting as attributes 'enesembl transcript id', and the three 'go term accession'
*  
*  This file used for input for GO annotation of RNAseq using topGO
*/

------------------[ Annotation datafile mapping transcripts => GO ]------------

select read_table($$ file:'C:/Downloads/mart_export.txt', dest_table:'go_annotation_ensembl_sscrofa57', header:True, limit:-1, overwrite:True $$)
select names($$ t:'go_annotation_ensembl_sscrofa57', rename:True, apply:['.lower()', '.replace(" ", "_")', 're.sub("\(|\)", "", x)', "re.sub('^go_term_accession$', 'go_term_accession_mf', x)"], echo:False  $$)
select names($$ t:'go_annotation_ensembl_sscrofa57', quote:'' $$)

create temp table go_annotation AS(
  select ensembl_transcript_id, array_to_string(array_agg(go_list), ', ') AS go_list
  from (select distinct ensembl_transcript_id, go_term_accession_bp AS go_list from go_annotation_ensembl_sscrofa57
        union select distinct ensembl_transcript_id, go_term_accession_cc from go_annotation_ensembl_sscrofa57
        union select distinct ensembl_transcript_id, go_term_accession_mf from go_annotation_ensembl_sscrofa57) AS t1 
  where go_list is not null
  group by ensembl_transcript_id
  ); -- 14740 annotated transcripts

-- Move this to U:/Documents/GO_annotation
copy go_annotation to 'C:/Tritume/go_annotation_ensembl_sscrofa57.txt' with csv delimiter E'\t';

----------------------------------[ Gene universe ]----------------------------
select names($$ t:'tmp_edger_transcript'  $$); "transcript_id", "logConc", "logFC", "PValue", "FDR"

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

-- Gene universe
copy (
  select tmp_edger_transcript.transcript_id, "logFC", "FDR", "PValue", fpkm_lps, fpkm_ctrl
  from tmp_edger_transcript 
    -- Restrict to transctipts having GO annotation
    inner join go_annotation on tmp_edger_transcript.transcript_id = ensembl_transcript_id
    -- Add FPKM
    inner join rnaseq_expr on rnaseq_expr.transcript_id = tmp_edger_transcript.transcript_id
  -- Remove transcripts not expressed in both libs
  where fpkm_lps is not null and fpkm_ctrl is not null
  order by "PValue"
  ) to 'C:/Tritume/rnaseq_gene_universe.txt' with csv header delimiter E'\t'; -- Move to U:/Documents/GO_annotation

-- Gene universe alterenative selection. Transcripts where both libraries combined make >n FPKM
copy (
  select tmp_edger_transcript.transcript_id, "logFC", "FDR", "PValue", fpkm_lps, fpkm_ctrl, (fpkm_lps + fpkm_ctrl) AS combined
  from tmp_edger_transcript 
    -- Restrict to transctipts having GO annotation
    inner join go_annotation on tmp_edger_transcript.transcript_id = ensembl_transcript_id
    -- Add FPKM
    inner join rnaseq_expr on rnaseq_expr.transcript_id = tmp_edger_transcript.transcript_id
  -- Remove transcripts where the combined expression is less than x
  where fpkm_lps + fpkm_ctrl > 5
  order by combined
  ) to 'C:/Tritume/rnaseq_gene_universe_fpkm_5.txt' with csv header delimiter E'\t';

-- Move to GO_annotation. Cygwin shell:
-- mv C:/Tritume/rnaseq_gene_universe_fpkm_5.txt U:/Documents/GO_annotation/

-------------------------------[ Tritume ]-------------------------------------
select * from tmp_edger_transcript_annotated inner join tmp_edger_transcript on tmp_edger_transcript.transcript_id = tmp_edger_transcript_annotated.ensembl_transcript_id 
where "go_term_name_(bp)" like 'negative regulation of organismal physiological process'



select go_term_name_bp from ensembl_annotation_sscrofa9_57 where  like 'GO:0051241';

select tmp_edger_transcript.*, count_pattern(go_list, 'GO:') AS no_go_terms from go_annotation inner join tmp_edger_transcript on transcript_id = ensembl_transcript_id
where go_list  like '%GO:0070667%'
order by no_go_terms desc

select distinct tmp_edger_transcript.* from go_annotation_ensembl_sscrofa57 inner join tmp_edger_transcript on transcript_id = ensembl_transcript_id
where go_term_accession_bp like 'GO:0006955'

select distinct tmp_edger_transcript.* from tmp_edger_transcript order by transcript_id


go_annotation_ensembl_sscrofa57



limit 10

select * from go_annotation_ensembl_sscrofa57 limit 10;

"ensembl_transcript_id, go_term_accession_bp, go_term_accession_cc, go_term_accession_mf"