--select read_table($$ file:'C:/Tritume/Porcine.na30.annot.csv', dest_table: 'porcine_na30_annot', sep:',', comment_char:'#', header:True, limit: -1 $$);
--select names($$ t:'porcine_na30_annot', apply: ['.lower()', ".replace(' ', '_')"], rename:True $$);

--select read_table($$ file:'C:/Tritume/mart_export.txt', dest_table: 'porcine_na30_ensembl', header:True, sep:'\t', limit: -1, limit:-1 $$)
--select names($$ t:'tmp_porcine_na30_ensembl', apply: ['.lower()', ".replace(' ', '_')"], rename:['probe_set_id', 'transcript_id'], echo:False $$);

--select read_table($$ file:'C:/Tritume/lwl_bmdm_lps_long.txt', dest_table:'affyarray_bmdm', sep:'\t', header:True, limit:-1, overwrite:True $$) -- C:/Tritume/pig_bmdm_filtered_probes.txt

--select read_table ($$ file: 'C:/Tritume/pig_bmdm_array_design.txt', dest_table:'affyarray_design', sep:'\t', header:True, limit:-1, overwrite:True $$)

------------------------------------[ Prepare table for limma ]---------------------------
/*
*  This part transforms the array data from long format to cross-tab and takes the mean
*  of each tissue culture replicate.
*  The output from the cross-tab table has column order different from the desing matrix
*  so the order has been changed manually.
*/

create table tmp_affyarray_bmdm_avg AS(
  select probe_set_id, pig || '_' || time_point || 'h'::text AS "Target", avg(intensity) as avg_intensity, pig, time_point
  from 
    affyarray_bmdm inner join affyarray_design on affyarray_bmdm.slide_id = affyarray_design.slide_id
  group by probe_set_id, "Target", pig, time_point
  order by pig, time_point
  )
-- Cross-tab from long format. NB column order changed manually
select cross_tab('select probe_set_id, "Target", avg_intensity::text from tmp_affyarray_bmdm_avg', 'tmp_affyarray_bmdm_avg_ct');
-- 'C:/Tritume/lwl_bmdm_lps_ct.txt' used for limma. R script 20100413_affyarray_RNAseq.R

/* This step probably no longer necessary
-- This is the file sent to limma after having changed the column order
copy tmp_affyarray_bmdm_avg_ct to 'C:/Tritume/lwl_bmdm_lps_avg_ct.txt' with csv delimiter E'\t' header;
*/

-- drop table tmp_affyarray_bmdm_avg;
-- drop table tmp_affyarray_bmdm_avg_ct;


--------------------------------[ Differential expression ]----------------------------

-- Read back the output from limma's topTable()
-- Comparing 0 and 7h
select read_table($$ overwrite:True, file:'C:/data/lwl_bmdm_lps_avg_ct_0h_7h.limma', dest_table:'limma_toptable', header:True, dataset_id:'avgrepl_contr_0h_7h', sep:'\t' $$)
select names($$ t:'limma_toptable', apply:['.lower()', ".replace('.', '_')"], rename: True, echo:False $$) -- ["dataset_id", "id", "logfc", "aveexpr", "t", "p_value", "adj_p_val", "b"]

select count(*) from limma_toptable where adj_p_val < 0.01;

-- Make cross-tab of RNAseq expression for CTRL and LPS libraries assembled (cufflinks) using GTF annotation
select cross_tab($$ 
  select transcript_id, source, fpkm 
  from cufflinks_transcript_gtf 
  where source in ('20100317_RNAseq_CTRL', '20100317_RNAseq_LPS') and 
    feature like 'transcript'$$, 
  'tmp_cufflinks_transcript_gtf_ct');
select names($$ t:'tmp_cufflinks_transcript_gtf_ct', which:[2,3], rename:['fpkm_ctrl', 'fpkm_lps'] $$);
select * from tmp_cufflinks_transcript_gtf_ct;

-- Make cross-tab of array data
-- drop table tmp_affyarray_bmdm_ct;
drop table tmp_affyarray_bmdm_ct;
select cross_tab($$
  select probe_set_id, slide_id, intensity
  from affyarray_bmdm $$,
  'tmp_affyarray_bmdm_ct')

-- Annotate toptable()
drop table affy_rnaseq;
create temp table affy_rnaseq AS(
select distinct
  u_transcripts.transcript_id,
  CASE WHEN porcine_na30_annot.gene_symbol LIKE '---' THEN '' ELSE porcine_na30_annot.gene_symbol END AS gene_symbol,
  -- Remove this field if importing to R since it causes problems to read.table()                                                                               
--  CASE WHEN porcine_na30_annot.gene_title LIKE '---' THEN '' ELSE porcine_na30_annot.gene_title END AS gene_title,  <--| Nota Bene   
  fpkm_ctrl, fpkm_lps,                                                                                             
  logfc AS affy_log2fc, adj_p_val AS affy_fdr,
  "logFC" AS rnaseq_log2fc, "FDR" AS rnaseq_fdr,
  tmp_affyarray_bmdm_ct.*
from 
  -- Union of transcripts in affy array and RNAseq
  (select porcine_na30_ensembl.transcript_id AS transcript_id from porcine_na30_ensembl union select tmp_edger_transcript.transcript_id from tmp_edger_transcript) as u_transcripts
  -- Add  to each ensembl transcript_id the affy probe_id
  left join porcine_na30_ensembl on porcine_na30_ensembl.transcript_id = u_transcripts.transcript_id -- probe_set_id = limma_toptable.id
  -- Add RNAseq expression data fromedgeR
  left join tmp_edger_transcript on u_transcripts.transcript_id = tmp_edger_transcript.transcript_id
  -- Add RNAseq FPKM
  left join tmp_cufflinks_transcript_gtf_ct on u_transcripts.transcript_id = tmp_cufflinks_transcript_gtf_ct.transcript_id
  -- toptable() probes
  left join limma_toptable on porcine_na30_ensembl.probe_set_id = limma_toptable.id
  -- Add array intensity data
  left join tmp_affyarray_bmdm_ct on tmp_affyarray_bmdm_ct.probe_set_id = porcine_na30_ensembl.probe_set_id
  -- Add annotation from Affymetrix porcinelayout
  left join porcine_na30_annot on porcine_na30_annot.probe_set_id = porcine_na30_ensembl.probe_set_id
  -- Restrict to top 100 transcripts found diff. expr. in RNAseq ranked by logFC
-- where porcine_na30_ensembl.transcript_id in (select transcript_id from (select * from tmp_edger_transcript where "FDR" < 0.01 order by abs("logFC") desc limit 100) AS t1)
order by u_transcripts.transcript_id
);
select * from affy_rnaseq limit 10;

-- Array layout
select count(distinct probe_set_id) from porcine_na30_annot;
select count(distinct probe_set_id) from porcine_na30_annot;
select count(distinct transcript_id) from porcine_na30_ensembl  -- ensembl transcripts on array
select count(distinct probe_set_id) from porcine_na30_ensembl  --  num probes matching transcripts
select count(*), count(distinct gene_symbol) from porcine_na30_annot where gene_symbol not like '---' and gene_symbol is not null;
select distinct gene_symbol from porcine_na30_annot where gene_symbol not like '---' and gene_symbol is not null;

copy affy_rnaseq to E'F:/data/20100413/affy_rnaseq_db.txt' with delimiter E'\t' null '' csv header;


-- Differentially expressed transcripts with an Affy probe associated
select distinct transcript_id, count(probe_set_id) AS no_probes 
from porcine_na30_ensembl 
where transcript_id in 
  (select transcript_id from tmp_edger_transcript where "FDR" < 0.01)
group by transcript_id order by transcript_id;

---------------------------------[ Triutme ]-----------------------------------
select distinct transcript_id from affy_rnaseq union select transcript_id from cufflinks_transcript_gtf where source in ('20100317_RNAseq_CTRL', '20100317_RNAseq_CTRL') and feature like 'transcript';

-- drop table affy_rnaseq ;
create table affy_rnaseq AS(
  select 
    porcine_na30_ensembl.transcript_id, affyarray_bmdm.probe_set_id, avg(intensity) AS avg_intensity, 
    stddev(intensity) AS sd_intensity, 
    100*(stddev(intensity)/avg(intensity)) AS cv_intensity, 
    "logConc", "logFC", "PValue", "FDR", pig, time_point
  from affyarray_bmdm 
    -- Add to the arrays the correpsonding ensembl transcript_id
    inner join porcine_na30_ensembl on affyarray_bmdm.probe_set_id = porcine_na30_ensembl.probe_set_id
    -- Add the expression value from RNAseq for transcripts in common between array and RNAseq
    inner join tmp_edger_transcript on tmp_edger_transcript.transcript_id = porcine_na30_ensembl.transcript_id
    -- Add array experimental layout
    inner join affyarray_design on affyarray_bmdm.slide_id = affyarray_design.slide_id
  where
    -- Restrict to differentially expressed RNAseq transcripts
    porcine_na30_ensembl.transcript_id in (select transcript_id from tmp_edger_transcript where "FDR" < 0.01)
    -- Restrict to time-points 0h and 7h
    and time_point in (0, 7)
  group by 
    -- Take the average for each tissue culture replicate
    porcine_na30_ensembl.transcript_id, affyarray_bmdm.probe_set_id, "logConc", "logFC", "PValue", "FDR", pig, time_point
  order by transcript_id, time_point, pig
  );

select * from affy_rnaseq limit 100;
-- Write this table out and go to R
copy affy_rnaseq to 'C:/Tritume/affy_rnaseq.txt' csv header delimiter E'\t';

select * from affy_rnaseq where transcript_id like 'ENSSSCT00000001533'; 

affy_rnaseq

select names($$t:'affyarray_bmdm', quote:''$$) -- "probe_set_id, slide_id, intensity"
select names($$t:'affyarray_design', quote:''$$) -- "slide_id, pig, time_point, cell_culture_replicate"

drop table tmp_porcine_na30_annot;
select names($$ t:'tmp_edger_transcript', quote:'"' $$) -- ""transcript_id", "logConc", "logFC", "PValue", "FDR""
select names($$ t:'affyarray_design', quote:'' $$) -- "slide_id, pig, time_point, cell_culture_replicate"
select names($$ t:'affyarray_bmdm', quote:'' $$) -- "probe_set_id, slide_id, intensity"

select * from porcine_na30_annot where probe_set_id like 'AFFX-r2-Ec-bioC-3_at';


select "Chromosomal Location", count(*)
from tmp_porcine_na30_annot
group by "Chromosomal Location";

select count(distinct "RefSeq Transcript ID") from tmp_porcine_na30_annot where "RefSeq Transcript ID" not like '---';
select distinct ensembl from tmp_porcine_na30_annot where ensembl like 'ENSSSCT00000001533';
select count(distinct "Representative Public ID") from tmp_porcine_na30_annot where "Representative Public ID" not like '---';
select count(distinct "Ensembl Transcript ID") from tmp_porcine_na30_ensembl;

select * from tmp_porcine_na30_annot where gene_symbol like 'IDO1' limit 100;
select * from tmp_porcine_na30_annot limit 100;
select count(*) from tmp_porcine_na30_annot where entrez_gene not like '---';

select names($$ t:'tmp_test_rt', apply:("re.sub('[dt]', 'x', x)", '.capitalize()'), which:[4], echo:True, rename:True $$)

select * from tmp_test_rt;

select * from tmp_edger_transcript where transcript_id like 'ENSSSCT00000011156';
select * from cufflinks_transcript_gtf where transcript_id like 'ENSSSCT00000011156';