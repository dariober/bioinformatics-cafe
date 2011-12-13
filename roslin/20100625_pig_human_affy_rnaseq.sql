/*
   
  Pig microarray probes corresponding to human genes induced by LPS mapped against RNAseq
  pig probes -> pig genes -> ortholog human genes induced by LPS
*/

-- Pig probes mapping corresponding to induced human genes
select read_table($$ file:'C:/Tritume/pig_7hr_human_top.txt', dest_table:'tmp_pig_7hr_human_top', overwrite:False, header:True $$);

select count(distinct probe_set_id) from tmp_pig_7hr_human_top;

select distinct limma_toptable.* 
   from limma_toptable inner join tmp_pig_7hr_human_top on
   tmp_pig_7hr_human_top.probe_set_id = limma_toptable.id
where dataset_id like 'avgrepl_contr_0h_7h' and 
   adj_p_val < 0.01;

select * from limma_toptable where adj_p_val < 0.01;
select count(*) from limma_toptable ;


-- Link to these probes the corresponding expression level (limma), the ensembl transcript_id, and the corresponding RNAseq expression.
create temp table pig_human AS(
SELECT 
  limma_toptable.dataset_id AS limma_toptable_dataset,
  tmp_pig_7hr_human_top.probe_set_id, 
  limma_toptable.logfc AS logfc_array, 
  limma_toptable.adj_p_val AS fdr_array, 
  --  limma_toptable.id, 
  edger_toptags.dataset_id AS edger_toptags_dataset, 
  porcine_na30_ensembl.transcript_id, 
  edger_toptags."logFC" AS logfc_rnaseq, 
  edger_toptags."FDR" AS fdr_rnaseq
FROM 
  "Pigs".tmp_pig_7hr_human_top, 
  "Pigs".porcine_na30_ensembl, 
  "Pigs".limma_toptable, 
  "Pigs".edger_toptags
WHERE 
  tmp_pig_7hr_human_top.probe_set_id = porcine_na30_ensembl.probe_set_id AND
  porcine_na30_ensembl.probe_set_id = limma_toptable.id AND
  porcine_na30_ensembl.transcript_id = edger_toptags.transcript_id AND
  limma_toptable.dataset_id like 'avgrepl_contr_0h_7h' AND
  edger_toptags.dataset_id like '20100408_LPSvsCTRL_GTF'
  );

select * from pig_human;
select count(distinct transcript_id) from pig_human;
select count(distinct transcript_id) from pig_human where fdr_rnaseq < 0.01;

COPY pig_human to 'C:/Tritume/pig_human_array_affy_rnaseq.txt' with csv delimiter E'\t' header;

	