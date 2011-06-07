/*
  Comparison Affymetrix vs RNAseq
*/

-- Annotate Affymetrix probes with ensembl gene_id
-- Select 
drop table if exists affy_annot;
create temp table affy_annot AS(
select time_point, avg(log2_intensity) as log2_intensity, ensembl_gene_id, count(*)
  -- Link probe to gene_name
  from affyarray_bmdm_avg inner join annotation_ctuggle on 
      affyarray_bmdm_avg.probe_set_id = annotation_ctuggle.probe_set_id
  -- Link gene_name to ense,mbl gene_id
  inner join (select distinct gene_name, ensembl_gene_id from annotation_ctuggle_ensembl) AS annotation_ctuggle_ensembl on 
      annotation_ctuggle.gene_name = annotation_ctuggle_ensembl.gene_name
where pig = 'LWL1' and time_point in (0, 7) -- LWL1 is the one used for rnaseq
group by time_point, ensembl_gene_id
order by ensembl_gene_id, time_point
);

select * from affy_annot limit 100;
select * from deseq_nbinom where id = 'ENSSSCG00000017892'

copy (
select 
  affy_0.ensembl_gene_id, 
  affy_0.log2_intensity AS log2intensity_0h, 
  affy_7.log2_intensity AS log2intensity_7h, 
  affy_7.log2_intensity - affy_0.log2_intensity AS log2fc_affy, 
  case when deseq_nbinom.log2foldchange is null then 0 
       when deseq_nbinom.log2foldchange = 'Infinity' then 20
       when deseq_nbinom.log2foldchange = '-Infinity' then -20 
       else deseq_nbinom.log2foldchange END AS log2fc_rnaseq
from       (select * from affy_annot where time_point = 0) as affy_0
inner join (select * from affy_annot where time_point = 7) as affy_7 on
    affy_7.ensembl_gene_id = affy_0.ensembl_gene_id
left join (select * from deseq_nbinom where deseq_nbinom.comparison = 'bmdm') as deseq_nbinom on 
    deseq_nbinom.id = affy_7.ensembl_gene_id
) to 'D:/Tritume/affy_rnaseq.txt' with csv header delimiter E'\t';
