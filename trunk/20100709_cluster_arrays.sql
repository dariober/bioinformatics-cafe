select * from porcine_na30_annot where gene_symbol like '%TNF%';

drop table go_terms;
create temp table go_terms AS(
  select distinct ensembl_transcript_id AS transcript_id, array_agg(go_term_name_bp) AS go_term_bp
  from (SELECT DISTINCT ensembl_transcript_id, go_term_name_bp FROM ensembl_annotation_sscrofa9_57 ORDER BY ensembl_transcript_id, go_term_name_bp) AS t1
  group by ensembl_transcript_id
  );
select * from go_terms limit 100;


create temp table gene_cluster AS(
-- Link clustered probes to the affymetrix annotation and to the differential expression analysis from limma
select 
       array_agg(porcine_na30_ensembl.transcript_id) AS transcript_ids,
       array_agg(tmp_hopach.probe_set_id) AS probe_set_ids,
       -- Cluster info
       tmp_hopach.label AS cluster_id,  -- Coordinates of the larger clusters in the tree
       tmp_hopach.label_final AS cluster_final_id,  -- All the coordinates of each cluster
       -- Annotation from Affymetrix/ensembl
       porcine_na30_annot.gene_symbol,
       -- From limma toptable
       avg(logfc) AS avg_logfc, 
       avg(aveexpr) AS avg_expr,
       avg(adj_p_val) AS avg_fdr,
       CASE WHEN avg(adj_p_val) > 0.1 THEN '' 
            WHEN avg(adj_p_val) < 0.1 and avg(adj_p_val) >= 0.05 THEN '+' 
            WHEN avg(adj_p_val) < 0.05 and avg(adj_p_val) >= 0.01 THEN '*' 
            WHEN avg(adj_p_val) < 0.01 and avg(adj_p_val) >= 0.001 THEN '**' 
            WHEN avg(adj_p_val) < 0.001 and avg(adj_p_val) >= 0.0001 THEN '***'
            WHEN avg(adj_p_val) < 0.0001 THEN '****' END AS diff_expr_signf_code,
       -- GO terms 
       replace(
           replace(porcine_na30_annot.gene_ontology_biological_process,  ' // inferred from electronic annotation', ''),
           ' // ', ': '
           ) AS go_term_bp
  from tmp_hopach 
      inner join porcine_na30_annot on 
          tmp_hopach.probe_set_id = porcine_na30_annot.probe_set_id
      inner join limma_toptable on 
          id = porcine_na30_annot.probe_set_id
      left join porcine_na30_ensembl on 
          porcine_na30_ensembl.probe_set_id = porcine_na30_annot.probe_set_id
  -- Restrict to probes with gene_symbol or transcript_id
  where porcine_na30_annot.gene_symbol not like '---' or porcine_na30_ensembl.transcript_id is not null  /* or gene_ontology_biological_process not like '---' */
  group by tmp_hopach.label, 
           porcine_na30_annot.gene_symbol,
           porcine_na30_annot.gene_ontology_biological_process,
           tmp_hopach.label_final
  order by cluster_final_id, avg_logfc
);
COPY gene_cluster TO 'D:/Tritume/gene_cluster.txt' with csv header delimiter E'\t';

select substring(cluster_final_id::text from 1 for 3), 
       count(*) AS no_features,
       sum(CASE WHEN go_term_bp like '%0006955: %' THEN 1 ELSE 0 END) AS "no_immuneResponse",
       round((sum(CASE WHEN go_term_bp like '%0006955: %' THEN 1 ELSE 0 END) / count(*)::numeric)*100, 2) AS "frac_immuneResponse"
from gene_cluster 
group by substring(cluster_final_id::text from 1 for 3)
having count(*) > 10
order by "frac_immuneResponse" desc;



select names('tmp_hopach');
"feature_order", "label", "label_final", "feature", "probe_set_id"

select names('limma_toptable');
"dataset_id", "id", "logfc", "aveexpr", "t", "p_value", "adj_p_val", "b"


select * from tmp_hopach where feature like '%PGHS%';
select * from porcine_na30_annot where probe_set_id like 'Ssc.7314.1.A1_at';
select * from porcine_na30_ensembl where probe_set_id like 'Ssc.7314.1.A1_at';