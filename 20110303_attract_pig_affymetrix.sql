/*
  Parse and analyze the output of R/attract.
  See 20110303_attract_pig_affymetrix.R for how the input table have been produced.
*/

select * from affymetrix.attract_syncorr limit 10;
select * from affymetrix.attract_synexpression;

-- Reformat table synexpression to have one probe per line. Then annotate probes with gene symbols.
drop table if exists synexpression_long;
create temp table synexpression_long AS(
    select keggid, pathway, cluster_group, (string_to_array(syn_groups, ', '))[i] AS probe_set_id
    from affymetrix.attract_synexpression, 
      generate_series(1, ( select max(array_upper(string_to_array(syn_groups, ', '), 1))  from affymetrix.attract_synexpression) ) i
    where (string_to_array(syn_groups, ', '))[i] is not null
    order by pathway, cluster_group, probe_set_id
    );

-- Annotation:
copy (
select pathway, 
    case when cluster_group = 1 then pathway else NULL END AS pathway_print, -- This column is for pretty printing only
    cluster_group, count(*) AS no_genes, array_to_string(array_agg(gene_symbol), ', ' ) as genes from 
   -- This inner query removes duplicate gene names (different probes mapping to the same gene in the same cluster group)
   (select distinct pathway, cluster_group, gene_symbol from synexpression_long inner join affymetrix.affymetrix_consensus_annot
   on synexpression_long.probe_set_id = affymetrix.affymetrix_consensus_annot.probe_set_id
   order by pathway, cluster_group, gene_symbol) as t1
group by pathway, cluster_group
order by pathway, cluster_group
) to 'D:/Tritume/synexpression_clusters.txt' with csv header delimiter E'\t';

-------------------------------------------------------------------------------
/* Gene expression levels
   Assing to each gene in the above clusters their gene expression level */
-------------------------------------------------------------------------------

drop table if exists gene_expression;
create temp table gene_expression AS(
  -- Convert probe expression to gene expression using the annotation file
  select array_id, pig, time_point, avg(log2_intensity) as log2_intensity, count(*) as no_probes, gene_symbol
  from affyarray_bmdm_avg inner join affymetrix_consensus_annot on 
      affyarray_bmdm_avg.probe_set_id = affymetrix_consensus_annot.probe_set_id
  group by array_id, pig, time_point, gene_symbol
);
select * from gene_expression limit 10;
create unique index ind_gene_expression on gene_expression (gene_symbol, time_point, pig );

drop table if exists gene_expression_cluster;
create table public.gene_expression_cluster AS(
select distinct synexpression_long.pathway, cluster_group, gene_expression.*
from synexpression_long 
  -- Link probes in synexpression groups to gene symbols:
    inner join affymetrix_consensus_annot on 
        affymetrix_consensus_annot.probe_set_id = synexpression_long.probe_set_id
  -- Link gene symbols to expression:
    inner join gene_expression on
        affymetrix_consensus_annot.gene_symbol = gene_expression.gene_symbol
order by pathway, cluster_group, gene_symbol, time_point, pig
);
select count(*) from gene_expression_cluster where gene_symbol = 'B2M'; 

/*
  Cross table gene expression suitable for plotting clusters
*/
drop table if exists tmp_cluster_ct;
select cross_tab($$ 
    select distinct
   --     pathway, 
   --     cluster_group, 
        gene_symbol, array_id, log2_intensity from gene_expression_cluster 
        where pathway = 'Apoptosis' and cluster_group >= 1
        $$, 'tmp_cluster_ct');
alter table tmp_cluster_ct set schema public;

copy tmp_cluster_ct to 'D:/Tritume/cluster_attract.expression' with csv header delimiter E'\t';

select distinct pathway from tmp_cluster_ct order by pathway limit 100;