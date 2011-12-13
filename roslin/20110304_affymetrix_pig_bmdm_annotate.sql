/*
  Affymetirx porcine chips.
  Rank genes from most D.E. and intersect those that are shared by more
  than one comparison (0 vs 2 h, 0 vs 7 h, 0 vs 24 h)
*/

-- Produce a consunsus annotation by merging the gene_symbols from affymetrix with
-- those extracted from humans:
create table affymetrix_consensus_annot AS(
select probe_set_id, gene_symbol, gene_title, 'affymetrix_porcine_na30_annot'::text AS annotation_source from porcine_na30_annot where gene_symbol not like '---'
union 
select probe_set_id, gene_symbol, description, 'human_affyarray_human_annotation'
  from affyarray_human_annotation 
  where probe_set_id not in (select probe_set_id from porcine_na30_annot where gene_symbol not like '---') and
  gene_symbol is not null
  );
comment on table affymetrix_consensus_annot is 'Annotation for affymetrix porcine chip 24K. See 20110304_affymetrix_pig_bmdm_annotate.sql';
alter table affymetrix_consensus_annot set schema affymetrix;
select count(*), annotation_source from affymetrix_consensus_annot group by annotation_source;
/*
 Annotated the output of LIMMA toptable. With the annotation above.       
*/ 

drop table if exists limma_toptable_annotated;
create table public.limma_toptable_annotated AS(
    select gene_symbol,
           comparison,
           avg(logfc) AS avg_logfc, 
           avg(aveexpr) AS avg_expr, 
           avg(adj_p_val) AS avg_adj_p_val, -- Is it correct to take the avg of pvalue? Probably this could be improved... 
           round(avg(logfc)::decimal, 1) AS logfc_round,
              -- This rounding is to give the same diff. expr. rank to probes with very similar fold change.
              -- E.g. by rounding we say that a probe with logfc=2.00001 is equal to one with logfc= 2.00002).
           array_agg(probe_set_id) AS probes
    from limma_toptable inner join affymetrix_consensus_annot on id = probe_set_id
    group by gene_symbol, comparison
    order by comparison
    );

select * from limma_toptable_annotated limit 1000;

-- Up/down-regulated in treated compared to time point 0:
drop table if exists affymetrix_bmdm_genes_de_ct;
select cross_tab($$
  select 
      CASE WHEN avg_logfc > 0 THEN 'upregulated' WHEN avg_logfc < 0 THEN 'downregulated' ELSE 'flat' END::text as regulation, 
      gene_symbol, comparison, avg_logfc
--     dense_rank() over (partition by comparison order by logfc_round desc) AS rank
     from limma_toptable_annotated 
     where avg_adj_p_val < 0.01
  $$, 'affymetrix_bmdm_genes_de_ct');
alter table affymetrix_bmdm_genes_de_ct set schema public;
comment on table affymetrix_bmdm_genes_de_ct is 'See 20110304_affymetrix_pig_bmdm_annotate.sql';

select * from affymetrix_bmdm_genes_de_ct limit 100;

-- Rank genes within comparison groups. Rank genes D.E. (but note that some probes map to multiple genes)
create table public.gene_rank AS(
  select 'upregulated'::text as regulation, gene_symbol, comparison, avg_logfc, avg_adj_p_val, logfc_round,
    dense_rank() over (partition by comparison order by logfc_round desc) AS rank_bin,
    dense_rank() over (partition by comparison order by avg_logfc desc) AS rank,
    probes
    from limma_toptable_annotated 
    where avg_adj_p_val < 0.01 and avg_logfc > 0
  UNION
  select 'downregulated'::text as regulation, gene_symbol, comparison, avg_logfc, avg_adj_p_val, logfc_round,
    dense_rank() over (partition by comparison order by abs(logfc_round) desc) AS rank,
    dense_rank() over (partition by comparison order by abs(avg_logfc) desc) AS rank,
    probes
    from limma_toptable_annotated 
    where avg_adj_p_val < 0.01 and avg_logfc < 0
  order by regulation , comparison, rank
);





---------------------------------[ TRITUME ]-----------------------------------

select * from porcine_na30_annot where probe_set_id like 'Ssc.24282.1.S1_at';

