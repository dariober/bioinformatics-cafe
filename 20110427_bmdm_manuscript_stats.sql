/*
   Various stats for BMDM manuscript RK
*/
-- Summary of gene expression: Probes & genes up and down regulated at each time point.
select 'LPS induced' as change, coef, count(id) as no_probes, count(distinct gene_name) as no_genes
from limma_toptable left join (select * from annotation_ctuggle where gene_name is not null) AS annotation_ctuggle on 
    limma_toptable.id = annotation_ctuggle.probe_set_id
where adjpval < 0.01 and logfc > 0
group by coef 
union 
select 'LPS repressed' as change, coef, count(id), count(distinct gene_name) as no_genes
from limma_toptable left join (select * from annotation_ctuggle where gene_name is not null) AS annotation_ctuggle on 
    limma_toptable.id = annotation_ctuggle.probe_set_id
where adjpval < 0.01 and logfc < 0 group by coef
order by coef, change;

-- Fold change for some genes of interest
-- 2 hours
select gene_name, coef, avg(logfc) as fold_change, count(id)
from limma_toptable inner join annotation_ctuggle on id = probe_set_id
where gene_name in ('IL1B', 'TNF', 'CXCL9', 'IL1A', 'CXCL2') and coef = '2vs0'
group by gene_name, coef
order by gene_name, coef;

-- 7 hours
select gene_name, coef, avg(logfc) as fold_change, count(id)
from limma_toptable inner join annotation_ctuggle on id = probe_set_id
where gene_name in ('CXCL11', 'CXCL9', 'TNFSF10', 'IL1B', 'CXCL13') and coef = '7vs0'
group by gene_name, coef
order by gene_name, coef;

-- 24 hours
select gene_name, coef, avg(logfc) as fold_change, count(id)
from limma_toptable inner join annotation_ctuggle on id = probe_set_id
where gene_name in ('CXCL9', 'CXCL13') and coef = '24vs0'
group by gene_name, coef
order by gene_name, coef;

select gene_name, coef, avg(logfc) as fold_change, count(id)
from limma_toptable inner join annotation_ctuggle on id = probe_set_id
where gene_name in ('GADD45A') and coef = '2vs0'
group by gene_name, coef
order by gene_name, coef;

-- Top down regulated:
select gene_name, coef, avg(logfc) as fold_change, avg(adjpval) as adjpval, count(id)
from limma_toptable inner join annotation_ctuggle on id = probe_set_id
where coef = '7vs0' and adjpval < 0.01
group by gene_name, coef
order by fold_change
limit 100;

select gene_name, coef, avg(logfc) as fold_change, avg(adjpval) as adjpval, count(id)
from limma_toptable inner join annotation_ctuggle on id = probe_set_id
where coef = '24vs0' and adjpval < 0.01
group by gene_name, coef
order by fold_change
limit 100;

-- Supplemenatry Table: All probes and genes
select gene_name, limma_toptable.* 
from limma_toptable left join annotation_ctuggle on id = probe_set_id 
order by coef, adjpval
limit 100;

------------------------[ KEGG pathway table ]---------------------------------

-- For convenience, make a table with the pathways of interest
drop table if exists selected_pathways;
create temp table selected_pathways (
    pathway_name text
);
insert into selected_pathways values ('Antigen processing and presentation');
insert into selected_pathways values ('Cell adhesion molecules (CAMs)');
insert into selected_pathways values ('Chemokine signaling pathway');
insert into selected_pathways values ('Cytokine-cytokine receptor interaction');
insert into selected_pathways values ('Endocytosis');
insert into selected_pathways values ('Lysosome');
insert into selected_pathways values ('MAPK signaling pathway');
insert into selected_pathways values ('NOD-like receptor signaling pathway');
insert into selected_pathways values ('Phagosome');
insert into selected_pathways values ('RIG-I-like receptor signaling pathway');
insert into selected_pathways values ('TGF-beta signaling pathway');
insert into selected_pathways values ('Toll-like receptor signaling pathway');
select * from selected_pathways;

-- Table of KEGG pathway names and gene names:
create or replace view affymetrix.attract_pathways AS(
select distinct attract_synexpr_groups_kegg.keggid, pathway, cluster_group, annotation_ctuggle.gene_name, ensembl_gene_id
from
    -- Add kegg pathway name 
    attract_synexpr_groups_kegg inner join attract_ranked_pathways on
        attract_synexpr_groups_kegg.keggid = attract_ranked_pathways.keggid
    -- Add gene_name from C.Tuggle annotation
    inner join annotation_ctuggle on
        attract_synexpr_groups_kegg.probe_set_id = annotation_ctuggle.probe_set_id
    -- Add ensembl_gene_id if known:
    left join annotation_ctuggle_ensembl on 
        annotation_ctuggle.gene_name = annotation_ctuggle_ensembl.gene_name
where
-- Restirct to probes annotated to KEGG (1:yes, 0:No)?
kegg_probe = 1 and
-- Restrict to probes annotated with gene_name and/or ensembl_gene_id?
annotation_ctuggle.gene_name is not null
-- Restrict to pathways of interest?
-- pathway in (select * from selected_pathways)
order by pathway, cluster_group
);
comment on view attract_pathways is 'Gene clusters produced by attract within KEGG pathways and annotated with gene names. See 20110427_bmdm_manuscript_stats.sql';

select pathway, cluster_group || ' (' || count(gene_name) || ')' as cluster, array_to_string(sort(array_agg(gene_name || CASE WHEN ensembl_gene_id is not null THEN '*' ELSE '' END)), ', ') as gene_name
from attract_pathways
where pathway in (select * from selected_pathways)
group by pathway, cluster_group
order by pathway, cluster_group;

-- KEGG pathways significantly enriched after LPS as (approx.) lm(gene_expression ~ time_point):
select * from attract_ranked_pathways;
select * from attract_ranked_pathways where pvalue < 0.01;
---------------------------------[ TRITUME ]-----------------------------------

select count(*) from porcine_na30_annot where probe_set_id like 'Ssc.24282.1.S1_at';
select count()



-----------------------------[ TRITUME ]---------------------------------------

select * from porcine_na30_annot where gene_symbol like '%CXC%'

select * from probe2transcript where ensembl_transcript_id like 'ENSSSCT00000009833'

select * from annotation_ctuggle where gene_name like '%CXC%'