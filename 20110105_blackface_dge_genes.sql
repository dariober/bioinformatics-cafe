SET search_path='blackface_dge', 'public', 'Pigs';

CREATE TABLE gtf_coverage
(
  library_id text,
  rname text,
  source text,
  feature text,
  f_start integer,
  f_end integer,
  score text,
  strand text,
  frame text,
  attributes text,
  no_tags int,
  no_nonzero_cov int,
  feature_length int,
  cov_ratio double precision
);
comment on table gtf_coverage is 'Coverage file produced by coverageBed see 09/01/2011.'
truncate table gtf_coverage;
select read_table($$ file:'F:/data/20101223_blackface/LN_all.nonzero.cov.gz', table:'gtf_coverage', append:True $$);
alter table gtf_coverage set tablespace hdd_free_agent;
select * from gtf_coverage limit 10;

drop table transcript_coverage;
create table gene_coverage AS (
  select library_id, min(rname) AS rname, min(f_start) AS f_start, max(f_end) AS f_end, min(strand) AS strand, 
         gtf_attribute(attributes, 'transcript_id') AS gene_id, sum(no_tags) AS coverage
  from gtf_coverage
  where feature like 'exon'  group by library_id, gtf_attribute(attributes, 'transcript_id')
  );
select count(distinct gene_id) from gene_coverage;
select library_id, count(distinct gene_id) AS no_genes from gene_coverage group by library_id;
select * from gene_coverage order by gene_id, library_id limit 200;

drop table gene_coverage_ct;
select cross_tab('select gene_id, library_id, coverage from gene_coverage', 'gene_coverage_ct');
select * from gene_coverage_ct limit 100;

-- Gene expression analysis done in R/edgeR and output of toptags()to edger_toptags 20110112_blackface_dge_gene_de.R


