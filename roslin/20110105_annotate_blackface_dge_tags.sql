-- Annotate DGE tags using Btau GTF file

-- Generate a table of all the tags positions:
-- drop table dge_tags;
create table dge_tags AS(
  select distinct tag_id, rname, pos, (CASE WHEN strand like '16' THEN '-' WHEN strand like '0' THEN '+' END)::text AS strand from tag_expression
  );
select * from dge_tags limit 10;
create unique index ind_dge_tags on dge_tags(pos, rname, strand);
alter table dge_tags add constraint pk_dge_tags primary key (tag_id);
vacuum analyze dge_tags;

create table dge_tags_gtf AS(
  select tag_id, pos, bos_taurus_btau_4_0_60_gtf.* from 
     bos_taurus_btau_4_0_60_gtf inner join dge_tags on 
         bos_taurus_btau_4_0_60_gtf.rname = dge_tags.rname and
         bos_taurus_btau_4_0_60_gtf.strand = dge_tags.strand
  where (dge_tags.strand = '+' and dge_tags.pos    between bos_taurus_btau_4_0_60_gtf.f_start and bos_taurus_btau_4_0_60_gtf.f_end) OR
        (dge_tags.strand = '-' and dge_tags.pos+21 between bos_taurus_btau_4_0_60_gtf.f_start and bos_taurus_btau_4_0_60_gtf.f_end) 
  ); -- 6914297 ms. (115.24 min, rows: 677898 dge_tags x 536297 bos_taurus..gtf)
comment on table dge_tags_gtf is 'DGE tags annotated with GTF file bos_taurus_btau_4_0_60_gtf. A tag is annotated with a GTF feature when it maps on the same chr and strand and its 5-p end is contained within the feature  mapped on the same chromosome and strand of a GTF feature'

/* -- EXPLAIN the create ... select:
Merge Join  (cost=116666.55..170378302.45 rows=653910928 width=196)
  Merge Cond: (((bos_taurus_btau_4_0_60_gtf.strand)::text = dge_tags.strand) AND ((bos_taurus_btau_4_0_60_gtf.rname)::text = dge_tags.rname))
  Join Filter: (((dge_tags.strand = '+'::text) AND (dge_tags.pos >= bos_taurus_btau_4_0_60_gtf.f_start) AND (dge_tags.pos <= bos_taurus_btau_4_0_60_gtf.f_end)) OR ((dge_tags.strand = '-'::text) AND ((dge_tags.pos + 21) >= bos_taurus_btau_4_0_60_gtf.f_start) AND ((dge_tags.pos + 21) <= bos_taurus_btau_4_0_60_gtf.f_end)))
  ->  Sort  (cost=57833.70..59174.44 rows=536297 width=178)
        Sort Key: bos_taurus_btau_4_0_60_gtf.strand, bos_taurus_btau_4_0_60_gtf.rname
        ->  Seq Scan on bos_taurus_btau_4_0_60_gtf  (cost=0.00..6797.87 rows=536297 width=178)
  ->  Sort  (cost=58832.86..60103.91 rows=508421 width=22)
        Sort Key: dge_tags.strand, dge_tags.rname
        ->  Seq Scan on dge_tags  (cost=0.00..10645.57 rows=508421 width=22)
              Filter: ((strand = '+'::text) OR (strand = '-'::text))
*/

limit 10;
select * from dge_tags_gtf limit 1000;
show work_mem;