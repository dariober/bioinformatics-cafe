/*
  Import and analyze dge tags from Blackface project.
*/

-- Import tags mapped and grouped by position within libraries.
drop table tag_expression;
select read_table($$ file:'F:/data/20101223_blackface/blackface_tag_count.csv.gz', dest_table: 'tag_expression', header:True, data_type:'text, text, int, int, int' $$)

-- Rename these libraries to have the format of the IDs all the same
update tag_expression set library_id = 'LN_114' where library_id = 'LN114';
update tag_expression set library_id = 'LN_124' where library_id = 'LN124';

-- Generate tag IDs
alter table tag_expression add column tag_id text;
update tag_expression set tag_id = rname || '_' || pos || '_' || strand;

-- Convert raw expression to Tags per Million:
alter table tag_expression add column tag_tpm double precision;
update tag_expression set tag_tpm = (tag_count::numeric/tot_library_count)*1000000
from 
    (select library_id, sum(tag_count) AS tot_library_count from tag_expression group by library_id) AS tot_expr
where tag_expression.library_id = tot_expr.library_id;

alter table tag_expression add constraint pk_tag_lib primary key (tag_id, library_id);
-- vacuum analyze tag_expression;
select * from tag_expression limit 1000;
   
-- Collect tags expressed in at least one group with more than a given tpm threshold:
drop table tmp_sig_tags;
create table tmp_sig_tags AS(
  select distinct tag_expression.tag_id
--       treatment_groups.treatment_group, 
--       sum(tag_tpm) AS group_tpm
  from treatment_groups inner join tag_expression on tag_expression.library_id = treatment_groups.library_id
  group by treatment_groups.treatment_group, tag_expression.tag_id
  having sum(tag_tpm) > 2
  );
select * from tmp_sig_tags limit 100;
select count(*) from tmp_sig_tags;
alter table tmp_sig_tags add constraint pk_sig_tags primary key (tag_id);
-- vacuum analyze sig_tags;

-- Generate cross-tab (good for edgeR,  DEseq)
drop table sig_tags_ct;
select cross_tab('SELECT tag_expression.tag_id, library_id, tag_count  
                  FROM tag_expression INNER JOIN tmp_sig_tags ON tag_expression.tag_id = tmp_sig_tags.tag_id', 'sig_tags_ct');
comment on table sig_tags_ct is 'Tag count in each of the DGE libraries in cross-tab format. Only tags having expression greater than a 2 tpm (cumulated by treatment group) included. See script 20101224_blackface_dge.sql'
select * from sig_tags_ct limit 10;

------------------------[ Tritume ]--------------------------
select tag_expression.tag_id,
     treatment_groups.treatment_group, 
     sum(tag_tpm) AS group_tpm
from treatment_groups inner join tag_expression on tag_expression.library_id = treatment_groups.library_id
group by treatment_groups.treatment_group, tag_expression.tag_id
limit 100;


select min(treatment_group) AS treatment_group, tag_expression.library_id, sum(tag_count) AS total_expression, count(tag_id) AS no_tags 
from tag_expression inner join treatment_groups on tag_expression.library_id = treatment_groups.library_id
group by tag_expression.library_id
order by treatment_group, no_tags;



copy (select * from "qrySolexaLambs") to 'D:/Tritume/lambs.txt' with csv header delimiter E'\t';

drop table treatment_groups;
select read_table($$ file:'D:/Tritume/lambs.txt', dest_table:'treatment_groups', data_type:'text, text, int, int, double precision, int, text', header:True $$)
alter table treatment_groups add constraint pk_treatment_group primary key (library_id);
select * from treatment_groups;
-- vacuum analyze treatment_groups;
