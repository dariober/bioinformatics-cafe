select read_sam($$ file:'D:/Tritume/20101026_cage_clones.sam', dest_table: 'tmp_inserts_oct10', overwrite:False $$)
alter table tmp_inserts_oct10 add column tss integer;
update tmp_inserts_oct10 set tss = CASE WHEN flag = 0 THEN pos WHEN flag = 16 THEN pos + length(seq) END;

select * from tmp_inserts_oct10 where flag = 0;

select no_hits, count(no_hits) as multi_mappers, round(count(no_hits)::numeric/ (select count(distinct qname) from tmp_inserts_oct10 where rname != '*'), 3) AS percent from 
( 
  select qname, count(qname) AS no_hits 
  from tmp_inserts_oct10
  where rname not like '*'
  group by qname
  order by no_hits desc
  ) AS t
group by no_hits
order by no_hits;


select *, count(*) OVER(PARTITION BY no_mismatches) AS no_reads_in_err_stratum from (
    -- Extract number and position of mismatches. Also add number of hits of each tag.
    select *, 
           sam_get_tagvalue(tags, 'NM') AS no_mismatches, 
           sam_get_tagvalue(tags, 'MD') AS mismatch_pos, 
           count(*) OVER(PARTITION BY qname) as no_hits
    from tmp_inserts_oct10 
    where rname not like '*'
    ) AS t
where no_hits =1
    order by no_mismatches;

drop table mism_pos;
create temp table mism_pos AS(
  select *, 
    sam_get_tagvalue(tags, 'MD') AS mismatch_string,
    CASE WHEN sam_get_tagvalue(tags, 'NM')::int = 0 THEN ARRAY[NULL]
         WHEN sam_get_tagvalue(tags, 'NM')::int = 1 THEN 
             ARRAY[(regexp_split_to_array(sam_get_tagvalue(tags, 'MD'), 'A|C|G|T'))[1]]
         WHEN sam_get_tagvalue(tags, 'NM')::int = 2 THEN 
             (regexp_split_to_array(sam_get_tagvalue(tags, 'MD'), 'A|C|G|T'))[1:2] END AS mism_pos,
    regexp_split_to_array(sam_get_tagvalue(tags, 'MD'), 'A|C|G|T')       
  from tmp_inserts_oct10 
  where sam_get_tagvalue(tags, 'NM'):: int > 0
);
select * from mism_pos;

select *, unnest(mism_pos) from mism_pos;

drop table tss_gtf;
create temp table tss_gtf AS(
    -- Table of most 5' exons in each gene.
    -- This SELECT gets the + strand, where the tss is the start of the feature
    select *, gtf_attribute(attributes, 'transcript_id') AS transcript_id, f_start AS tss
    from sus_scrofa_sscrofa9_59_gtf
    where gtf_attribute(attributes, 'exon_number')::integer = 1 and 
          feature like 'exon' and 
          strand like '+'
          
    -- Get the tss of the genes on the - strand (f_end)
    UNION SELECT *, gtf_attribute(attributes, 'transcript_id') AS transcript_id, f_end AS tss
    from sus_scrofa_sscrofa9_59_gtf
    where gtf_attribute(attributes, 'exon_number')::integer = 1 and 
          feature like 'exon' and 
          strand like '-'
    );

select * from tss_gtf order by transcript_id limit 120;

select dense_rank() OVER(PARTITION BY 1 ORDER BY qname) AS read_index_no,
       count(*) OVER(PARTITION BY qname) AS no_transcripts,
       tmp_inserts_oct10.qname, 
       tmp_inserts_oct10.flag,
       tmp_inserts_oct10.tss AS cage_tss,
       gtf_attribute(attributes, 'gene_name') AS gene_name,
       tss_gtf.*,
       CASE WHEN strand like '+' THEN tss_gtf.tss - tmp_inserts_oct10.tss
            WHEN strand like '-' THEN tmp_inserts_oct10.tss - tss_gtf.tss END AS tss_diff
from tmp_inserts_oct10 inner join tss_gtf on tmp_inserts_oct10.rname = tss_gtf.rname AND
                                             (CASE WHEN tmp_inserts_oct10.flag = 0 THEN '+' WHEN flag = 16 THEN '-' END) = tss_gtf.strand
where tmp_inserts_oct10.tss between tss_gtf.tss - 1000 and tss_gtf.tss + 1000
order by qname;

-- Number of features in gtf files:
select 'Pig' AS species, round(count(*)::numeric/1e6, 2) AS no_features from sus_scrofa_sscrofa9_59_gtf AS no_features
union select 'H_sapiens', round(count(*)::numeric/1e6, 2) from homo_sapiens_grch37_59_gtf
union select 'Mouse', round(count(*)::numeric/1e6, 2) from mus_musculus_37_59_gtf;

select * from sus_scrofa_sscrofa9_59_gtf limit 100;

select distinct feature, source from sus_scrofa_sscrofa9_59_gtf;

select * from sus_scrofa_sscrofa9_59_gtf where source like 'rRNA' limit 100;
select * from sus_scrofa_sscrofa9_59_gtf where source like 'Mt_tRNA' limit 100;


select tmp_inserts_oct10.*, sus_scrofa_sscrofa9_59_gtf.*,
       dense_rank() OVER(PARTITION BY 1 ORDER BY qname) AS read_index_no,
       count(*) OVER(PARTITION BY qname) AS no_transcripts
from tmp_inserts_oct10 inner join sus_scrofa_sscrofa9_59_gtf on 
    tmp_inserts_oct10.rname = sus_scrofa_sscrofa9_59_gtf.rname AND
    (CASE WHEN tmp_inserts_oct10.flag = 0 THEN '+' WHEN flag = 16 THEN '-' END) = sus_scrofa_sscrofa9_59_gtf.strand
where (tmp_inserts_oct10.tss between sus_scrofa_sscrofa9_59_gtf.f_start - 1000 and sus_scrofa_sscrofa9_59_gtf.f_start + 1000) OR
      (tmp_inserts_oct10.tss between sus_scrofa_sscrofa9_59_gtf.f_end - 1000 and sus_scrofa_sscrofa9_59_gtf.f_end + 1000)
