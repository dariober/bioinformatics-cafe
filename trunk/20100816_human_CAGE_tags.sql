/*
  Import and parse human CAGE tags from FANTOM4
*/


-- Import
select read_table($$  file:'D:/Tritume/fantom4_human_cage_tags.txt', sep:'\t', header: True, dest_table: 'cage_tags_human', limit: -1, overwrite: False, data_type: 'varchar, varchar, bigint' $$)
comment on table cage_tags_human is $$ All human CAGE tags downloaded from http://fantom.gsc.riken.jp/4/download/Tables/human/CAGE/tag_counts/ and concatenated with fantom4_cage_tags_concatenate-1.0.py. Column tag_count is the total tag count. See LabBook 16/08/2010 $$;
comment on column cage_tags_human.tag_count IS 'Total tag count (last column from in each concatenated file)';

alter table cage_tags_human set schema "FANTOM4";

select * from cage_tags_human limit 120;

-- All unique tags (this file sent to Bowtie after having converted it to FASTA)

copy (select distinct tag from cage_tags_human) to 'D:/Tritume/human_cage_tags_unique.txt';

/*
----------------[  Convert D:/Tritume/human_cage_tags_unique.txt to FASTA format with the following python ] --------------------

f= open('D:/Tritume/human_cage_tags_unique.txt', 'r')
fout= open('D:/Tritume/human_cage_tags_unique.fa', 'w')

for line in f:
    fout.write('>' + line)
    fout.write(line)
    
f.close()
fout.close()

*/


-- Library sizes and standardization (tags per million reads)

drop table library_size;
create temp table library_size AS (
  select dataset_id, sum(tag_count) AS library_size from cage_tags_human group by dataset_id
  );
select * from library_size;

alter table cage_tags_human add column tpm double precision;

update cage_tags_human set tpm = t1.tpm 
-- This takes ages
from ( 
    select cage_tags_human.dataset_id, cage_tags_human.tag, (tag_count / library_size) * 10^6 AS tpm   
    from cage_tags_human inner join library_size on 
    cage_tags_human.dataset_id = library_size.dataset_id
    ) as t1
where t1.dataset_id = cage_tags_human.dataset_id and 
      t1.tag = cage_tags_human.tag;

-- Upload SAM file from bowtie (human CAGE tags against pig). Mapped reads only. Single mappers
-- Put this in Pig schema
select read_sam($$ file:'D:/Tritume/20100816_human_cage_vs_sscrofa9.sam.aln', dest_table: 'sam_bowtie_human_vs_pig', dataset_id: '20100816_human_cage_vs_sscrofa9' $$)

-- Upload SAM file from bowtie (human CAGE tags against pig). All reads (also not mapped). Multi-mappers (-m 10)
-- Put this in Pig schema

-- drop table sam_bowtie_human_vs_pig;
select read_sam($$ file:'D:/Tritume/20100817_human_cage_vs_sscrofa9.sam', dest_table: 'sam_bowtie_human_vs_pig', dataset_id: '20100817_human_cage_vs_sscrofa9' $$)
select names('sam_bowtie_human_vs_pig');


drop table cage_tags;
create temp table cage_tags AS(
  select "qname", 
          rname, 
          CASE WHEN flag = 0 THEN '+' WHEN flag = 16 THEN '-' END AS strand,
          pos AS "start", 
          pos + length(seq) AS "end", 
          (abs(pos - (pos + length(seq)))/2) + pos AS mid
  from sam_bowtie_human_vs_pig where flag = 0  AND or flag = 16
  limit 100000
  );

copy cage_tags to 'D:/Tritume/cage_tags.txt' with csv delimiter E'\t';

drop index cage_fstrand;
create index cage_midrname on cage_fstrand (mid, rname);

drop table gtf;
create table gtf AS (
    select rname, source, f_start, f_end, score, strand, frame, attributes, 
       gtf_attribute(attributes, 'transcript_id') AS transcript_id,
       gtf_attribute(attributes, 'exon_number') AS exon_number,
       gtf_attribute(attributes, 'transcript_id') || '_ex_' || gtf_attribute(attributes, 'exon_number') AS exon_id
    from sus_scrofa_sscrofa9_56_gtf where feature like 'exon'
    );
copy gtf to 'D:/Tritume/sscrofa9_56gtf.txt' with csv delimiter E'\t';
drop table gtf;

create index gtf_f_startrname on gtf_f (f_start, rname);

select cage_fstrand.*, gtf_f.* 
from cage_fstrand, gtf_f
where 
     gtf_f.f_start = cage_fstrand.mid and 
     cage_fstrand.rname like '18' and
     gtf_f.rname like '18';

-- Number of mismatches
select sam_get_tagvalue(tags, 'NM') AS n_mismatches, count(*) from sam_bowtie_human_vs_pig group by sam_get_tagvalue(tags, 'NM') limit 10 ;

-- Tag position relative to ENSEMBL transcripts


create temp table ensembl
select distinct sam_bowtie_human_vs_pig.rname,
       qname,
       pos AS tag_start,
       pos + length(seq) AS tag_end,
       gtf.strand,
       gtf_attribute(attributes, 'transcript_id') AS transcript_id,
       gtf_attribute(attributes, 'exon_number') AS exon_number,
       gtf.f_start,
       gtf.f_end
  from (select * from sam_bowtie_human_vs_pig limit 1000) AS sam_bowtie_human_vs_pig inner join gtf ON
    sam_bowtie_human_vs_pig.rname = gtf.rname and
    CASE WHEN flag = 0 THEN '+' WHEN flag = 16 THEN '-' END = gtf.strand
  where
    (gtf.strand like '+' AND (gtf.f_start - (pos + length(seq))   between -1000 and +1000 OR 
                              gtf.f_start - pos between -1000 and +1000) 
                              ) OR
    (gtf.strand like '-' AND (gtf.f_end - (pos + length(seq))    between -1000 and +1000 OR
                              gtf.f_end - pos  between -1000 and +1000) )


----------------------------[ Cluster CAGE tags ]------------------------------

select cluster_reads($$ select rname, 
                               CASE WHEN flag = 0 THEN '+' WHEN flag = 16 THEN '-' ELSE '.' END AS strand, 
                               pos AS read_start,
                               pos + length(seq) AS read_end 
                        from sam_bowtie_human_vs_pig $$,
                        'tmp_cage_hvsp_clustered'
                    );

select count(distinct tag_cluster) from tmp_cage_hvsp_clustered;

-- drop table cluster_cage_hvsp;
create temp table cluster_cage_hvsp AS(
  select tag_cluster, 
       rname, 
       strand, 
       min(read_start) AS cluster_start, 
       max(read_end) AS cluster_end, 
       max(read_end) - min(read_start) AS cluster_length, 
       count(*) AS no_clustered_tags 
  from tmp_cage_hvsp_clustered
  where rname not like 'gi|555853|gb|U13369.1|HSU13369' 
  group by tag_cluster, rname, strand
  );

select * from cluster_cage_hvsp order by cluster_length desc limit 200;

select count(*) from cluster_cage_hvsp;


------------------------------[ Clusters vs RNAseq: Position ]-----------------

create table tmp_tss_rnaseq AS(
  select distinct cluster_cage_hvsp.rname,
       tag_cluster,
       cluster_start,
       cluster_end,
       gtf.strand,
       transcript_id,
       gtf.start,
       gtf.end
  from cluster_cage_hvsp inner join (select * from cufflinks_transcript_gtf where feature like 'transcript' and source in ('20100602_LPS2_gtf', '20100602_CTRL2_gtf') ) AS gtf ON
    cluster_cage_hvsp.rname = gtf.rname and
    cluster_cage_hvsp.strand = gtf.strand
  where
    (gtf.strand like '+' AND (gtf.start - cluster_end   between -1000 and +1000 OR 
                              gtf.start - cluster_start between -1000 and +1000) 
                              ) OR
    (gtf.strand like '-' AND (gtf.end - cluster_end    between -1000 and +1000 OR
                              gtf.end - cluster_start  between -1000 and +1000) )
  );
comment on table tmp_tss_rnaseq is 'Human CAGE tags (single-mappers) mapped to Sscrofa9, clustered, and mapping near ensembl/RNAseq transcripts. See 20100816_human_CAGE_tags.sql' 

select * from tmp_tss_rnaseq where transcript_id like 'ENSSSCT00000007031';
select *, "start" - cluster_start from tmp_tss_rnaseq where transcript_id like 'ENSSSCT00000001533';



select * from sam_bowtie_human_vs_pig where rname like '4' and flag = 16 and pos between 94928458-1000 and 94956044+1000 order by pos;

select * from cage_tags_human where tag in (select qname from sam_bowtie_human_vs_pig where rname like '4' and flag = 0 and pos between 94928458-1000 and 94956044+1000 order by pos);

select names('sam_bowtie_human_vs_pig');
select distinct source from cufflinks_transcript_gtf order by source;

------------------------------------[ Tritume ]--------------------------------

select * from sam_bowtie_human_vs_pig limit 20;
select distinct rname from sam_bowtie_human_vs_pig;
select count(*) from sam_bowtie_human_vs_pig where rname like 'gi|555853|gb|U13369.1|HSU13369';

select 16 & 128;

