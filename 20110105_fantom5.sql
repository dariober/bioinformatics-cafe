SET search_path= fantom5, cage, "Pigs", blackface_dge, public;
show search_path;

select read_table($$ 
  header:['rname', 'f_start', 'f_end', 'ctss_tag', 'tag_count', 'strand'],
  data_type:['text', 'int', 'int', 'text', 'int', 'text' ],
  wdir: 'D:/Tritume',
  file:'F:/data/20110105_FANTOM5/Macrophage%20-%20monocyte%20derived%20donor1.CNhs10861.11232-116C8.hg19.ctss.bed.gz', table:'fantom5.macrophage_monocyte_derived_hg19_ctss_bed'
  $$); -- 1,816,577 rows
comment on table macrophage_monocyte_derived_hg19_ctss_bed is 
  $$Position and expression of CAGE trascription start sites from human macrophages monocyte derived mapped against hg19. 
Data downloaded from FANTOM5 https://fantom5-collaboration.gsc.riken.jp/files/data/shared/UPDATE_007/f5pipeline/human.primary_cell.hCAGE/Macrophage%2520-%2520monocyte%2520derived%2520donor1.CNhs10861.11232-116C8.hg19.ctss.bed.gz
Coordinates are 0-based.
$$
select * from macrophage_monocyte_derived_hg19_ctss_bed limit 100;

-- Add column with tag TPM
alter table macrophage_monocyte_derived_hg19_ctss_bed add column tpm double precision;
update macrophage_monocyte_derived_hg19_ctss_bed set tpm = t1.tpm
    FROM (
        -- Generate a table of ctss_tag and tpm
        select ctss_tag, (tag_count / (select sum(tag_count)::numeric from macrophage_monocyte_derived_hg19_ctss_bed))*10^6 AS tpm
        from macrophage_monocyte_derived_hg19_ctss_bed
        ) AS t1
    where macrophage_monocyte_derived_hg19_ctss_bed.ctss_tag = t1.ctss_tag


--------------------------[ hg19 TSS mapped to pig ]---------------------------

drop table bed_mphage_mcyte_ctss_sscrofa9;
-- This table was obtained by mapping human tags to pig (see labbook 05/01/2011) and converting the SAM output to BED.
select read_table($$ header: ['rname', 'f_start', 'f_end', 'feature', 'mapq', 'strand'], 
                     data_type: ['text', 'int', 'int', 'text', 'int', 'text'], 
                     file:'F:/data/20110105_FANTOM5/mphage_mcyte_f5_hg19_21bp.nordna.single.bed',  
                     table:'fantom5.bed_mphage_mcyte_ctss_sscrofa9' $$) -- 1244124 rows
COMMENT ON TABLE bed_mphage_mcyte_ctss_sscrofa9 IS 'Table was obtained by mapping human tags to pig (see labbook 05/01/2011) and converting the SAM output to BED';

-- Add a column of expression levels as raw counts and TPM found the human samples (NOTE: It might be better to use TPM by scaling by the number of tags mapped to pig rather than the total in humans.)
alter table bed_mphage_mcyte_ctss_sscrofa9 add column tpm double precision;
alter table bed_mphage_mcyte_ctss_sscrofa9 add column tag_count int;

update bed_mphage_mcyte_ctss_sscrofa9 set tpm = macrophage_monocyte_derived_hg19_ctss_bed.tpm, 
                                          tag_count = macrophage_monocyte_derived_hg19_ctss_bed.tag_count
       from macrophage_monocyte_derived_hg19_ctss_bed
       where ctss_tag = feature;
select * from bed_mphage_mcyte_ctss_sscrofa9 limit 10;

-- Produce a table of clustered transcription start sites identified by human TSS mapped to pig.
-- Group human ctss mapping to the same position and same strand. Sum their tpm and total counts (NOTE: that these tpm come from total tags mapped in humans).
-- IMPORTANT: bed file is 0-based coordinate. Here we convert to 1-based.
-- For (+)TSS add 1 to the feature start (leftmost). For (-)TSS add 1 to feature end (rightmost).
-- drop table mphage_mcyte_hg19_on_ss9_ctss;
create table fantom5.mphage_mcyte_hg19_on_ss9_ctss AS(
  -- Forward TSSs:
  select min('TSS_' || rname || '_' || f_start+1 || '_0') AS tss_id,
    rname, f_start + 1 AS ctss_pos, min(strand) AS strand, count(feature) AS no_human_ctss, sum(tpm) AS human_tpm, sum(tag_count) AS human_tag_count
    from bed_mphage_mcyte_ctss_sscrofa9
    where strand = '+'
    group by rname, f_start+1
  -- Reverse TSSs:
  union select min('TSS_' || rname || '_' || f_end+1 || '_16') AS tss_id,
    rname, f_end + 1 AS ctss_pos, min(strand), count(feature) AS no_human_ctss, sum(tpm) AS human_tpm, sum(tag_count) AS human_tag_count
    from bed_mphage_mcyte_ctss_sscrofa9
    where strand = '-'
    group by rname, f_end+1
  ); -- 1235916 rows
comment on table mphage_mcyte_hg19_on_ss9_ctss is 'Clustered TSS obtained by mapping human CTSS, macrophages monocyte derived, to the pig genome Sscrofa 9. Coordinates are 1-based. See 20110105_fantom5.sql' 
select * from mphage_mcyte_hg19_on_ss9_ctss limit 10;

select count(*) from mphage_mcyte_hg19_on_ss9_ctss limit 10;
select sum(human_tpm) from mphage_mcyte_hg19_on_ss9_ctss; -- 613281.290411746 

-- Annotate human-derived TSS with Sus GTF:
-- drop table mphage_mcyte_hg19_on_ss9_ctss_gtf;
create table fantom5.mphage_mcyte_hg19_on_ss9_ctss_gtf AS(
  -- Collect all TSSs that map within (100) bases from the beginning or end of any annotated pig exon.
  select 
       mphage_mcyte_hg19_on_ss9_ctss.tss_id,
       mphage_mcyte_hg19_on_ss9_ctss.rname, 
       mphage_mcyte_hg19_on_ss9_ctss.ctss_pos,
       mphage_mcyte_hg19_on_ss9_ctss.strand,
       mphage_mcyte_hg19_on_ss9_ctss.human_tpm,
       mphage_mcyte_hg19_on_ss9_ctss.human_tag_count,
       attributes,
       feature,
       f_start,
       f_end,
       (CASE WHEN pig_exons.strand = '+' THEN (ctss_pos - f_start) 
             WHEN pig_exons.strand = '-' THEN (f_end - ctss_pos) END)::int AS tss_exon_dist
  from mphage_mcyte_hg19_on_ss9_ctss inner join (select * from sus_scrofa_sscrofa9_59_gtf where feature like 'exon') AS pig_exons ON
    mphage_mcyte_hg19_on_ss9_ctss.rname = pig_exons.rname AND mphage_mcyte_hg19_on_ss9_ctss.strand = pig_exons.strand
  where 
    (pig_exons.strand = '+' and ctss_pos between f_start - 1000 and f_start + 1000) OR
    (pig_exons.strand = '-' and ctss_pos between f_end - 1000 and f_end + 1000)
  ); -- 7716963 ms (128 min)
comment on table mphage_mcyte_hg19_on_ss9_ctss_gtf is 'TSSs from human macrophages FANTOM5 mapping near the 5-prime end of pig exons. See 20110105_fantom5.sql.';
select * from mphage_mcyte_hg19_on_ss9_ctss_gtf limit 10;

select 
    (select count(distinct tss_id) from mphage_mcyte_hg19_on_ss9_ctss) AS no_tss,
    (select count(distinct gtf_attribute(attributes, 'gene_id')) from sus_scrofa_sscrofa9_59_gtf) AS no_ensembl_59_genes,
    (select count(distinct tss_id) from mphage_mcyte_hg19_on_ss9_ctss_gtf) AS no_annotated_tss, 
    (select count(distinct tss_id) from mphage_mcyte_hg19_on_ss9_ctss_gtf where abs(tss_exon_dist) <= 100) AS no_annotated_tss_100bp,
    (select count(distinct gtf_attribute(attributes, 'gene_id')) from mphage_mcyte_hg19_on_ss9_ctss_gtf) AS no_tagged_genes,
    (select count(distinct gtf_attribute(attributes, 'gene_id')) from mphage_mcyte_hg19_on_ss9_ctss_gtf where abs(tss_exon_dist) <= 100) AS no_tagged_genes_100bp;

select sum(human_tag_count) from
    (select distinct tss_id, human_tag_count from mphage_mcyte_hg19_on_ss9_ctss_gtf) AS expr; -- 809129 <- Total expression (tag count) mapped to pig exons
select sum(human_tag_count) from mphage_mcyte_hg19_on_ss9_ctss; -- 3676024 Total expression.

-- Generate transcription clusters. Group TSS separated by no mnore than 50 bp
select cluster_reads($$ select rname, strand, ctss_pos AS tss_start, ctss_pos AS ctss_end, tss_id from mphage_mcyte_hg19_on_ss9_ctss $$, 'public.hg_ss_ctss', -50);

-- bind column of trascription clusters to mphage_mcyte_hg19_on_ss9_ctss:
alter table mphage_mcyte_hg19_on_ss9_ctss add column ctss_id text;
update mphage_mcyte_hg19_on_ss9_ctss set ctss_id = tag_cluster from hg_ss_ctss
  where mphage_mcyte_hg19_on_ss9_ctss.tss_id = hg_ss_ctss.tss_id
drop table public.hg_ss_ctss;
comment on column mphage_mcyte_hg19_on_ss9_ctss.ctss_id is 'ID for transcription start clusters produced by grouping TSS separated by no more than 50 bp. See 20110105_fantom5.sql'

select * from mphage_mcyte_hg19_on_ss9_ctss limit 10;

-- How many genes have been tagged by human TSS tags AND by pig CAGE (050810) tags?
-- TSS identified by both pig CAGE and human TSS on pig:
create table public.common_tss AS (
  select mphage_mcyte_hg19_on_ss9_ctss.tss_id, mphage_mcyte_hg19_on_ss9_ctss.human_tag_count,
                                             cage050810_tss.tss_count
  from mphage_mcyte_hg19_on_ss9_ctss inner join cage050810_tss on 
    mphage_mcyte_hg19_on_ss9_ctss.tss_id = cage050810_tss.tss_id
  );
select * from common_tss;


select distinct gtf_attribute(attributes, 'gene_id') from fantom5.mphage_mcyte_hg19_on_ss9_ctss_gtf
intersect
select distinct gtf_attribute(attributes, 'gene_id') from cage.cage050810_gtf where abs(tss_exon_dist) <= 100;


select * from mphage_mcyte_hg19_on_ss9_ctss_gtf where attributes like '%TNF%'limit 100;
select * from mphage_mcyte_hg19_on_ss9_ctss where ctss_pos between 3868537 - 100 and 3868537 + 100 and rname like '3';
select * from cage050810 where pos between 3868537 - 100 and 3868537 + 100 and rname like '3';
select * from promoters where promoter_start between 3868537 - 1000 and 3868537 + 1000 and rname like '3';


---------------------------------[ TRITUME ]-----------------------------------

select array_agg(mapq)::double precision[] from bed_mphage_mcyte_ctss_sscrofa9 limit 10;

CREATE TABLE sam_mphage_mcyte_ctss_sscrofa9
(
  qname text, -- Query pair NAME if paired; or Query NAME if unpaired
  flag integer, -- Sum of all applicable flags. See SAM format  for details (http://samtools.sourceforge.net/SAM1.pdf)
  rname text, -- Reference sequence NAME
  pos integer, -- 1-based leftmost POSition/coordinate of the clipped sequence
  mapq integer, -- MAPping Quality (phred-scaled posterior probability that the mapping position of this read is incorrect)
  cigar text, -- extended CIGAR string
  mrnm text, -- Mate Reference sequence NaMe; = if the same as <RNAME>
  mpos integer, -- 1-based leftmost Mate POSition of the clipped sequence
  isize integer, -- inferred Insert SIZE
  seq text, -- query SEQuence; = for a match to the reference; n/N/. for ambiguity; cases are not maintained
  qual text, -- query QUALity; ASCII-33 gives the Phred base quality
  tags text -- Optional tags
);



select read_sam($$ file:'F:/data/20110105_FANTOM5/mphage_mcyte_f5_hg19_21bp.nordna.single.sam', dest_table: 'sam_mphage_mcyte_ctss_sscrofa9', overwrite:True $$);
alter table sam_mphage_mcyte_ctss_sscrofa9 set schema fantom5;
select * from sam_mphage_mcyte_ctss_sscrofa9 limit 100;


select * from sam_mphage_mcyte_ctss_sscrofa9 where qname like 'chr14:69865018..69865019,-';




select current_schema();


select distinct rname from macrophage_monocyte_derived_hg19_ctss_bed;
select * from macrophage_monocyte_derived_hg19_ctss_bed where rname like 'chr6' and f_start between 31543344-100 and 31543344+100 -- TNFA; ENST00000449264  chr6:31543344-31546113
select * from macrophage_monocyte_derived_hg19_ctss_bed where rname like 'chr7' and f_start between 5570231-200 and 5570231+200 order by f_start-- ACTB; ENSG00000075624; ENST00000331789 7:5566782-5570340 strand:-1


select read_table($$ file:'D:/Tritume/out_actb.txt', table: 'public.actb', allow_rugged:True, overwrite:True $$)
select read_table($$ file:'D:/Tritume/out_tnf.txt', table: 'public.tnf', allow_rugged:True, overwrite:True $$)
select * from actb where v4 + length(v10)= 5570233;
select * from tnf;
select v4 + length(v10), count(*) from actb group by v4 + length(v10) order by v4 + length(v10);
select v4, count(*) from tnf group by v4 order by v4;
select * from tnf where v4 = 31543338 order by v5;

ACTB tss peak: 'chr7:5570231..5570232,-'


SELECT ctss_tag, rname, CASE WHEN strand = '+' THEN (f_start + 1)
                       WHEN strand = '-' THEN (f_start + 1) - 49 END AS f_start,
                  CASE WHEN strand = '+' THEN (f_start + 1) + 49
                       WHEN strand = '-' THEN (f_start + 1)      END AS f_end,
           strand,
           tpm
    FROM macrophage_monocyte_derived_hg19_ctss_bed
    limit 10



	