-- drop table bowtie_tophat_accepted_hits_sam;
-- select read_sam($$ dataset_id:'CTRL_20100122', file:'C:/Tritume/CTRL_20100122_accepted_hits.sam', dest_table:'bowtie_tophat_accepted_hits_sam' $$);
-- select read_sam($$ dataset_id:'LPS_20100122', file:'C:/Tritume/LPS_20100122_accepted_hits.sam', dest_table:'bowtie_tophat_accepted_hits_sam', append:True $$);

----------------[ Overview of number of mismatches ]---------------------------

-- Pairs with both mates mapped
select dataset_id, edit_distance, count(*) from
(select dataset_id, sam_get_tagvalue(tags, 'NM') AS edit_distance
from bowtie_tophat_accepted_hits_sam 
where (CASE WHEN flag & 4 = 4 THEN False WHEN flag & 8 = 8 THEN False ELSE True END)::boolean = True) AS s1
group by dataset_id, edit_distance;

-- Pairs with one mate only mapped
select dataset_id, edit_distance, count(*) from
(select dataset_id, sam_get_tagvalue(tags, 'NM') AS edit_distance
from bowtie_tophat_accepted_hits_sam 
where (CASE WHEN flag & 4 = 4 THEN False WHEN flag & 8 = 8 THEN False ELSE True END)::boolean = False) AS s1
group by dataset_id, edit_distance;

drop table flag_8;
create temp table flag_8 AS(
  select * from bowtie_tophat_accepted_hits_sam 
  where (flag & 8) = 8
  ); -- 142744 ms.

-- These two queries written out to a file are source for the hist graph on labbook 19/3/2010 
-- Reads with flag 8 on chr 1
select dataset_id, qname, flag, rname, pos
from flag_8
where rname like '1';	
-- Reads with flag 2 on chr 1. To compare with flag 8
select dataset_id, qname, flag, rname, pos, mrnm, mpos 
from bowtie_tophat_accepted_hits_sam 
where flag & 2 = 2 AND rname like '1'; --245287 ms.

-------------------------------------------------------------------------------
-- Num. of reads mapped per dataset: All
-- LPS:6574017; CTRL:5863496
select dataset_id, count(distinct qname) from bowtie_tophat_accepted_hits_sam group by dataset_id;

-- Pairs properly mapped (flag & 2):
create temp table prop_mapped AS(
  select dataset_id, qname, flag, rname, pos, mrnm, mpos 
  from bowtie_tophat_accepted_hits_sam 
  where flag & 2 = 2
  ); -- 290489 ms.
select * from prop_mapped limit 10;
  
select dataset_id, count(*), count(distinct qname) from prop_mapped group by dataset_id; -- 769814 ms.
+--------------+---------+---------+
 CTRL_20100122 | 8431080 | 4215486 |
 LPS_20100122  | 9255954 | 4627913 |
+--------------+---------+---------+

-----------------------------[ Junctions ]-------------------------------------

-- Extract junction datasets from 22/01/2010. Those ran with GFF file
-- Also add junction positions (defined as per TopHat manual)
-- drop table junx;
create temp table junx AS(
  select *, 
    block_start + ((string_to_array(block_sizes, ','))[1])::int -1 as j_left,   -- Left position of the junctions as per tophat *.juncs specs
    block_start + ((string_to_array(block_starts, ','))[2])::int   as j_right   -- Ditto, right
  from tophat 
  where tophat_job like '20100122_CTRL_gff' or 
        tophat_job like '20100122_LPS_gff'
 );
-- Add identifier for each junction:
alter table junx add column junction_id varchar;
update junx set junction_id = chr || '_' || strand || '_' || j_left || '_' || j_right;

select * from junx where tophat_job like '%LPS%' order by chr, block_start limit 10;

/* Add junction_id column to sscrofa_juncs
  alter table sscrofa_juncs add column junction_id varchar;
  update sscrofa_juncs set junction_id = chr || '_' || strand || '_' || j_left || '_' || j_right;
*/

-- All junctions found in LPS, CTRL, Sscrofa 9.56:
-- drop table all_junx;
create temp table all_junx AS(
select distinct 
  chr || '_' || strand || '_' || j_left || '_' || j_right AS junction_id,
  chr, strand, j_left, j_right
  from junx
union select 
  junction_id, chr, strand, j_left, j_right
  from sscrofa_juncs
  ); -- 148852 rows

-- All junctions with lps and ctrl scores (where present) and junctions in pig gff in a sort of cross-tab table
-- drop table junx_count;
create temp table junx_count AS(
  select distinct 
    all_junx.junction_id, 
    ctrl.score::numeric AS ctrl_count,  -- Junction score will be coverted to tpm
    lps.score::numeric  AS lps_count, 
    gff.j AS gff_junction
  from all_junx left join
  -- Junctions in ctrl library
    (select junction_id, score from junx where tophat_job like '20100122_CTRL_gff') as ctrl on
      all_junx.junction_id = ctrl.junction_id
  left join 
  -- Junctions in lps library
    (select junction_id, score from junx where tophat_job like '20100122_LPS_gff') as lps on
      all_junx.junction_id = lps.junction_id
  left join 
  -- Add a tag where junction is documented in gff file (table sscrofa_juncs)
    (select junction_id, 'gff'::varchar AS j from sscrofa_juncs) AS gff on
      all_junx.junction_id = gff.junction_id
  ); -- 148852 rows

-- Export to txt file junctions that are not null in either library, then go to R:
copy (select * from junx_count where (ctrl_count is not null) or (lps_count is not null)) to 'C:/Tritume/junx_count.txt' with csv header delimiter E'\t';


-- Normalized to 1 million if necessay:
update junx_count set ctrl_count = (ctrl_count::numeric/(select sum(ctrl_count) from junx_count))*10^6;
update junx_count set lps_count = (lps_count::numeric/(select sum(lps_count) from junx_count))*10^6;


-- Novel junctions
select 
  (select count(*) from junx_count where gff_junction is null) as all_novel, 
  (select count(*) from junx_count where (gff_junction is null) and (ctrl_count is not null) and (lps_count is null)) as ctrl_novel, -- found in ctrl and not in gff or lps
  (select count(*) from junx_count where (gff_junction is null) and (ctrl_count is null) and (lps_count is not null)) as lps_novel,
  (select count(*) from junx_count where (gff_junction is null) and (ctrl_count is not null) and (lps_count is not null)) as novel_in_both -- Not on gff but present in both libs

-- Differentially used junctions:
select * from tmp_junction_de order by "PValue" limit 100;

---------------------------[ Cufflinks ]---------------------------------------

-- Import .gtf files
-- select read_table($$ file:'C:/Tritume/20100317_RNAseq_CTRL_noGTF/transcripts.gtf.split', skip:1, dest_table:'cufflinks_transcript_gtf', append:True, limit:-1$$)
-- select read_table($$ file:'C:/Tritume/20100317_RNAseq_LPS_noGTF/transcripts.gtf.split', skip:1, dest_table:'cufflinks_transcript_gtf', append:True, limit:-1$$)

/*
Prepare data for diff expr analysis in R. Using cufflinks supplied with GTF 
file and sam file from TopHat 22/01/2010 
*/

-- A sort of cross tab table with [transcript_id; ctrl_fpkm; lps_fpkm]
select trans_id.transcript_id, fpkm_ctrl, fpkm_lps FROM
  -- All trascripts found in selected libs, then FPKM for LPS and CTRL libs.
  (select distinct transcript_id from cufflinks_transcript_gtf where source like '20100317_RNAseq_CTRL' or source like '20100317_RNAseq_LPS') as trans_id
  LEFT JOIN
  (select transcript_id, fpkm AS fpkm_ctrl from cufflinks_transcript_gtf where feature like 'transcript' and source like '20100317_RNAseq_CTRL') as ctrl
    ON trans_id.transcript_id = ctrl.transcript_id 
  LEFT JOIN
  (select transcript_id, fpkm AS fpkm_lps from cufflinks_transcript_gtf where feature like 'transcript' and source like '20100317_RNAseq_LPS') as lps
    ON trans_id.transcript_id = lps.transcript_id;

	select * from tmp_edger_transcript where transcript_id like 'ENSSSCT00000001533'; -- TNF-alpha
select * from tmp_edger_transcript where transcript_id like 'ENSSSCT00000013865'; -- HPRT
select * from tmp_edger_transcript where transcript_id like 'ENSSSCT00000008324'; -- ACTB

-- Send these transcript_id's to biomart and retrieve annotation
select transcript_id from tmp_edger_transcript where "PValue" < 0.01 order by "PValue" ;
-- Load annotation here:
select read_table($$ file:'C:/Downloads/RNAseq_transcript_de_mart_export.txt', dest_table:'tmp_edger_transcript_annotated', header:True $$)

select ensembl_transcript_id, associated_gene_name, associated_transcript_name, array_agg("go_term_name_(bp)"), array_agg("go_term_name_(cc)"), array_agg("go_term_name_(mf)"), description
from tmp_edger_transcript_annotated inner join 
tmp_edger_transcript on transcript_id = ensembl_transcript_id
group by ensembl_transcript_id, associated_gene_name, associated_transcript_name, description, "PValue"
order by "PValue";

(select "go_term_name_(bp)", count(distinct ensembl_transcript_id) AS go_bp_count from tmp_edger_transcript_annotated group by "go_term_name_(bp)" order by go_bp_count desc)
(select "go_term_name_(cc)", count(distinct ensembl_transcript_id) AS go_cc_count from tmp_edger_transcript_annotated group by "go_term_name_(cc)" order by go_cc_count desc)
(select "go_term_name_(mf)", count(distinct ensembl_transcript_id) AS go_mf_count from tmp_edger_transcript_annotated group by "go_term_name_(mf)" order by go_mf_count desc)

-- Number of transcripts in ensembl transcriptome for pig
select count(distinct ensembl_transcript_id) from mart_tome where genome like 'ss9';
-- Number of transcripts from cufflinks (run  with GTF file)
select count(distinct transcript_id) from cufflinks_transcript_gtf 
where feature like 'transcript' and (source like '20100317_RNAseq_LPS' or source like '20100317_RNAseq_CTRL'); 
-- Number of genes and transrcipts
select source, count(distinct transcript_id), count(distinct gene_id) from cufflinks_transcript_gtf group by source;

-- Cufflinks run with GTF file assembles transcripts that are the same size as the the reference transcripts.
-- Check this here:
select ensembl_transcript_id, chr, ss_tome.strand, transcript_start, "start", (transcript_start - "start") AS start_diff, transcript_end, "end", (transcript_end - "end") AS end_diff
  from (select * from mart_tome where genome like 'ss9') AS ss_tome 
  inner join 
  (select * from cufflinks_transcript_gtf where feature like 'transcript') AS cuff_trans on
  ensembl_transcript_id = transcript_id
where (transcript_start - "start")  != 0 or (transcript_end - "end") != 0 -- Are there any cufflinks transcripts different from refernce?
limit 100;

select source, count(distinct transcript_id) AS transcript_count, count(distinct gene_id) AS gene_count from cufflinks_transcript_gtf group by source;


----------------------------------[ Import cuffcompare data ]---------------------

/* 
   Import .refmap and .tmap
   tmap reports for each assembled transcript the closest reference transcript.
   How many reference transcripts are hit by assembled transcripts?
*/
-- select read_table($$ file:'C:/Tritume/RNAseq_LPS_noGTF_transcripts.refmap', dest_table:'cuffcompare_refmap', dataset_id:'RNAseq_LPS_noGTF_transcripts', header:True, limit:-1 $$)
-- select read_table($$ file:'C:/Tritume/RNAseq_CTRL_noGTF_transcripts.refmap', dest_table:'cuffcompare_refmap', dataset_id:'RNAseq_CTRL_noGTF_transcripts', skip:1, append:True, limit:-1 $$)
-- select read_table($$ file:'C:/Tritume/RNAseq_LPS_noGTF_transcripts.tmap', dest_table:'cuffcompare_tmap', dataset_id:'RNAseq_LPS_noGTF_transcripts', header:True, limit:-1 $$)
-- select read_table($$ file:'C:/Tritume/RNAseq_CTRL_noGTF_transcripts.tmap', dest_table:'cuffcompare_tmap', dataset_id:'RNAseq_CTRL_noGTF_transcripts', skip:1, limit:-1, append:True $$)

-- Number of reference transcripts closest to assembled transcripts.
select dataset_id, count(distinct ref_id) AS num_ref_trans, count(distinct cuff_gene_id) AS num_cuff_genes, count(distinct cuff_id) AS num_cuff_trans from cuffcompare_tmap group by dataset_id;
-- Number of reference transcripts hit by assembled transcripts
select dataset_id, count(ref_id) AS num_ref_transcripts, count(distinct ref_id) AS unique_trans from cuffcompare_refmap group by dataset_id;

select distinct cuff_id, ensembl_transcript_id, chr, ref.strand, transcript_start AS ref_start, transcript_end AS ref_end, "start" AS q_start, "end" AS q_end, 
  -- Negative difference: the assembled is shorter than the refernce
  (transcript_start - "start") AS start_diff, ("end" - transcript_end) AS end_diff,
  -- Reference and query length:
  (transcript_end - transcript_start) AS ref_length, ("end" - "start") AS q_length, cuffcompare_tmap.len
from
  (select * from mart_tome where genome like 'ss9') as ref
  inner join cuffcompare_tmap on ensembl_transcript_id = ref_id
  inner join cufflinks_transcript_gtf on transcript_id = cuff_id
where feature like 'transcript'
limit 10;

/* How well are the ensembl transcripts covered by the assembled transcripts? */

-------[ Get the length of the ensembl transcripts (by summing exon lengths) ]-------

-- Get transcript lengths by summing the length of each exon (+1 for each exon) in the GTF annotation file.
-- Important: transcript id is expected to be 18 chars long and starts after: 'transcript_id ' (transcript_id plus backspace). E.g.
-- " gene_id ENSSSCG00000000001; transcript_id ENSSSCT00000000001; exon_number 1; gene_name CELSR1;"
create or replace view view_ss9_56_transcript_length AS(
select substring(attributes from position('transcript_id ' in attributes)+14 for 18) AS transcript_id, 
  sum((f_end - f_start)+1) AS transcript_length 
  from sus_scrofa_sscrofa9_56_gtf 
  where feature like 'exon'
  group by transcript_id
  ); 

-- Total length of ensembl transcripts hit by assembled transcripts 
-- (i.e. length of ensembl trascripts in cuffcompare_refmap)
select dataset_id, sum(transcript_length) AS tot_trans_length, count(transcript_id) AS num_transcripts
from view_ss9_56_transcript_length inner join (select distinct dataset_id, ref_id from cuffcompare_refmap) AS ensembl_trans on
  ref_id = transcript_id
group by dataset_id;

-- Total length of assembled transcripts:
-- Note: Some (not many) transcripts overlap, therefore the total trascriptome length is inflated.
select dataset_id, count(cuff_id), sum(len) 
from
  (select distinct
    cuffcompare_tmap.dataset_id, cuff_id, len
  from cuffcompare_tmap 
    inner join cuffcompare_refmap on
      cuffcompare_tmap.ref_id = cuffcompare_refmap.ref_id
    ) AS refmap_assembled_transcript
group by dataset_id;

select dataset_id, count(cuff_id), count(distinct cuff_id) from cuffcompare_tmap;
  
-- Length of reference transcripts for which the matching assembled transcripts have FPKM >= n
select dataset_id, count(transcript_id) AS num_transcripts, sum(transcript_length) AS tot_length, avg(transcript_length) AS avg_length, stddev(transcript_length) AS stdev_length
from
(select distinct dataset_id, transcript_id, transcript_length
  from cuffcompare_tmap inner join view_ss9_56_transcript_length on
    transcript_id = ref_id
  where "FPKM" >= 5) AS reftrans -- 10089 rows
group by dataset_id; 

-- Length of assembled transcripts where FPKM >= n
select dataset_id, count(distinct ref_id) AS num_ref_trans, count(cuff_id) AS num_assembled_trans, sum(len) AS tot_assembled_trans, avg(len) AS avg_length, stddev(len) AS stdev_length
from 
(select distinct dataset_id, ref_id, cuff_id, len from cuffcompare_tmap where "FPKM" >= 5 and ref_id not like '-%') AS assembled_trans
group by dataset_id;

-- Intersection between reference transcripts and assembled transcripts:
-- Transform reference transcriptome to list of intervals. All exons on each chr in a big list
-- Send to file and go to python for intersection
drop table ref_exon;
create temp table ref_exon AS ( 
select distinct cuffcompare_refmap.dataset_id, gtf.transcript_id,
  'interval.interval(' || array_to_string(array_agg(distinct '[' || f_start || ', ' || f_end || ']'), ', ') || ')' AS exon_boundaries
  from (select 
        substring(attributes from position('transcript_id ' in attributes)+14 for 18) AS transcript_id, * 
        from sus_scrofa_sscrofa9_56_gtf 
        order by transcript_id, f_start) AS gtf
  inner join cuffcompare_refmap on
    gtf.transcript_id = cuffcompare_refmap.ref_id
  inner join cuffcompare_tmap on 
    cuffcompare_tmap.ref_id = cuffcompare_refmap.ref_id and
    cuffcompare_tmap.dataset_id = cuffcompare_refmap.dataset_id
  where feature like 'exon' and 
        cuffcompare_refmap.dataset_id like 'RNAseq_LPS_noGTF_transcripts' and 
        -- Filter for FPKM
        "FPKM" >= 5
  group by cuffcompare_refmap.dataset_id, transcript_id
  order by cuffcompare_refmap.dataset_id, transcript_id
  );
select count(distinct transcript_id) from ref_exon;
  
-- Transform assembled gtf file
drop table assembled_exon;
create temp table assembled_exon AS (
  select distinct dataset_id, ref_id, 'interval.interval(' || array_to_string(array_agg(distinct '[' || "start" || ', ' || "end" || ']'), ', ')  || ')' AS exon_boundaries
  from cuffcompare_tmap inner join cufflinks_transcript_gtf on
    cuffcompare_tmap.cuff_id = cufflinks_transcript_gtf.transcript_id
  where feature like 'exon' and 
        ref_id not like '-%' and 
        dataset_id like 'RNAseq_LPS_noGTF_transcripts' and 
        cuffcompare_tmap."FPKM" >= 5
  group by dataset_id, ref_id
  );
-- Num selected transcripts
select source, count(distinct transcript_id) from cufflinks_transcript_gtf where fpkm >= 5 group by source;

-- Pair reference exons and assembled exons
COPY (
select distinct assembled_exon.*, ref_exon.exon_boundaries AS ref_exon_boundaries
from assembled_exon inner join ref_exon on
  assembled_exon.ref_id = ref_exon.transcript_id
  ) TO 'C:/Tritume/exon_intervals.txt' with csv delimiter E'\t';

-- Comparing some gene of interest: 
-- ACTB. How well covered is it?
select * from sus_scrofa_sscrofa9_56_gtf where attributes like '%ENSSSCT00000008324%' order by f_start; -- ACTB: ENSSSCT00000008324
select * from cuffcompare_refmap where ref_id like 'ENSSSCT00000008324'; -- LPS.115774.0 "CTRL.122412.0"
select * from cufflinks_transcript_gtf where transcript_id like 'LPS.115774.0' order by "start";

-- GAPDH
select * from sus_scrofa_sscrofa9_56_gtf where attributes like '%ENSSSCT00000000756%' order by f_start; -- GAPDH: ENSSSCT00000000756
select * from cuffcompare_refmap where ref_id like 'ENSSSCT00000000756'; 
select * from cufflinks_transcript_gtf where transcript_id like 'LPS.152114.0' or transcript_id like 'CTRL.160181.0' order by source, start;

select * from cuffcompare_refmap where ref_id like 'ENSSSCT00000008324'; -- LPS.115774.0 
select * from cuffcompare_tmap where ref_id like 'ENSSSCT00000008325';

-- TNFA_PIG
select * from sus_scrofa_sscrofa9_56_gtf where attributes like '%ENSSSCT00000001533%' order by f_start;
select * from cuffcompare_refmap where ref_id like 'ENSSSCT00000001533'; 
select * from cufflinks_transcript_gtf where transcript_id in ('LPS.175619.0', 'CTRL.184978.0', 'CTRL.184981.0') order by source, start;

select * from cuffcompare_refmap where ref_id like 'ENSSSCT00000008324'; -- LPS.115774.0 
select * from cuffcompare_tmap where ref_id like 'ENSSSCT00000008325';
"LPS.175619|LPS.175619.0"
"CTRL.184978|CTRL.184978.0,CTRL.184981|CTRL.184981.0"

select * from sus_scrofa_sscrofa9_56_gtf where attributes like '%ENSSSCT00000008324%' and feature like 'exon' order by f_start; -- ACTB: ENSSSCT00000008324, HPRT1 ENSSSCT00000013865
select * from cufflinks_transcript_gtf where transcript_id like 'ENSSSCT00000013865' limit 100;
select * from cuffcompare_refmap where ref_id like 'ENSSSCT00000008324'; -- LPS.115774.0 
select * from cuffcompare_tmap where ref_id like 'ENSSSCT00000013865';
select dataset_id, count(distinct ref_id) from cuffcompare_tmap where ref_id not like '-%' group by dataset_id;

------------------------[ Allele specific expression (SNP calling) ]-----------

select read_table($$ file:'C:/Tritume/20100406_RNAseq_LPS_sscrofa9.56.pileup.snp', dest_table:'tmp_pileup_snp', limit:-1, overwrite:True, 
  header: True, data_type:'varchar, varchar, int, varchar, varchar, int, int, int, int, varchar, int, varchar, int' $$);
select read_table($$ file:'C:/Tritume/20100406_RNAseq_CTRL_sscrofa9.56.pileup.snp', dest_table:'tmp_pileup_snp', limit:-1, append:True, 
  header: True $$);

select * from tmp_pileup_snp limit 1000;
select dataset_id, count(dataset_id) from tmp_pileup_snp group by dataset_id;

-- SNPs in common btw LPS and CTRL
-- drop table snp;
create temp table snp AS(
  select 
    lps.rname || '_' || lps.pos AS snp_id, lps.rname, lps.pos, lps.allele_1 AS lps_1, lps.count_allele_1 AS lps_count_1, lps.allele_2 AS lps_2, lps.count_allele_2 AS lps_count_2,
    ctrl.allele_1 AS ctrl_1, ctrl.count_allele_1 AS ctrl_count_1, ctrl.allele_2 AS ctrl_2, ctrl.count_allele_2 AS ctrl_count_2,
    log(2, lps.count_allele_1::numeric/lps.count_allele_2) AS lps_logfold, 
    log(2, ctrl.count_allele_1::numeric/ctrl.count_allele_2) AS ctrl_logfold,
    log(2, (lps.count_allele_1::numeric / (lps.count_allele_1 + lps.count_allele_2))/
           (ctrl.count_allele_1::numeric/ (ctrl.count_allele_1 + ctrl.count_allele_2))
           ) AS allele_ratio
  from (select * from tmp_pileup_snp where dataset_id like '20100406_RNAseq_LPS') as lps inner join 
      (select * from tmp_pileup_snp where dataset_id like '20100406_RNAseq_CTRL') as ctrl on
      lps.rname= ctrl.rname and lps.pos= ctrl.pos
  -- Remove snps where the genotype of ctrl and lps disagree
  where (lps.allele_1 like ctrl.allele_1) or (lps.allele_2 like ctrl.allele_2) -- SNPs removed (but these could be interesting if confirmed?!)
  order by allele_ratio, lps.rname, lps.pos
  );
select count(distinct snp_id) from snp;
select * from snp limit 10;
------------------[ Are there SNP in transcript of interest (e.g. differentailly expressed?) ]--------------

-- Transcripts passed to edgeR analysis and their exon locations
-- drop table edger;
create temp table edger AS (
  select distinct transcript_id, rname, feature, f_start, f_end, strand, attributes
  from tmp_edger_transcript inner join sus_scrofa_sscrofa9_56_gtf on
  transcript_id = substring(attributes from position('transcript_id ' in attributes)+14 for 18)
  where feature like 'exon' -- "PValue" < 0.01 
  );
select * from edger limit 20;
select count(*) from edger;
select count(distinct transcript_id) from edger;

-- SNPs in the exon extracted above
-- drop table snp_transcripts;
create temp table snp_transcripts AS (
  select distinct transcript_id, snp_id, snp.rname, pos, f_start, f_end, strand, allele_ratio
  from snp inner join edger on
    snp.rname = edger.rname and
    snp.pos between edger.f_start and edger.f_end
    where transcript_id in (select transcript_id from tmp_edger_transcript where "PValue" < 0.01 order by "logFC" desc limit 20) -- (select transcript_id from tmp_edger_transcript order by random() desc limit 20) -- 
  );
select * from snp_transcripts order by allele_ratio;
select transcript_id, count(*) AS n_snp from snp_transcripts group by transcript_id order by n_snp desc limit 100;

/* 
Note: table tmp_edger_transcript_annotated obtained by feeding biomart with the transcript_id's from 
tmp_edger_transcript. 
Attributes: 
['ensembl_transcript_id', 'associated_gene_name', 'associated_transcript_name', 'transcript_count', 'go_term_name_(bp)', 'go_term_name_(cc)', 'go_term_name_(mf)', 'description']
*/

-- Annotate top n transcripts with largest fold change and pvalue < 0.01
select distinct transcript_id, "logFC", "PValue", associated_gene_name, array_to_string(array_agg("go_term_name_(mf)"), ', '), description
from 
  (select * from tmp_edger_transcript where "PValue" < 0.01) as edger
inner join 
  (select distinct ensembl_transcript_id, associated_gene_name, "go_term_name_(mf)", description from tmp_edger_transcript_annotated) as ensembl on
transcript_id = ensembl_transcript_id
where description is not null
group by transcript_id, associated_gene_name, description, "logFC", "PValue"
order by "logFC" desc, "PValue" asc
limit 100;

---------------------------[ Notes ]----------------------------------------

select * from cufflinks_transcript_gtf where feature like 'transcript';	
select distinct source from cufflinks_transcript_gtf;

select distinct
    cuffcompare_tmap.dataset_id, cuffcompare_tmap.ref_id, cuff_id, len, "start", "end", feature, transcript_start, transcript_end 
from cuffcompare_tmap 
  inner join cuffcompare_refmap on
    cuffcompare_tmap.ref_id = cuffcompare_refmap.ref_id
  inner join cufflinks_transcript_gtf on
    cufflinks_transcript_gtf.transcript_id = cuffcompare_tmap.cuff_id
  inner join mart_tome on
    ensembl_transcript_id = cuffcompare_refmap.ref_id
where feature like 'transcript'
order by cuffcompare_tmap.dataset_id, cuffcompare_tmap.ref_id, "start"
limit 1000;


select  get_column_names('cuffcompare_refmap');
"['dataset_id', 'ref_gene_id', 'ref_id', 'class_code', 'cuff_id_list']"

select  get_column_names('mart_tome');
"['genome', 'ensembl_gene_id', 'ensembl_transcript_id', 'chr', 'transcript_start', 'transcript_end', 'strand']"

select get_column_names('cuffcompare_tmap');
"['dataset_id', 'ref_gene_id', 'ref_id', 'class_code', 'cuff_gene_id', 'cuff_id', 'FMI', 'FPKM', 'FPKM_conf_hi', 'FPKM_conf_lo', 'cov', 'len', 'major_iso_id']"

select get_column_names('cufflinks_transcript_gtf');
"['rname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'gene_id', 'transcript_id', 'exon_number', 'fpkm', 'frac', 'conf_lo', 'conf_hi', 'cov']"

select * from cuffcompare_refmap limit 10;
select * from cuffcompare_tmap limit 10;
select * from cufflinks_transcript_gtf where transcript_id like 'CTRL.100034.0';

select * from mart_tome where ensembl_gene_id like 'ENSSSCG00000000005'; "5";126757;134285;
select * from cufflinks_transcript_gtf where transcript_id like 'LPS.144798.0'; "5";133990;134251

select * from mart_tome where ensembl_transcript_id like 'ENSSSCT00000004429';

---------------------------[ Tritume? ]----------------------------------------

select read_table($$ file:'C:/Tritume/20100406_RNAseq_LPS_sscrofa9.56.pileup.snp', dest_table:'tmp_pileup_snp', limit:-1, overwrite:True, 
  header: True, data_type:'varchar, varchar, int, varchar, varchar, int, int, varchar, int, varchar, int' $$)
select read_table($$ file:'C:/Tritume/20100406_RNAseq_CTRL_sscrofa9.56.pileup.snp', dest_table:'tmp_pileup_snp', limit:-1, append:True, 
  header: True $$)


-- select set_column_names('sus_scrofa_sscrofa9_56_gtf', $$['rname', 'source', 'feature', 'f_start', 'f_end', 'score', 'strand', 'frame', 'attributes']$$)


select distinct v2 from sus_scrofa_sscrofa9_56_gtf;

select source, count(distinct transcript_id) from cufflinks_transcript_gtf group by source;

" gene_id ENSSSCG00000004698; transcript_id ENSSSCT00000005188; exon_number 1; gene_name SERINC4;"

select * from cufflinks_transcript_gtf limit 10;

ALTER TABLE "Pigs".tmp_edger_transcript ADD CONSTRAINT pk_edger_transcript PRIMARY KEY (transcript_id);

select array_agg(job_id)
from bowtie_jobs

select array_agg(column_name::text)
from information_schema.columns



"ENSSSCT00000001533";-14.5449422799917;3.89426540450006;5.24923311012038e-016;5.22088725132573e-014

drop table tmp_edger_transcript ;
alter table tmp_junction_edger  rename to tmp_edger_junction;

drop index ind_junx;
alter table all_junx drop constraint pk_junction_id;
-- add primary key
alter table all_junx add constraint pk_junction_id primary key (junction_id);
-- add index
create index ind_junx on junx (junction_id);

select 10^6;

-- Import .juncs files containing junction positions.
-- select read_table($$ file:'C:/Tritume/20100122_CTRL_junctions.juncs', dest_table:'tophat_juncs', dataset_id:'CTRL_20100122', header:['dataset_id', 'chr', 'j_left', 'j_right', 'strand'], data_type:'varchar, varchar, int, int, varchar' $$);
-- select read_table($$ append:True, file:'C:/Tritume/20100122_LPS_junctions.juncs', dest_table:'tophat_juncs', dataset_id:'LPS_20100122', header:['dataset_id', 'chr', 'j_left', 'j_right', 'strand'], data_type:'varchar, varchar, int, int, varchar' $$);

select * from tophat_juncs limit 10;

select * from tmp_20100122_rnaseq_lps_coverage_wig order by rname, start limit 5000;

create temp table ctrl_hits AS (
  select qname, flag, rname, pos, mrnm, mpos from bowtie_tophat_accepted_hits_sam where dataset_id like 'CTRL_20100122'
)
select * from ctrl_hits limit 10;

create temp table lps_hits AS (
  select qname, flag, rname, pos, mrnm, mpos from bowtie_tophat_accepted_hits_sam where dataset_id like 'LPS_20100122'
)
select * from lps_hits limit 10;

select count(distinct qname) from ctrl_hits; -- 5863496 (311837 ms)



select * from bowtie_tophat_accepted_hits_sam limit 10;

-- drop table accepted_hits;
create temp table accepted_hits AS(
select count(*)
from bowtie_tophat_accepted_hits_sam
where (CASE WHEN flag & 4 = 4 THEN False WHEN flag & 8 = 8 THEN False ELSE True END)::boolean = False
group by qname
);
select count(count), sum(count) from accepted_hits limit 10;
select count, count(*) AS num_hits from accepted_hits group by count order by count;


-- Reject reads whose mate is not mapped
create temp table accepted_hits AS (
select dataset_id, qname, rname, pos, mrnm, mpos
from bowtie_tophat_accepted_hits_sam 
where (CASE WHEN flag & 4 = 4 THEN False WHEN flag & 8 = 8 THEN False ELSE True END)::boolean = True
);

select dataset_id, count(distinct qname) from accepted_hits group by dataset_id;



create temp table accepted_hits AS(
select count(*)
from bowtie_tophat_accepted_hits_sam
where (CASE WHEN flag & 4 = 4 THEN False WHEN flag & 8 = 8 THEN False ELSE True END)::boolean = True
group by qname
);

select *
from ctrl_accepted_hits
where count = 1
order by count;

select * from tmp_ctrl_accepted_hits_sam where qname like 'EBRI093151_0001:8:15:242:1330#NNNNNN';
select * from tmp_ctrl_accepted_hits_sam  limit 10;

select read_table($$ file:'C:/Tritume/20100122_RNAseq_LPS_coverage.wig', dest_table:'tmp_20100122_rnaseq_lps_coverage_wig', limit:-1, skip:1,
    dataset_id: '20100122_RNAseq_LPS_coverage', data_type:'varchar, varchar, integer, integer, integer', header: ['dataset_id', 'rname', 'start', 'end', 'value'] $$);
drop table tmp_20100122_rnaseq_lps_coverage_wig;

select * from tmp_20100122_rnaseq_lps_coverage_wig limit 10; --8316034 rows

select (select 'tmp_20100122_rnaseq_lps_coverage_wig') from tmp_20100122_rnaseq_lps_coverage_wig;


select * from accepted_hits_sam;

select flag, count(*) AS flag_count
from tmp_ctrl_accepted_hits_sam
group by flag
;

select 129 & 2;


select * from bowtie_tophat_accepted_hits_sam where flag & 8 = 8 limit 10;

drop table tmp_sam;
select read_table($$file:'C:/Tritume/LPS_20100122_accepted_hits.sam', dest_table:'tmp_sam', limit:100000, overwrite:True, allow_rugged:True$$);
select * from tmp_sam order by v1 limit 1000;

select v1, count(*) from tmp_sam 
group by v1 
having count(*) > 2
order by count desc;
select * from tmp_sam where v1 in (select v1 from tmp_sam group by v1 order by count(*) desc limit 100) order by v1, v2;


select * from tmp_20100122_rnaseq_lps_coverage_wig 
where rname like '3' and ("start" between 38956245 - 10 and 38956245 + 10) or ("end" between 38957024 - 10 and 38957024 + 10)
limit 100;

select * from cufflinks_transcript_gtf where transcript_id like 'ENSSSCT00000006150';
drop table tmp_edger_transcript;
where 3_-_38956245_38957024
