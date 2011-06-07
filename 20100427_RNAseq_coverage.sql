/*
   Explore read coverage from RNAseq across known transcripts extracted from 
   GTF file.

   Pileup dataset comes from 26/4/2010 (RNAseq CTRL and LPS combined). Tophat
   run from 22/1/2010 (selected reads aligned with the aid of GTF file)
   
*/

------------------------[ Import pileup-genotype file ]------------------------

-- Table describing dataset ID(s)
select read_table($$ file:'C:/Tritume/pileup_genotype.txt', dest_table: 'pileup_genotype_dataset', header:True, limit: -1, data_type: 'varchar, date, varchar, text, text' $$)
select index('pileup_genotype_dataset', 'dataset_id', 'p');
comment on table pileup_genotype_dataset is 'Description of datasets in table pileup_genotype';
select * from pileup_genotype_dataset;

-- Actual data
-- Tophat output from CTRL and LPS concatenated and sent to 'samtools pileup'
select read_table($$ file:'C:/Tritume/20100426_RNAseq_LPS-CTRL_conc_sscrofa9.56.pileup.snp', dest_table: 'pileup_genotype', header:True,
  limit: -1, overwrite: False,
  data_type: 'varchar, varchar, int, varchar, varchar, int, int, int, int, varchar, int, varchar, int'
  $$)
select read_table($$ file:'C:/Tritume/	', dest_table: 'pileup_genotype', header:True,
  limit: -1, overwrite: False,
  data_type: 'varchar, varchar, int, varchar, varchar, int, int, int, int, varchar, int, varchar, int'
  $$)

Comment on table pileup_genotype is 'Output from "samtools pileup -s -c -f" converted to genotypes by "yyyymmdd_pileup_to_genotypes.py"'

--


-- drop table tmp_ref_trans;
create table tmp_ref_trans AS(
    -- Format annotation in GTF file
    select 
      array_to_string(regexp_matches(attributes, '.*transcript_id (.*?);.*'), '|') AS transcript_id,
      array_to_string(regexp_matches(attributes, '.*exon_number(.*?);.*'), '|')::int AS exon_number,
      array_to_string(regexp_matches(attributes, '.*gene_name (.*?);.*'), '|') AS gene_name,      
      rname, f_start, f_end,
      0::integer AS nreads, 
      (f_end - f_start) AS feature_length 
    from sus_scrofa_sscrofa9_56_gtf 
    where feature like 'exon'
    -- Uncomment any of these lines to select only some transcripts
    AND (select array_to_string(regexp_matches(attributes, '.*gene_name (.*?);.*'), '|')) in ('TNFA_PIG', 'ACTB_PIG', 'HPRT_PIG', 'IL8_PIG', 'IL23A_PIG')
    -- AND (select array_to_string(regexp_matches(attributes, '.*transcript_id (.*?);.*'), '|')) in (select transcript_id from tmp_edger_transcript)
    order by transcript_id, exon_number, f_start
    );
/*
select index('tmp_ref_trans', 'rname, f_start, f_end', 'i');
select index('tmp_ref_trans', 'rname', 'i');
select index('tmp_ref_trans', 'f_start', 'i');
select index('tmp_ref_trans', 'f_end', 'i');
*/

select * from tmp_ref_trans limit 100;

-- drop table tmp_coverage;
create table tmp_coverage AS(
-- Link coverage to annotation
select transcript_id, gene_name, exon_number, rname, pos, max(nreads) AS coverage, feature_length FROM (
      -- Annotate position effectively sequenced
      select tmp_ref_trans.transcript_id, tmp_ref_trans.gene_name, exon_number, tmp_ref_trans.rname, pos, pileup_genotype.nreads, feature_length 
      from pileup_genotype inner join tmp_ref_trans on
        tmp_ref_trans.rname = pileup_genotype.rname and
        pileup_genotype.pos >= tmp_ref_trans.f_start and
        pileup_genotype.pos <= tmp_ref_trans.f_end
        -- Make sure you are using the right dataset!
        where dataset_id like '20100426_RNAseq_conc'
     -- Add exon starts and ends from GTF
     union select transcript_id, gene_name, exon_number, rname, f_start AS pos, nreads, feature_length from tmp_ref_trans
     union select transcript_id, gene_name, exon_number, rname, f_end AS pos, nreads, feature_length from tmp_ref_trans
   ) AS t1
  -- Group by to remove start and end positions ('union select') that have been sequenced
  group by transcript_id, gene_name, exon_number, rname, pos, feature_length
  order by transcript_id, exon_number, pos
    );
select * from tmp_coverage limit 100;

-- drop table tmp_exons;
create table tmp_exons AS(
-- Format the coverage table for plotting in R
select t1.*, 
  -- Cumulative exon length. I.e. as if all exon are stiched together starting from position 1. 
  -- CASE ... is to make the first exon cumulate from 1 (WHEN clause) and to add a spacer between exons (ELSE clause)
  sum(CASE WHEN tmp_ref_trans.feature_length is null THEN 1 ELSE tmp_ref_trans.feature_length + 0 END) + exon_pos AS exon_offset
  from tmp_ref_trans right join 
      (
       -- Subquery t1: Add a column for the exon position. I.e. make the position of each exon to restart from 0
       select tmp_coverage.*, (pos - f_start) AS exon_pos
       from tmp_coverage 
       inner join tmp_ref_trans 
       -- Link transcript to transcript
       on tmp_coverage.transcript_id = tmp_ref_trans.transcript_id and
       -- Link exon number ot exon number
       tmp_coverage.exon_number = tmp_ref_trans.exon_number
      order by tmp_coverage.transcript_id, tmp_coverage.exon_number, pos
      ) AS t1 
  -- Link transcript in pileup to reference transcript
  on t1.transcript_id = tmp_ref_trans.transcript_id and
     -- To make cumulative, link to previuos exon number
     t1.exon_number >= (tmp_ref_trans.exon_number+1) 
group by t1.transcript_id, t1.gene_name, t1.exon_number, t1.rname, t1.pos, t1.coverage, t1.feature_length, t1.exon_pos
order by t1.transcript_id, t1.exon_number, t1.pos
);

select * from tmp_exons limit 100;

drop table tmp_ref_trans;
drop table tmp_coverage;
drop table tmp_exons;

----------------------------------[ Tritume ]----------------------------------
CASE WHEN (exon_offset + tmp_ref_trans.feature_length) is null THEN exon_offset ELSE  (exon_offset + tmp_ref_trans.feature_length) + 1 END AS cum_exon
select names($$ t:'tmp_coverage', quote:''  $$) 
select * from describe('ref_trans');
"dataset_id", "rname", "pos", "rbase", "consensus", "consensus_phred", "snp_phred", "rms", "nreads", "allele_1", "count_allele_1", "allele_2", "count_allele_2"

select * from ensembl_annotation_sscrofa9_57 where associated_gene_name in ('TNFA_PIG', 'ACTB', 'HPRT_PIG', 'IL8', 'IL23A');
select count(*) from pileup_genotype where nreads < 6;

select tmp_ref_trans.transcript_id, tmp_ref_trans.gene_name, exon_number, tmp_ref_trans.rname, pos, pileup_genotype.nreads, feature_length 
      from pileup_genotype inner join tmp_ref_trans on
        tmp_ref_trans.rname = pileup_genotype.rname and
        pileup_genotype.pos >= tmp_ref_trans.f_start and
        pileup_genotype.pos <= tmp_ref_trans.f_end
limit 3100;

drop index "index_tmp_ref_trans.f_end";

select distinct * from (select *, array_to_string(regexp_matches(attributes, '.*gene_name (.*?);.*'), '|') as gene from sus_scrofa_sscrofa9_56_gtf ) as t1
where gene in ('TNFA_PIG', 'ACTB_PIG', 'HPRT_PIG', 'IL8_PIG', 'IL23A_PIG')

select * from sus_scrofa_sscrofa9_56_gtf  order by attributes, rname, f_start limit 100 ;
