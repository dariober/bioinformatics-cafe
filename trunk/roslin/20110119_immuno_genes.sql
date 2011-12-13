show search_path;
-------------------------------------------------------------------------------
-- Attributes of the genes/transcripts of interest
-------------------------------------------------------------------------------
-- This data produced by R/biomaRt see 20110119_immuno_genes_biomart.R

-- TRANSCRIPTS
select read_table($$ file:'D:/Tritume/transcript_attributes.txt', table: 'immunogenes.transcript_attributes', header:True, overwrite:True $$);
alter table transcript_attributes alter column chromosome_name type text;

-- Add column of with the gene 5'end useful for later (this is the predicted TSS):
alter table transcript_attributes add column transcript_5p_end int;
update transcript_attributes set transcript_5p_end = (CASE WHEN strand = '+' THEN transcript_start ELSE transcript_end END);

create index ind_transcr_attr on transcript_attributes ( transcript_5p_end, chromosome_name, species, strand );
vacuum analyze transcript_attributes;

select * from transcript_attributes;

/* Deprecated use transcript_attributes instead:
 GENES (the info here is duplicated in transcript above)
select read_table($$ file:'D:/Tritume/gene_attributes.txt', table: 'immunogenes.gene_attributes', header:True, overwrite:True $$);
alter table gene_attributes alter column chromosome_name type text;

-- Add column of with the gene 5'end useful for later (this is the predicted TSS):
alter table gene_attributes add column gene_5p_end int;
update gene_attributes set gene_5p_end = (CASE WHEN strand = '+' THEN start_position ELSE end_position END);

select * from gene_attributes;
*/

-------------------------------------------------------------------------------
-- CAGE TSSs in three species, coordinates 1-based and chromosomes w/o 'chr'
-------------------------------------------------------------------------------
drop table immunogenes.tss_bed;
create table immunogenes.tss_bed AS (
    SELECT -- H. sapiens
      regexp_replace(rname, '^chr', '') AS chr,
      CASE WHEN strand like '+' THEN f_start + 1 ELSE f_end + 1 END AS tss_pos,
      ctss_tag AS tss_id,
      strand AS strand,
      tpm AS tpm,
      'hsapiens'::text AS species,
      'macrophage_monocyte_derived_hg19'::text AS rna
      FROM macrophage_monocyte_derived_hg19_ctss_bed

    UNION SELECT -- M. musculus
      regexp_replace(chrom, '^chr', '') AS chr,
      tss_pos + 1 AS tss_pos, -- Make 1-based coordinates
      tss_id AS tss_id,
      strand AS strand,
      tpm_ex_ribo AS tpm,
      'mmusculus' AS species,
      'CGO macrophage unstimulated sample' AS rna
      FROM "FANTOM4".mouse_macrophage_tss_bed 
      WHERE rna like 'CGO macrophage unstimulated sample' -- Use this library because is the largest unstimulated (1542181 tags).
      
    UNION SELECT -- S. scrofa
      rname,
      tss_pos,
      tss_id,
      strand,
      tss_tpm AS tpm,
      'sscrofa' AS species,
      'Alveolar macrophages unstimulated' AS rna
      FROM cage050810_tss
      );
comment on table tss_bed is 'Position and quantification of TSS in different species. Coordinates are 1-based. See 20110119_immuno_genes.sql';
create index ind_tss_species on tss_bed ( tss_pos, chr, species, strand );

select * from tss_bed limit 100;
vacuum analyze tss_bed;

-------------------------------------------------------------------------------
-- Fetch TSSs near transcripts of interest
-------------------------------------------------------------------------------

-- drop table tss_immunogenes ;
create table tss_immunogenes AS (
    select distinct 
           tss_bed.species, 
           reference_gene_name, 
           ensembl_gene_id,
           ensembl_transcript_id, 
           external_transcript_id, 
           tss_bed.chr, 
           transcript_start, 
           transcript_end, 
           tss_bed.strand,
           tss_bed.tss_id, 
           tss_bed.tss_pos,
           transcript_5p_end,
           CASE WHEN tss_bed.strand = '+' THEN (tss_pos - transcript_5p_end) ELSE (transcript_5p_end - tss_pos) END AS tss_transcript_dist,  -- Distance between CAGE tss and ensembl transcript 5'end. (-)ve means TSS is upstream.
           tpm,
           dense_rank() OVER (PARTITION BY ensembl_transcript_id ORDER BY tpm DESC)::int AS tpm_rank_transcript, -- TPM rank of each TSS within transcripts. 1 for the most expressed TSS.
           dense_rank() OVER (PARTITION BY ensembl_gene_id ORDER BY tpm DESC, abs(tss_pos - transcript_5p_end) ASC)::int AS tpm_rank_gene, -- TPM rank of each TSS within genes. 1 for the most expressed TSS.
           rna
    from tss_bed inner join transcript_attributes on
        tss_bed.species = transcript_attributes.species and
        tss_bed.chr = transcript_attributes.chromosome_name and
        transcript_attributes.strand = tss_bed.strand
    where tss_pos between transcript_5p_end - 1000 and transcript_5p_end + 1000
    order by reference_gene_name, tss_bed.species, tss_transcript_dist
    );
comment on table tss_immunogenes is 'Collect all TSSs in the vicinity of the 5-end of each transcript of interest. Use "SELECT * FROM tss_immunogenes WHERE tpm_rank_gene = 1;" to have the TSS peak for each gene. See 20110119_immuno_genes.sql';
comment on column tss_immunogenes.tpm_rank_transcript is 'TPM rank of each TSS within transcripts. 1 for the most expressed TSS.';
comment on column tss_immunogenes.tpm_rank_gene is 'TPM rank of each TSS within genes. Ties (same peak shared by multiple transcripts) are separated by the distance between transcript 5end and peak. 1 for the most expressed TSS and nearest transcript.';
comment on column tss_immunogenes.tss_transcript_dist is 'Distance between CAGE tss and ensembl transcript 5-end as (tss - transcr) if (+)strand or (trans - tss) if (-)strand. (-)ve means TSS is upstream.';
select * from tss_immunogenes order by reference_gene_name, species, ensembl_transcript_id, tss_transcript_dist;
select * from tss_immunogenes where reference_gene_name like 'ACTB' and species like 'hsapiens' order by tss_pos, tpm_rank_gene;


create temp table tss_peaks  AS (select * from tss_immunogenes where tpm_rank_gene = 1 order by reference_gene_name, species); -- select * from tss_peaks;
drop table tss_peaks_immunogenes;
create table tss_peaks_immunogenes AS(
select distinct 
    tss_peaks.species,
    tss_peaks.reference_gene_name,
    tss_peaks.ensembl_gene_id,
    tss_peaks.ensembl_transcript_id AS reference_transcript, -- The closest transcript to the peak
    tss_peaks.chr,
    tss_peaks.strand,
    tss_peaks.tss_id tss_peak_id,
    tss_peaks.tss_pos AS tss_peak_pos,
    tss_peaks.transcript_5p_end AS ref_trans_5p_end,
    tss_bed.tss_pos, 
    tss_bed.tpm,
    CASE WHEN tss_bed.strand = '+' THEN (tss_bed.tss_pos - tss_peaks.transcript_5p_end) ELSE (tss_peaks.transcript_5p_end - tss_bed.tss_pos) END::int AS tss_reftrans_dist -- Distance between TSS and reference transcript. (-)ve means TSS is upstream
from tss_peaks right join tss_bed on
    tss_peaks.chr = tss_bed.chr and tss_peaks.strand = tss_bed.strand and tss_peaks.species = tss_bed.species
where tss_bed.tss_pos between tss_peaks.tss_pos - 1000 and tss_peaks.tss_pos + 1000
order by tss_peaks.reference_gene_name, tss_peaks.species, tss_peaks.chr, tss_bed.tss_pos
);
comment on table tss_peaks_immunogenes is 'Collect TSS around the highest TSS peak of each gene and also report the distance between peak and closest transcript.';



select tss_immunogenes.*, (tss_immunogenes.tss_pos - tss_peaks.transcript_5p_end) AS dist_ref_transcript
from tss_immunogenes inner join tss_peaks on tss_immunogenes.ensembl_gene_id = tss_peaks.ensembl_gene_id
where (tss_immunogenes.tss_pos - tss_peaks.transcript_5p_end) between -100 and 100  and tss_immunogenes.reference_gene_name = 'ACTB' and tss_immunogenes.species = 'hsapiens'
order by tss_immunogenes.reference_gene_name, dist_ref_transcript, tss_immunogenes.tss_pos;

select * from tss_immunogenes where tpm_rank_gene = 1 order by reference_gene_name, species;


-------------------------------------------------------------------------------
-- Representative TSS peaks
-------------------------------------------------------------------------------

-- Fetch peaks
create table tss_genes AS(
select distinct species,
       reference_gene_name,
       ensembl_gene_id,
       chr,
       strand,
       tss_id,
       tss_pos,
       tpm,
       tss_transcript_dist,
       transcript_5p_end,
       tpm_rank_gene
from tss_immunogenes where tpm_rank_gene = 1 order by reference_gene_name, species
);
comment on table tss_genes is 'Putative transcription start site of each gene. That is, the TSS with the highest tpm in the region +/- 1000 bp from any 5p-end of each transcript within gene.'

-------------------------------------------------------------------------------
-- SEQUENCES from biomaRt
-------------------------------------------------------------------------------

-- Import sequences 300 bp upstream of coding region 	
/* SELECT read_table($$ file: 'D:/Tritume/coding_gene_flank_up300_hs.txt', table: 'coding_gene_flank', header:True, overwrite:True, CSV:False $$);
SELECT read_table($$ file: 'D:/Tritume/coding_gene_flank_up300_mm.txt', table: 'coding_gene_flank', header:True, append:True, CSV:False $$);
SELECT read_table($$ file: 'D:/Tritume/coding_gene_flank_up300_ss.txt', table: 'coding_gene_flank', header:True, append:True, CSV:False $$); */
comment on table coding_gene_flank is '300 bp flanking upstream the coding sequence. Sequence is reverse complementet when on the (-)strand. Extracted with biomaRt setting getSequence(... seqType= "coding_gene_flank" ). 20110119_immuno_genes-2.R ';
select * from coding_gene_flank;

-- Import full gene sequence +300 bp upstream
/* SELECT read_table($$ file: 'D:/Tritume/gene_exon_intron_up300_hs.txt', table: 'gene_exon_intron', header:True, overwrite:True, CSV:False $$);
SELECT read_table($$ file: 'D:/Tritume/gene_exon_intron_up300_mm.txt', table: 'gene_exon_intron', header:True, append:True, CSV:False $$);
SELECT read_table($$ file: 'D:/Tritume/gene_exon_intron_up300_ss.txt', table: 'gene_exon_intron', header:True, append:True, CSV:False $$); */
comment on table gene_exon_intron is '300 bp flanking upstream the coding sequence. Sequence is reverse complementet when on the (-)strand. Extracted with biomaRt setting getSequence(... seqType= "gene_exon_intron". See 20110119_immuno_genes-2.R ) ';
select *, length(gene_exon_intron) from gene_exon_intron order by species, ensembl_gene_id;

-- Region flanking the gene 5'end +300/-100 bp. Let's call this part "promoter region":
select 
    gene_attributes.*,
    CASE WHEN strand = 1  THEN start_position -300 
         WHEN strand = -1 THEN end_position   -100 END AS promoter_start, -- leftmost position, regardless of strand
    CASE WHEN strand = 1  THEN start_position +100 
         WHEN strand = -1 THEN end_position   +300 END AS promoter_end, -- rightmost position, regardless of strand
    CASE WHEN strand = 1  THEN start_position ELSE end_position END AS gene_5p_end -- Position on the genome of gene's 5'end.
from gene_attributes;

-- Tables w/ Transcription Start Sites (position and count):
/*
Mouse: mouse_macrophage_tss_bed
Pig:   cage050810_tss 
Human: macrophage_monocyte_derived_hg19_ctss_bed; SELECT * from macrophage_monocyte_derived_hg19_ctss_bed limit 1000;
*/



-- Write out FASTA sequences for later use:
select write_fasta($$
  select reference_gene_name || '|' || gene_exon_intron.species || '|' ||  gene_exon_intron.ensembl_gene_id AS fasta_name, 
--  gene_exon_intron,
  substring(gene_exon_intron from 0 for 500 ),
  gene_exon_intron.species 
  from gene_exon_intron inner join reference_genes on reference_genes.ensembl_gene_id = gene_exon_intron.ensembl_gene_id 
--   where reference_gene_name like 'TNF'
  order by fasta_name, species
  $$, 80);



---------------------------------[ Tritume ]-----------------------------------

select * from transcript_attributes; 3 868 535
select * from tss_bed where chr like '17' and tss_pos between 9681681-10000 and 9681681+10000;

select * from cage050810 where rname like '17' and pos between 9681681-10000 and 9681681+10000;

select * from pileup_rnaseq where rname like '17' and pos between 9681681-100 and 9681681+100;


/*
select read_table($$ file:'D:/Tritume/reference_genes_ct.txt', table: 'immunogenes.reference_genes_ct', header:True, overwrite:True $$);
comment on table reference_genes_ct is 'Reference name of genes of immunological interest.';
select * from reference_genes_ct order by reference_gene_name;

select read_table($$ file:'D:/Tritume/reference_genes_norm.txt', table: 'immunogenes.reference_genes', header:True, overwrite:True $$);
comment on table reference_genes is 'Reference name of genes of immunological interest.';
select * from immunogenes.reference_genes order by reference_gene_name;
*/



select distinct rna from "FANTOM4".mouse_macrophage_tss_bed ;

-- drop table sequences_coding_flank_genes;
create table sequences_coding_flank_genes AS(
  select coding_gene_flank_up300_hs.ensembl_gene_id,  
       300::int AS upstream_bp,
       100::int AS downstream_bp,
       'hsapiens' AS species,
       coding_gene_flank_up300_hs.coding_gene_flank || coding_gene_flank_down100_hs.coding_gene_flank AS sequence
  from coding_gene_flank_up300_hs inner join coding_gene_flank_down100_hs on coding_gene_flank_up300_hs.ensembl_gene_id = coding_gene_flank_down100_hs.ensembl_gene_id

  union select coding_gene_flank_up300_mm.ensembl_gene_id,  
       300::int AS upstream_bp,
       100::int AS downstream_bp,
       'mmusculus' AS species,
       coding_gene_flank_up300_mm.coding_gene_flank || coding_gene_flank_down100_mm.coding_gene_flank AS sequence
  from coding_gene_flank_up300_mm inner join coding_gene_flank_down100_mm on coding_gene_flank_up300_mm.ensembl_gene_id = coding_gene_flank_down100_mm.ensembl_gene_id

  union select coding_gene_flank_up300_ss.ensembl_gene_id,  
       300::int AS upstream_bp,
       100::int AS downstream_bp,
       'sscrofa' AS species,
       coding_gene_flank_up300_ss.coding_gene_flank || coding_gene_flank_down100_ss.coding_gene_flank AS sequence
  from coding_gene_flank_up300_ss inner join coding_gene_flank_down100_ss on coding_gene_flank_up300_ss.ensembl_gene_id = coding_gene_flank_down100_ss.ensembl_gene_id
  );
select count(*), species from sequences_coding_flank_genes group by species;
select '>'||ensembl_gene_id, "sequence", species from sequences_coding_flank_genes order by species;


SELECT read_table($$ file: 'D:/Tritume/coding_gene_flank_down100_hs.txt', temp:True, header:True, overwrite:True $$);
SELECT read_table($$ file: 'D:/Tritume/coding_gene_flank_down100_mm.txt', temp:True, header:True, overwrite:True $$);
SELECT read_table($$ file: 'D:/Tritume/coding_gene_flank_down100_ss.txt', temp:True, header:True, overwrite:True $$);

coding_gene_flank_up300_hs.txt
coding_gene_flank_down100_hs.txt
coding_gene_flank_up300_mm.txt
coding_gene_flank_down100_mm.txt
coding_gene_flank_up300_ss.txt
coding_gene_flank_down100_ss.txt


-- species, reference_gene_name, ensembl_gene_id, ensembl_transcript_id, external_transcript_id, chr, transcript_start, transcript_end, 
-- strand, tss_id, tss_pos, transcript_5p_end, tss_transcript_dist, tpm, rna

/* Identify for each gene where the highest peak is located. Then get the all the TSSs 
   in the region -1000/+1000 bp around each peak
*/

create temp table peaks AS (
    -- Max tpm within of each gene
    select ensembl_gene_id, max(tpm) AS tpm_peak from tss_immunogenes group by ensembl_gene_id
    ); 
-- drop table peak_pos;
create temp table peak_pos AS (
    -- Add position info to each tpm peak in "peaks"
    select distinct tss_id AS tss_peak_id, tss_immunogenes.ensembl_gene_id, tss_immunogenes.tpm, tss_immunogenes.tss_pos, species, reference_gene_name, strand, chr
    from tss_immunogenes inner join peaks on tss_immunogenes.ensembl_gene_id = peaks.ensembl_gene_id and tss_immunogenes.tpm= peaks.tpm_peak
  ); -- select * from peak_pos order by reference_gene_name, species;
-- Find the nearest transcript to each peak.

-- drop table immunogenes_peaks_annotated;
create table immunogenes_peaks_annotated AS(
select peak_pos.*, transcript_attributes.ensembl_transcript_id, transcript_5p_end,
  CASE WHEN transcript_attributes.strand = '+' THEN tss_pos - transcript_5p_end ELSE transcript_5p_end - tss_pos END AS transcript_peak_dist, -- Distance between the transcript 5'end and the TSS peak. (-)ve means TSS is upstream,
  dense_rank() OVER (PARTITION BY transcript_attributes.ensembl_gene_id ORDER BY abs(tss_pos - transcript_5p_end)) AS transcript_peak_dist_rank
from peak_pos inner join transcript_attributes on transcript_attributes.ensembl_gene_id = peak_pos.ensembl_gene_id
order by species, reference_gene_name, transcript_peak_dist
);
comment on table immunogenes_peaks_annotated is 'For each peak in each gene report the distance from each transcript. "Peak" is the position with the highest TPM within each gene.';
comment on column immunogenes_peaks_annotated.transcript_peak_dist IS 'Distance between the transcript 5p-end and the TSS peak. (-)ve means TSS is upstream';
comment on column immunogenes_peaks_annotated.transcript_peak_dist_rank IS 'Rank of the distance between the transcript 5p-end and the TSS peak. The nearest transcript has rank 1';

-- Annotate the region around each peak.
-- drop table tss_immunogenes_peaks;
create table tss_immunogenes_peaks AS(
    select distinct tss_immunogenes.species, 
                    tss_immunogenes.reference_gene_name, 
                    tss_immunogenes.ensembl_transcript_id, 
                    tss_immunogenes.transcript_5p_end, 
                    tss_immunogenes.ensembl_gene_id, tss_immunogenes.chr, 
                    tss_immunogenes.tss_pos, 
                    tss_immunogenes.tpm,
                    immunogenes_peaks_annotated.transcript_peak_dist_rank
    from tss_immunogenes inner join immunogenes_peaks_annotated on tss_immunogenes.ensembl_transcript_id = immunogenes_peaks_annotated.ensembl_transcript_id
    -- Change here the span around each peak:
    where tss_immunogenes.tss_pos between immunogenes_peaks_annotated.tss_pos - 100 and immunogenes_peaks_annotated.tss_pos + 100
    order by species, reference_gene_name, tss_immunogenes.tss_pos, transcript_5p_end
    );
comment on table tss_immunogenes_peaks is 'TSS peak -/+ bp in each gene See 20110119_immuno_genes.sql';
select distinct * from tss_immunogenes_peaks where reference_gene_name like 'ACP5' and species like '%' order by species, tss_pos, transcript_5p_end;

-- For plotting pourposes: Add blank rows around each peak if no tags have been detected. This is to be able to use plot(type= 'l'), rather than type='h'.

-- drop table genome_pos;
create temp table genome_pos AS(
  -- This table generates +/-100 bp around each peak. Genes not detected are included
  select distinct 
         transcript_attributes.ensembl_gene_id, 
         transcript_attributes.reference_gene_name, 
         transcript_attributes.species, 
         transcript_attributes.strand, 
         transcript_attributes.chromosome_name AS chr, 
         CASE WHEN tss_pos is null THEN 0 ELSE tss_pos END AS peak_pos,
         (CASE WHEN peak_pos.tss_pos is null THEN 0 - generate_series(-100, 100) ELSE peak_pos.tss_pos - generate_series(-100, 100) END )::int AS pos 
  from peak_pos right join transcript_attributes on transcript_attributes.ensembl_gene_id = peak_pos.ensembl_gene_id
  order by transcript_attributes.reference_gene_name, transcript_attributes.species
);

-- drop table tss_immunogenes_peaks_plot;
create table tss_immunogenes_peaks_plot AS (
select distinct 
   genome_pos.ensembl_gene_id, 
   genome_pos.pos,
   genome_pos.peak_pos,
   CASE WHEN genome_pos.strand like '+' THEN (genome_pos.peak_pos - genome_pos.pos) ELSE (genome_pos.pos - genome_pos.peak_pos) END AS peak_dist, -- Distance of this position from the peak. (-)ve means position is upstream
   genome_pos.reference_gene_name, genome_pos.species, genome_pos.chr,
   tss_immunogenes_peaks.ensembl_transcript_id, tss_immunogenes_peaks.transcript_start, tss_immunogenes_peaks.transcript_end, 
   tss_immunogenes_peaks.strand, tss_immunogenes_peaks.tss_id, tss_immunogenes_peaks.tss_pos, 
   tss_immunogenes_peaks.transcript_5p_end, tss_immunogenes_peaks.tss_transcript_dist, 
   CASE WHEN tss_immunogenes_peaks.tpm is null THEN 0 ELSE tss_immunogenes_peaks.tpm END, 
   tss_immunogenes_peaks.rna   
  from tss_immunogenes_peaks right join genome_pos on
     -- Set here the span around the peaks.
--    (select ensembl_gene_id, reference_gene_name, species, strand, chr, tss_pos AS peak_pos, (peak_pos.tss_pos - generate_series(-100, 100))::int AS pos from peak_pos) as genome_pos on
  tss_immunogenes_peaks.ensembl_gene_id = genome_pos.ensembl_gene_id and 
  tss_immunogenes_peaks.tss_pos = genome_pos.pos
  order by genome_pos.reference_gene_name, genome_pos.species, genome_pos.pos
);
select distinct * from tss_immunogenes_peaks_plot where reference_gene_name like 'ACTB' and species like 'mmusculus' order by pos, transcript_5p_end;
select * from transcript_attributes;

-- Produce some summary stats:
copy (
    select 
      transcript_attributes.species,
      min(transcript_attributes.reference_gene_name) AS reference_gene_name,
      transcript_attributes.ensembl_transcript_id,
      max(tpm) AS max_tpm,
      sum(tpm) AS sum_tpm,
      avg(tss_transcript_dist*tpm) AS tss_to_transcript_dist_weighted,
      count(distinct tss_id) AS no_distinct_tss
    from tss_immunogenes right join transcript_attributes on transcript_attributes.species = tss_immunogenes.species AND
                                   transcript_attributes.reference_gene_name = tss_immunogenes.reference_gene_name
    group by transcript_attributes.species, transcript_attributes.ensembl_transcript_id
    order by reference_gene_name, transcript_attributes.species
    ) TO 'D:/Tritume/tss_immunogenes_summary.txt' with csv header delimiter E'\t';

    select 
      min(tss_immunogenes.species) AS species,
      min(tss_immunogenes.reference_gene_name) AS reference_gene_name,
      tss_immunogenes.ensembl_transcript_id,
      max(tpm) AS max_tpm,
      sum(tpm) AS sum_tpm,
      avg(tss_transcript_dist) AS tss_to_transcript_dist_weighted,
      count(distinct tss_id) AS no_distinct_tss
    from tss_immunogenes group by tss_immunogenes.ensembl_transcript_id
    order by reference_gene_name, species



