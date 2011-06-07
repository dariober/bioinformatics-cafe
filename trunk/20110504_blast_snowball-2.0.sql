/*
   Annotation of the Snowball array
   See Labbook 04/05/2011
   Data in F:/data/20110504_snowball_annotation/
*/


-- Import gene_info from NCBI ftp://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz
select read_table($$ file:'F:/data/20110504_snowball_annotation/gene_info', table: 'public.ncbi_gene_info', skip: 1, overwrite:True, limit: -1,
    header:['tax_id', 'gene_id', 'symbol', 'locus_tag', 'synonyms', 'dbxrefs', 'chromosome', 'map_location', 'description', 'type_of_gene', 'symbol_from_nomenclature_authority', 'full_name_from_nomenclature_authority', 'nomenclature_status', 'other_designations', 'modification_date'],
    data_type: ['int', 'int', ['text'] * 13]
$$);
comment on table ncbi_gene_info is 'Table imported form NCBI';
create unique index ind_ncbi_gene_info on ncbi_gene_info(gene_id);
vacuum analyze ncbi_gene_info;
select * from ncbi_gene_info limit 10;


-------------------------------------------------------------------------------
--  IMPORT BLASTN RESULTS
-------------------------------------------------------------------------------

-- Import FASTA sequence names used for BLASTs.
-- BLAST output has the these names truncated to the first whte-space.
select read_fasta('D:/Tritume/non_psr_fastas_cat.fa', 'tmp_fa');
create table snowball.non_psr_fastas AS(
    select 
        substring(sequence_name from 1 for (position(' ' in sequence_name)-1 )) AS qseqid, * 
    from 
        tmp_fa
);
comment on table non_psr_fastas is 'Fasta files in non_PSR_fastas.zip concatenated and sequence names parsed to match the BLAST names (i.e. names trimmed to first white space found). See 20110504_blast_snowball-2.0.sql'
drop table tmp_fa;
select * from non_psr_fastas limit 10;


-- BLASTN oputput
drop table snowball.blastn_nonpsr_refseqrna;
select read_table($$ file:'D:/Tritume/non_psr_fastas_cat.blastnout', table: 'blastn_nonpsr_refseqrna_tmp',
    header:['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'], overwrite:True, temp:True $$)

-- Extract gi and Refseq ID
create table snowball.blastn_nonpsr_refseqrna AS(
    select get_fasta_attribute(sseqid, 'gi') AS gi, get_fasta_attribute(sseqid, 'ref') AS refseq_id, * from blastn_nonpsr_refseqrna_tmp
);
comment on table blastn_nonpsr_refseqrna is $$
    'Output from BLASTN: non_psr_fastas_cat.fa Database refseq_rna
    BLAST command was (input file splitted in 26 chunks): 
    Subject database: ftp://ftp.ncbi.nih.gov/blast/db/refseq_rna.tar.gz
    blastn -query non_psr_fastas_cat.fa.001 -db $blast_db -task blastn -evalue 0.001 -best_hit_overhang 0.2 -best_hit_score_edge 0.05 -word_size 11 -out non_psr_fastas_cat.fa.001.blastnout -outfmt 6
    Output from blastn. See labbook 05/05/2011 script 200110504_blast_snowball.sql'
  $$;
select * from blastn_nonpsr_refseqrna limit 10;
select count(*) AS no_hits, count(distinct qseqid) AS no_queries, count(distinct sseqid) AS no_subjects from blastn_nonpsr_refseqrna;


-- Output from MEGABLAST: non_psr_fastas_cat.fa Database refseq_rna
-- Subject database: ftp://ftp.ncbi.nih.gov/blast/db/refseq_rna.tar.gz
-- blastn -query non_psr_fastas_cat.fa -db $blast_db -task megablast -evalue 0.01 -out non_psr_fastas_cat_refseq_rna.blastout -outfmt 6
select read_table($$ file:'D:/Tritume/non_psr_fastas_cat_refseq_rna.blastout.gz', table:'snowball.megablast_nonpsr_refseqrna',
    header:['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'], overwrite:True  $$)
comment on table megablast_nonpsr_refseqrna is 
    'Output from megablast: FASTA in non_PSR_fastas.zip (concatenated) against blast db ftp://ftp.ncbi.nih.gov/blast/db/refseq_rna.tar.gz. See labbook 04/05/2011;
blastn -query non_psr_fastas_cat.fa -db $blast_db -task megablast -evalue 0.01 -out non_psr_fastas_cat.blastout -outfmt 6'
select * from megablast_nonpsr_refseqrna limit 10;
select count(*) AS no_hits, count(distinct qseqid) AS no_queries, count(distinct sseqid) AS no_subjects from megablast_nonpsr_refseqrna;

-------------------------------------------------------------------------------
-- ANNOTATION refseq_id
-------------------------------------------------------------------------------

-- Accession IDs to taxonomic id
-- Import this file: 
select read_table($$ file:'D:/Tritume/release46.accession2geneid.gz', overwrite:True, table:'public.ncbi_release46_accession2geneid', header: ['taxonomic_id', 'entrez_gene_id', 'transcript_accession_ver', 'protein_accession_ver' ]$$)
comment on table ncbi_release46_accession2geneid is $$
    'Mapping of accession ID to taxonomic ID from ftp://ftp.ncbi.nih.gov/refseq/release/release-catalog/release46.accession2geneid.gz'
$$;
select * from ncbi_release46_accession2geneid limit 10;

-- Catalog of refseq
select read_table($$ file:'D:/Tritume/RefSeq-release46.catalog', table:'public.ncbi_refseq_release46_catalog', overwrite:True $$)
select * from ncbi_refseq_release46_catalog limit 10;
create table ncbi_taxon_id2name AS(
    select distinct v1 as taxonomic_id, v2 as species_name from ncbi_refseq_release46_catalog
    );
alter table ncbi_taxon_id2name set schema snowball;
comment on table ncbi_taxon_id2name is 'See 20110504... Taxonimic ID from NCBI linked to species name as extracted from ftp://ftp.ncbi.nih.gov/refseq/release/release-catalog/RefSeq-release46.catalog.gz';
select * from ncbi_taxon_id2name limit 10;

-- Get the species names corresponding to each taxonomic id:
create temp table species2refseq AS(
    select ncbi_release46_accession2geneid.*, species_name 
    from ncbi_release46_accession2geneid inner join ncbi_taxon_id2name on ncbi_release46_accession2geneid.taxonomic_id = ncbi_taxon_id2name.taxonomic_id
);
select * from species2refseq limit 10;
select count(*) from species2refseq limit 10;

-- Get number of hits by species:
select count(*)
from blastn_nonpsr_refseqrna inner join species2refseq on
    transcript_accession_ver = refseq_id;

drop table blastn_nonpsr_refseqrna_species;
create table snowball.blastn_nonpsr_refseqrna_species AS (
  select 
      ncbi_release46_accession2geneid.entrez_gene_id,
      blastn_nonpsr_refseqrna.*, 
      species_name, 
      ncbi_release46_accession2geneid.taxonomic_id
  from 
      blastn_nonpsr_refseqrna 
  left join 
      species2refseq 
          on species2refseq.transcript_accession_ver = refseq_id -- to add species name
  left join 
      ncbi_release46_accession2geneid
          on ncbi_release46_accession2geneid.transcript_accession_ver = blastn_nonpsr_refseqrna.refseq_id -- To add gene_id
);

-- Add gene information to blast table with all hits (blastn_nonpsr_refseqrna_species).
-- Reminder: blastn_nonpsr_refseqrna_species is the same as blastn_nonpsr_refseqrna but with some columns 
-- added. Therefore blastn_nonpsr_refseqrna is redundant.
alter table blastn_nonpsr_refseqrna_species add column symbol text;
update blastn_nonpsr_refseqrna_species set symbol = ncbi_gene_info.symbol from ncbi_gene_info
where blastn_nonpsr_refseqrna_species.entrez_gene_id = ncbi_gene_info.gene_id;
comment on column blastn_nonpsr_refseqrna_species.symbol  is 'Gene symbol from NCBI gene_info correposponding to entrez_gene_id = gene_id';
-- Gene Symbol Description
alter table blastn_nonpsr_refseqrna_species add column description text;
update blastn_nonpsr_refseqrna_species set description = ncbi_gene_info.description from ncbi_gene_info
where blastn_nonpsr_refseqrna_species.entrez_gene_id = ncbi_gene_info.gene_id;
comment on column blastn_nonpsr_refseqrna_species.description  is 'Symbol description from NCBI gene_info correposponding to entrez_gene_id = gene_id';
-- Type of gene
alter table blastn_nonpsr_refseqrna_species add column type_of_gene text;
update blastn_nonpsr_refseqrna_species set type_of_gene = ncbi_gene_info.type_of_gene from ncbi_gene_info
where blastn_nonpsr_refseqrna_species.entrez_gene_id = ncbi_gene_info.gene_id;
comment on column blastn_nonpsr_refseqrna_species.type_of_gene  is 'type_of_gene from NCBI gene_info correposponding to entrez_gene_id = gene_id';

-- Type of gene
alter table blastn_nonpsr_refseqrna_species add column type_of_gene text;
update blastn_nonpsr_refseqrna_species set type_of_gene = ncbi_gene_info.type_of_gene from ncbi_gene_info
where blastn_nonpsr_refseqrna_species.entrez_gene_id = ncbi_gene_info.gene_id;
comment on column blastn_nonpsr_refseqrna_species.type_of_gene  is 'type_of_gene from NCBI gene_info correposponding to entrez_gene_id = gene_id';

-- Add probe_set_id:
alter table blastn_nonpsr_refseqrna_species add column probe_set_id text;
update blastn_nonpsr_refseqrna_species set probe_set_id = non_psr_fastas.probe_set_id from non_psr_fastas
where blastn_nonpsr_refseqrna_species.qseqid = non_psr_fastas.qseqid;

-- Hit rank within species, for each query
alter table blastn_nonpsr_refseqrna_species rename to blastn_nonpsr_refseqrna_species2;
alter table blastn_nonpsr_refseqrna_species2 drop column bitscore_rank_species;
create table snowball.blastn_nonpsr_refseqrna_species AS(
    select *, row_number() over (partition by qseqid, species_name order by bitscore desc, qseqid asc) as bitscore_rank_species from blastn_nonpsr_refseqrna_species2
);
drop table blastn_nonpsr_refseqrna_species2;

-- Hit rank within queries:
alter table blastn_nonpsr_refseqrna_species rename to blastn_nonpsr_refseqrna_species2;
alter table blastn_nonpsr_refseqrna_species2 drop column bitscore_rank;
create table snowball.blastn_nonpsr_refseqrna_species AS(
    select *, row_number() over (partition by qseqid order by bitscore desc, species_name asc) as bitscore_rank from blastn_nonpsr_refseqrna_species2
);
drop table blastn_nonpsr_refseqrna_species2;
comment on table blastn_nonpsr_refseqrna_species is 'Output of blastn blastn_nonpsr_refseqrna with added: entrez_gene_id (from ncbi_release46_accession2geneid) and species (from ncbi_taxon_id2name). This table replaces blastn_nonpsr_refseqrna. See 20110504_blast_snowball-2.0.sql';

select count(*) from blastn_nonpsr_refseqrna_species;
drop table blastn_nonpsr_refseqrna; -- This table replaced by blastn_nonpsr_refseqrna_species.
create index ind_blastn_nonpsr_refseqrna_species on blastn_nonpsr_refseqrna_species(entrez_gene_id); -- drop index ind_blastn_nonpsr_refseqrna_species;
vacuum analyze blastn_nonpsr_refseqrna_species;
select * from blastn_nonpsr_refseqrna_species where species_name like 'Sus%' limit 1000;

-- No. hits by species:
select species_name,
       count(qseqid) AS no_query_hits, 
       count(distinct qseqid) no_distinct_query_hits, 
       count(distinct refseq_id) as no_distinct_refseq 
from blastn_nonpsr_refseqrna_species
group by species_name
order by no_distinct_query_hits desc;

-- No. hits by best hit
select count(distinct qseqid) from blastn_nonpsr_refseqrna_species where evalue <= 1e-9;
select count()
select count(*), substring(refseq_id from 1 for 2) from blastn_nonpsr_refseqrna_species group by substring(refseq_id from 1 for 2);
limit 10;


-------------------------------------------------------------------------------
-- Get best hits, one per query sequence (for some definition of best)
-------------------------------------------------------------------------------

/*
   A query (~probe) hits multiple genes within the same species, for simplicity 
   we want only one hit per query sequence.
   If a query has multiple hits within the same species, get the best hit based on
   min evalue (whatever that is, see temp table best_hits). 
   To chose which species use the pipeline below where the precedence is given below
   (in general, start from non LOC genes).
*/

alter table blastn_nonpsr_refseqrna_species drop column approved_hit;
alter table blastn_nonpsr_refseqrna_species add column approved_hit boolean;
update blastn_nonpsr_refseqrna_species set approved_hit = False where approved_hit is null;
create unique index ind_blast on blastn_nonpsr_refseqrna_species (qseqid, approved_hit) where approved_hit is True;


/* ----------------------------------------------------------------------------
 This section generated with 20110523_snowball_update.py.
 Don't edit the SQL code here. Edit the python script instead!!!!
---------------------------------------------------------------------------- */
-- Not LOC annotation
update blastn_nonpsr_refseqrna_species set approved_hit= True
where species_name = 'Homo sapiens' and symbol not like 'LOC%' and evalue < 1e-9 and bitscore_rank_species = 1 and
qseqid not in (select distinct qseqid from blastn_nonpsr_refseqrna_species where approved_hit is true);

update blastn_nonpsr_refseqrna_species set approved_hit= True
where species_name = 'Sus scrofa' and symbol not like 'LOC%' and evalue < 1e-9 and bitscore_rank_species = 1 and
qseqid not in (select distinct qseqid from blastn_nonpsr_refseqrna_species where approved_hit is true);

update blastn_nonpsr_refseqrna_species set approved_hit= True
where species_name = 'Bos taurus' and symbol not like 'LOC%' and evalue < 1e-9 and bitscore_rank_species = 1 and
qseqid not in (select distinct qseqid from blastn_nonpsr_refseqrna_species where approved_hit is true);

update blastn_nonpsr_refseqrna_species set approved_hit= True
where species_name = 'Pan troglodytes' and symbol not like 'LOC%' and evalue < 1e-9 and bitscore_rank_species = 1 and
qseqid not in (select distinct qseqid from blastn_nonpsr_refseqrna_species where approved_hit is true);

update blastn_nonpsr_refseqrna_species set approved_hit= True
where species_name = 'Mus musculus' and symbol not like 'LOC%' and evalue < 1e-9 and bitscore_rank_species = 1 and
qseqid not in (select distinct qseqid from blastn_nonpsr_refseqrna_species where approved_hit is true);

update blastn_nonpsr_refseqrna_species set approved_hit= True
where species_name = 'Canis lupus familiaris' and symbol not like 'LOC%' and evalue < 1e-9 and bitscore_rank_species = 1 and
qseqid not in (select distinct qseqid from blastn_nonpsr_refseqrna_species where approved_hit is true);

update blastn_nonpsr_refseqrna_species set approved_hit= True
where species_name = 'Pongo abelii' and symbol not like 'LOC%' and evalue < 1e-9 and bitscore_rank_species = 1 and
qseqid not in (select distinct qseqid from blastn_nonpsr_refseqrna_species where approved_hit is true);

update blastn_nonpsr_refseqrna_species set approved_hit= True
where species_name = 'Equus caballus' and symbol not like 'LOC%' and evalue < 1e-9 and bitscore_rank_species = 1 and
qseqid not in (select distinct qseqid from blastn_nonpsr_refseqrna_species where approved_hit is true);

update blastn_nonpsr_refseqrna_species set approved_hit= True
where species_name = 'Rattus norvegicus' and symbol not like 'LOC%' and evalue < 1e-9 and bitscore_rank_species = 1 and
qseqid not in (select distinct qseqid from blastn_nonpsr_refseqrna_species where approved_hit is true);

update blastn_nonpsr_refseqrna_species set approved_hit= True
where species_name = 'Macaca mulatta' and symbol not like 'LOC%' and evalue < 1e-9 and bitscore_rank_species = 1 and
qseqid not in (select distinct qseqid from blastn_nonpsr_refseqrna_species where approved_hit is true);

-- All remaining not LOC genes:
update blastn_nonpsr_refseqrna_species set approved_hit= True
where species_name like '%' and symbol not like 'LOC%' and bitscore_rank = 1 and
qseqid not in (select distinct qseqid from blastn_nonpsr_refseqrna_species where approved_hit is true);


-- Use LOC annotation, whatever the evalue:
update blastn_nonpsr_refseqrna_species set approved_hit= True
where species_name = 'Homo sapiens' and bitscore_rank_species = 1 and
qseqid not in (select distinct qseqid from blastn_nonpsr_refseqrna_species where approved_hit is true);

update blastn_nonpsr_refseqrna_species set approved_hit= True
where species_name = 'Sus scrofa' and bitscore_rank_species = 1 and
qseqid not in (select distinct qseqid from blastn_nonpsr_refseqrna_species where approved_hit is true);

update blastn_nonpsr_refseqrna_species set approved_hit= True
where species_name = 'Bos taurus' and bitscore_rank_species = 1 and
qseqid not in (select distinct qseqid from blastn_nonpsr_refseqrna_species where approved_hit is true);

update blastn_nonpsr_refseqrna_species set approved_hit= True
where species_name = 'Pan troglodytes' and bitscore_rank_species = 1 and
qseqid not in (select distinct qseqid from blastn_nonpsr_refseqrna_species where approved_hit is true);

update blastn_nonpsr_refseqrna_species set approved_hit= True
where species_name = 'Mus musculus' and bitscore_rank_species = 1 and
qseqid not in (select distinct qseqid from blastn_nonpsr_refseqrna_species where approved_hit is true);

update blastn_nonpsr_refseqrna_species set approved_hit= True
where species_name = 'Canis lupus familiaris' and bitscore_rank_species = 1 and
qseqid not in (select distinct qseqid from blastn_nonpsr_refseqrna_species where approved_hit is true);

update blastn_nonpsr_refseqrna_species set approved_hit= True
where species_name = 'Pongo abelii' and bitscore_rank_species = 1 and
qseqid not in (select distinct qseqid from blastn_nonpsr_refseqrna_species where approved_hit is true);

update blastn_nonpsr_refseqrna_species set approved_hit= True
where species_name = 'Equus caballus' and bitscore_rank_species = 1 and
qseqid not in (select distinct qseqid from blastn_nonpsr_refseqrna_species where approved_hit is true);

update blastn_nonpsr_refseqrna_species set approved_hit= True
where species_name = 'Rattus norvegicus' and bitscore_rank_species = 1 and
qseqid not in (select distinct qseqid from blastn_nonpsr_refseqrna_species where approved_hit is true);

update blastn_nonpsr_refseqrna_species set approved_hit= True
where species_name = 'Macaca mulatta' and bitscore_rank_species = 1 and
qseqid not in (select distinct qseqid from blastn_nonpsr_refseqrna_species where approved_hit is true);

-- Use everything else:
update blastn_nonpsr_refseqrna_species set approved_hit= True
where bitscore_rank = 1 and
qseqid not in (select distinct qseqid from blastn_nonpsr_refseqrna_species where approved_hit is true);

/* ----------------------------------------------------------------------------
 End of section generated by 20110523_snowball_update.py 
---------------------------------------------------------------------------- */


select * from blastn_nonpsr_refseqrna_best limit 10000;
select count(qseqid) from blastn_nonpsr_refseqrna_best where symbol not like 'LOC%'; -- 29744
select species_name, count(*) from blastn_nonpsr_refseqrna_best group by species_name order by count desc;
select species_name, count(*) from blastn_nonpsr_refseqrna_best where symbol like 'LOC%' group by species_name order by count desc;
select count(*) from blastn_nonpsr_refseqrna_best where refseq_id is not null; -- 42255

/* ----------------------------------------------------------------------------
   AFFYMETRIX ARRAY
   Import affymetrix probeset ids and annotation from FIOS 
---------------------------------------------------------------------------- */

-- This Excel file from Tom F. It should be the annotation from FIOS
select read_table($$ file:'F:/data/20110504_snowball_annotation/Snowball_array_annotation_8-3-11.xls', table:'snowball.array_annotation_fios', header:True $$)
comment on table array_annotation_fios is 'Annotation from FIOS/Tom Freeman. Source file is F:/data/20110504_snowball_annotation/Snowball_array_annotation_8-3-11.xls. See also 20110504_blast_snowball-2.0.sql';
comment on column array_annotation_fios."Seqname" is 'Column Seqname and SNOWBALL_ID are identical. Test with SELECT * FROM array_annotation_fios';
comment on column array_annotation_fios."SNOWBALL_ID" is 'Column Seqname and SNOWBALL_ID are identical. Test with SELECT * FROM array_annotation_fios';
select * from array_annotation_fios limit 1000;

-- Try to reverse Excel dates for some gene symbols incorrectly converted.
-- These symbols found in Excel with "Find -> format -> custom -> 'mmm-yy'"
UPDATE array_annotation_fios SET "Symbol_1" =
    CASE WHEN "Symbol_1" like 'Sep-0%' THEN replace("Symbol_1", 'Sep-0', 'SEPT')
         WHEN "Symbol_1" like 'Mar-0%' THEN replace("Symbol_1", 'Mar-0', 'MARCH')
         WHEN "Symbol_1" like 'Sep-%' THEN replace("Symbol_1", 'Sep-', 'SEPT')
         WHEN "Symbol_1" like 'Mar-%' THEN replace("Symbol_1", 'Mar-', 'MARCH') END
where "Probeset" in ('SNOWBALL_011586_st', 'SNOWBALL_027725_st', 'SNOWBALL_031750_st', 'SNOWBALL_004470_st', 'SNOWBALL_018442_st', 'SNOWBALL_032796_st', 'SNOWBALL_006726_st', 'SNOWBALL_018223_st', 'SNOWBALL_000190_st', 'SNOWBALL_000191_st', 'SNOWBALL_017385_s_st', 'SNOWBALL_024598_st', 'SNOWBALL_028035_st', 'SNOWBALL_035960_st', 'SNOWBALL_041960_st', 'SNOWBALL_018725_st', 'SNOWBALL_001879_s_st', 'SNOWBALL_001880_s_st', 'SNOWBALL_001881_s_st', 'SNOWBALL_001882_st', 'SNOWBALL_001883_s_st', 'SNOWBALL_001884_s_st', 'SNOWBALL_001886_s_st', 'SNOWBALL_001887_st', 'SNOWBALL_001888_s_st', 'SNOWBALL_001889_s_st', 'SNOWBALL_001890_s_st', 'SNOWBALL_001891_st', 'SNOWBALL_035657_st', 'SNOWBALL_023846_st', 'SNOWBALL_019887_st', 'SNOWBALL_017393_st', 'SNOWBALL_024646_st', 'SNOWBALL_027764_st', 'SNOWBALL_045020_st', 'SNOWBALL_044274_st', 'SNOWBALL_007450_st', 'SNOWBALL_014173_st', 'SNOWBALL_047476_st', 'SNOWBALL_047477_s_st', 'SNOWBALL_047502_st', 'SNOWBALL_014360_st', 'SNOWBALL_034615_st', 'SNOWBALL_018541_st', 'SNOWBALL_011450_st', 'SNOWBALL_004889_st', 'SNOWBALL_004425_st', 'SNOWBALL_021021_s_st', 'SNOWBALL_022338_st')

-- Match FASTA sequence IDs to probe_set_ids:
select "Probeset", substring("Probeset" from length("Probeset")-2 for 3), substring("Probeset" from 1 for length("Probeset")-3) from array_annotation_fios limit 1000;



create temp table probeset_1 AS(
  select non_psr_fastas.qseqid, array_annotation_fios."Probeset" from non_psr_fastas inner join array_annotation_fios on qseqid = substring("Probeset" from 1 for length("Probeset")-3)
); --40165
create temp table probeset_2 AS(
  select non_psr_fastas.qseqid, array_annotation_fios."Probeset" from non_psr_fastas inner join array_annotation_fios on qseqid = substring("Probeset" from 1 for length("Probeset")-5)
  where non_psr_fastas.qseqid not in (select qseqid from probeset_1) and 
        array_annotation_fios."Probeset" not in (select "Probeset" from probeset_1)
); -- 7557
create temp table probeset as(
    select * from probeset_1 
    union 
    select *from probeset_2 
    );
select count(*) from array_annotation_fios where "Probeset" not in (select "Probeset" from probeset); -- 123
-- All probes have been linked to qseqid (i.e. to FASTA sequences):
-- 40165 + 7557 + 123 = 47845 = select count(*) array_annotation_fios 

-- Add probeset id to fasta file;
alter table non_psr_fastas add column probe_set_id text;
update non_psr_fastas set probe_set_id = "Probeset" from probeset where probeset.qseqid = non_psr_fastas.qseqid;
select * from non_psr_fastas where probe_set_id is not null limit 10;

-- Add probeset id to "best" annotation file;
alter table blastn_nonpsr_refseqrna_best add column probe_set_id text;
update blastn_nonpsr_refseqrna_best set probe_set_id = "Probeset" from probeset where probeset.qseqid = blastn_nonpsr_refseqrna_best.qseqid;
select * from blastn_nonpsr_refseqrna_best where probe_set_id is not null limit 10;

-------------------------------------------------------------------------------
-- Consensus annotation:
-------------------------------------------------------------------------------
create temp table blastn_nonpsr_refseqrna_best AS(select * from blastn_nonpsr_refseqrna_species where approved_hit is True);
drop table blastn_nonpsr_refseqrna_consensus;
create table snowball.blastn_nonpsr_refseqrna_consensus AS(
    select array_annotation_fios."Probeset" as probe_set_id,
       blastn_nonpsr_refseqrna_best.qseqid,    -- Identifier of the FASTA sequence from which probe set was designed
       blastn_nonpsr_refseqrna_best.refseq_id, -- Best refseq hit for this sequence
       blastn_nonpsr_refseqrna_best.entrez_gene_id, -- entrez id associated to this refseq_id (from NCBI release46.accession2geneid)
       blastn_nonpsr_refseqrna_best.species_name, -- Species from which this refseq comes from
       blastn_nonpsr_refseqrna_best.symbol AS symbol_ncbi, -- Gene symbol associated to this refseq_id (from NCBI gene_info table)
       array_annotation_fios."Symbol_1" as symbol_fios, -- Gene symbol associated to this sequence given by FIOS.
       CASE WHEN symbol like 'LOC%' and "Symbol_1" is not null THEN "Symbol_1"
            WHEN symbol is null and "Symbol_1" is not null THEN "Symbol_1"
            ELSE symbol END as symbol_consensus, -- Gene symbol for this probeset taking either symbol_ncbi, if available, or symbol_fios.
       blastn_nonpsr_refseqrna_best.description -- Gene description associated to this entrez_gene_id (from NCBI gene_info table)
    from array_annotation_fios left join blastn_nonpsr_refseqrna_best on array_annotation_fios."Probeset" = blastn_nonpsr_refseqrna_best.probe_set_id
);
comment on table blastn_nonpsr_refseqrna_consensus is 'Annotation of the snowball array using either the results from BLASTN vs refseq (blastn_nonpsr_refseqrna_best) or fios annotation. See 20110504_blast_snowball-2.0.sql';
comment on column blastn_nonpsr_refseqrna_consensus.qseqid is 'Identifier of the FASTA sequence from which probe set was designed';
comment on column blastn_nonpsr_refseqrna_consensus.refseq_id is 'Best refseq hit for this sequence';
comment on column blastn_nonpsr_refseqrna_consensus.entrez_gene_id is 'entrez id associated to this refseq_id (from NCBI release46.accession2geneid)';
comment on column blastn_nonpsr_refseqrna_consensus.species_name is 'Species from which this refseq comes from';
comment on column blastn_nonpsr_refseqrna_consensus.symbol_ncbi is 'Gene symbol associated to this refseq_id (from NCBI gene_info table)';
comment on column blastn_nonpsr_refseqrna_consensus.symbol_fios is 'Gene symbol associated to this sequence given by FIOS.';
comment on column blastn_nonpsr_refseqrna_consensus.symbol_consensus is 'Gene symbol for this probeset taking either symbol_ncbi, if available, or symbol_fios.';
comment on column blastn_nonpsr_refseqrna_consensus.description is 'Gene description associated to this entrez_gene_id (from NCBI gene_info table)';
create unique index indx_probesetid_consensus on blastn_nonpsr_refseqrna_consensus (probe_set_id);


/*---------------------------------------------------------------------------
  HGNC annotation of SNOWBALL array
  See labbook 18/05/2011
---------------------------------------------------------------------------*/

-- This table from HGNC complete HGNC dataset (http://www.genenames.org/cgi-bin/hgnc_downloads.cgi?title=HGNC+output+data&hgnc_dbtag=on&preset=all&status=Approved&status=Entry+Withdrawn&status_opt=2&level=pri&=on&where=&order_by=gd_app_sym_sort&limit=&format=text&submit=submit&.cgifields=&.cgifields=level&.cgifields=chr&.cgifields=status&.cgifields=hgnc_dbtag)
show client_encoding;
select read_table($$ file:'F:/data/HGNC/hgnc_downloads.txt', table: 'snowball.hgnc_dataset', 
    header:['hgnc_id', 'approved_symbol', 'approved_name', 'status', 'locus_type', 'locus_group', 'previous_symbols', 'previous_names', 
            'synonyms', 'name_synonyms', 'chromosome', 'date_approved', 'date_modified', 'date_symbol_changed', 'date_name_changed', 'accession_numbers', 
            'enzyme_ids', 'entrez_gene_id', 'ensembl_gene_id', 'mouse_genome_database_id', 'specialist_database_links', 'specialist_database_ids', 
            'pubmed_ids', 'refseq_ids', 'gene_family_tag', 'record_type', 'primary_ids', 'secondary_ids', 'ccds_ids', 'vega_ids', 
            'locus_specific_databases', 'gdb_id_mapped', 'entrez_gene_id_mapped', 'omim_id_mapped', 'refseq_id_mapped', 'uniprot_id_mapped', 'ensembl_id_mapped', 
            'ucsc_id', 'mouse_genome_database_id', 'rat_genome_database_id'], skip:1, overwrite:True $$);
comment on table hgnc_dataset is 'HGNC complete dataset from http://www.genenames.org/cgi-bin/hgnc_stats.pl. HGNC update 18/05/2011. See labbook 18/05/2011 and 20110504_blast_snowball-2.0.sql';
select * from hgnc_dataset limit 10;

-- Add a column to the consensus table for the HGNC symbol:
alter table blastn_nonpsr_refseqrna_consensus add column hgnc_approved_symbol text;
alter table blastn_nonpsr_refseqrna_consensus add column hgnc_approved_name text;

-- Populate hgnc column: Match names from symbol_consensus to approved_name case insensitive
update blastn_nonpsr_refseqrna_consensus set hgnc_approved_symbol = hgnc_dataset.approved_symbol, hgnc_approved_name = hgnc_dataset.approved_name
from hgnc_dataset where lower(blastn_nonpsr_refseqrna_consensus.symbol_consensus) = lower(hgnc_dataset.approved_symbol);

-- See if some missing hgnc names can be rescued by querying the fios annotation:
update blastn_nonpsr_refseqrna_consensus set hgnc_approved_symbol = hgnc_dataset.approved_symbol, hgnc_approved_name = hgnc_dataset.approved_name
from hgnc_dataset 
where lower(blastn_nonpsr_refseqrna_consensus.symbol_fios) = lower(hgnc_dataset.approved_symbol) and 
hgnc_approved_symbol is null;

-- Convert synonims column to row:
drop table if exists hgnc_synonyms;
create temp table hgnc_synonyms AS(
    select * from 
        (select distinct approved_symbol, approved_name, (string_to_array(synonyms, ', '))[i] AS synonyms
        from hgnc_dataset, generate_series(1, ( select max(array_upper(string_to_array(synonyms, ', '), 1)) from hgnc_dataset) ) i) as t1
    where synonyms is not null
);


-- Populate remaining hgnc column with consensus names that can be found in the synonyms column:
update blastn_nonpsr_refseqrna_consensus set hgnc_approved_symbol = hgnc_synonyms.approved_symbol, hgnc_approved_name = hgnc_synonyms.approved_name
from hgnc_synonyms
where lower(symbol_consensus) = lower(hgnc_synonyms.synonyms) and
  hgnc_approved_symbol is null;

update blastn_nonpsr_refseqrna_consensus set hgnc_approved_symbol = hgnc_synonyms.approved_symbol, hgnc_approved_name = hgnc_synonyms.approved_name
from hgnc_synonyms
where lower(symbol_fios) = lower(hgnc_synonyms.synonyms) and
  hgnc_approved_symbol is null;

drop table if exists consensus;
create temp table consensus AS(
select probe_set_id, qseqid, refseq_id, entrez_gene_id, species_name, symbol_ncbi, symbol_fios, hgnc_approved_symbol,
    CASE WHEN hgnc_approved_symbol is not null THEN hgnc_approved_symbol                -- 1) if there is annotation from HGNC
         WHEN (symbol_ncbi is not null AND symbol_ncbi not like 'LOC%') THEN symbol_ncbi  -- 2) if NCBI has symbol not starting with LOC
         WHEN symbol_fios is not null THEN symbol_fios                                 -- 3) if fios has annotation (not fios predicted genes start with //)
         WHEN symbol_ncbi is not null THEN symbol_ncbi                                  -- 4) Fill up with LOC... from NCBI
         ELSE NULL END AS symbol_consensus,              
    CASE WHEN hgnc_approved_name is not null THEN hgnc_approved_name
         WHEN description is not null THEN description
         ELSE NULL END AS name_consensus,
    CASE WHEN hgnc_approved_symbol is not null THEN 'HGNC'                                  -- 1) if there is annotation from HGNC
         WHEN (symbol_ncbi is not null AND symbol_ncbi not like 'LOC%') THEN 'dberaldi_ncbi'  -- 2) if NCBI has symbol not starting with LOC
         WHEN symbol_fios is not null THEN 'fios'                                          -- 3) if fios has annotation (not fios predicted genes start with //)
         WHEN symbol_ncbi is not null THEN 'dberaldi_ncbi'                                  -- 4) Fill up with LOC... from NCBI
         ELSE NULL END AS source_consensus
from blastn_nonpsr_refseqrna_consensus
);
drop table blastn_nonpsr_refseqrna_consensus;
create table snowball.blastn_nonpsr_refseqrna_consensus as (select * from consensus);
comment on table blastn_nonpsr_refseqrna_consensus is 'Annotation of the snowball array using annotation in order: HGNC, NCBI, FIOS. See 20110504_blast_snowball-2.0.sql';
comment on table blastn_nonpsr_refseqrna_consensus is 'Annotation of the snowball array using either the results from BLASTN vs refseq (blastn_nonpsr_refseqrna_best) or fios annotation. See 20110504_blast_snowball-2.0.sql';
comment on column blastn_nonpsr_refseqrna_consensus.qseqid is 'Identifier of the FASTA sequence from which probe set was designed';
comment on column blastn_nonpsr_refseqrna_consensus.refseq_id is 'Best refseq hit for this sequence';
comment on column blastn_nonpsr_refseqrna_consensus.entrez_gene_id is 'entrez id associated to this refseq_id (from NCBI release46.accession2geneid)';
comment on column blastn_nonpsr_refseqrna_consensus.species_name is 'Species from which this refseq comes from';
comment on column blastn_nonpsr_refseqrna_consensus.symbol_ncbi is 'Gene symbol associated to this refseq_id (from NCBI gene_info table)';
comment on column blastn_nonpsr_refseqrna_consensus.symbol_fios is 'Gene symbol associated to this sequence given by FIOS.';
comment on column blastn_nonpsr_refseqrna_consensus.symbol_consensus is 'Gene symbol for this probeset taking either symbol_ncbi, if available, or symbol_fios.';
create unique index indx_probesetid_consensus on blastn_nonpsr_refseqrna_consensus (probe_set_id);
select * from blastn_nonpsr_refseqrna_consensus limit 1000;

-- Update qseqid (name of query sequence) to the full name (there is some infomormation in there).
update blastn_nonpsr_refseqrna_consensus set qseqid = non_psr_fastas.sequence_name from non_psr_fastas
where blastn_nonpsr_refseqrna_consensus.probe_set_id = non_psr_fastas.probe_set_id;
alter table blastn_nonpsr_refseqrna_consensus rename column qseqid to sequence_name;


-------------------------------------------------------------------------------
--  EXPORT DATAFILES
-------------------------------------------------------------------------------

copy (select * from blastn_nonpsr_refseqrna_species order by qseqid) to 'F:/data/20110504_snowball_annotation/snowball_annotation_20110523/non_PSR_fastas.blastn.refseqrna-20110523.txt' with csv header delimiter E'\t';
copy(
  select 
    blastn_nonpsr_refseqrna_consensus.probe_set_id,
    blastn_nonpsr_refseqrna_species.qseqid,
    blastn_nonpsr_refseqrna_consensus.sequence_name as non_psr_fastas_sequence_name,
    blastn_nonpsr_refseqrna_consensus.entrez_gene_id,
    blastn_nonpsr_refseqrna_species.refseq_id,
    blastn_nonpsr_refseqrna_species.species_name,
    blastn_nonpsr_refseqrna_species.description,
    blastn_nonpsr_refseqrna_species.type_of_gene,
    blastn_nonpsr_refseqrna_species.pident,
    blastn_nonpsr_refseqrna_species.length,
    blastn_nonpsr_refseqrna_species.mismatch,
    blastn_nonpsr_refseqrna_species.gapopen,
    blastn_nonpsr_refseqrna_species.qstart,
    blastn_nonpsr_refseqrna_species.qend,
    blastn_nonpsr_refseqrna_species.sstart,
    blastn_nonpsr_refseqrna_species.send,
    blastn_nonpsr_refseqrna_species.evalue,
    blastn_nonpsr_refseqrna_species.bitscore,
    blastn_nonpsr_refseqrna_consensus.symbol_ncbi, 
    blastn_nonpsr_refseqrna_consensus.symbol_fios, 
    blastn_nonpsr_refseqrna_consensus.hgnc_approved_symbol, 
    blastn_nonpsr_refseqrna_consensus.symbol_consensus, 
    blastn_nonpsr_refseqrna_consensus.name_consensus, 
    blastn_nonpsr_refseqrna_consensus.source_consensus
  from blastn_nonpsr_refseqrna_consensus left join blastn_nonpsr_refseqrna_species on 
    blastn_nonpsr_refseqrna_consensus.probe_set_id = blastn_nonpsr_refseqrna_species.probe_set_id
  where approved_hit is True or approved_hit is null order by probe_set_id
) to 'F:/data/20110504_snowball_annotation/snowball_annotation_20110523/blastn_nonpsr_refseqrna_consensus-20110523.txt' with csv header delimiter E'\t';

copy (
select CASE WHEN source_consensus is null THEN 'Probes w/o gene symbol' ELSE source_consensus END AS "Annotation", 
       count(*) AS "No. probes",
       count(distinct symbol_consensus) AS "No. distinct genes",
       count(refseq_id) AS "No. RefSeq IDs",
       count(distinct refseq_id) AS "No. distinct RefSeq IDs"
from blastn_nonpsr_refseqrna_consensus group by source_consensus
union select 'All sources',
       count(*) AS "No. genes", 
       count(distinct symbol_consensus) AS "No. distinct genes",
       count(refseq_id) AS "No. RefSeq IDs",
       count(distinct refseq_id) AS "No. distinct RefSeq IDs"
from blastn_nonpsr_refseqrna_consensus
order by "No. probes"
) to 'F:/data/20110504_snowball_annotation/snowball_annotation_20110523/blastn_nonpsr_refseqrna_consensus_summary-20110523.txt' with csv header delimiter E'\t';

/*-----------------------------------------------------------------------------
  Some stats about sequence and annotation
-----------------------------------------------------------------------------*/

-- Lenght and count of query sequences.
select count(*), round(avg(length("sequence"))) as avg, round(stddev(length("sequence"))) as sd, min(length("sequence")) as min, max(length("sequence")) as max
from non_psr_fastas where probe_set_id is not null;

-- Blast hits evalue < 0.001
select count(*) as no_hits, count(distinct blastn_nonpsr_refseqrna_species.probe_set_id) as no_distinct_probes, count(distinct blastn_nonpsr_refseqrna_species.sseqid) no_distinct_subj 
from blastn_nonpsr_refseqrna_species 
where probe_set_id is not null;

-- Blast hits evalue < 1e-9
select count(*) as no_hits, count(distinct blastn_nonpsr_refseqrna_species.probe_set_id) as no_distinct_probes, count(distinct blastn_nonpsr_refseqrna_species.sseqid) no_distinct_subj 
from blastn_nonpsr_refseqrna_species 
where probe_set_id is not null and evalue < 1e-9;

-- Summary stats for the approved probes:
/*In R
library(RODBC); conn<- odbcConnect(dsn= 'pgVitelleschi')
evalue<- sqlQuery(conn, "select probe_set_id, evalue, bitscore from blastn_nonpsr_refseqrna_species where probe_set_id is not null and approved_hit is True;")
*/

-- Number of gene symbols from HGNC and from other sources
select source_consensus, count(probe_set_id) as no_probes, count(distinct symbol_consensus) as no_distinct_symbols from blastn_nonpsr_refseqrna_consensus group by source_consensus;
select count(*) as no_probes, count(distinct symbol_consensus) from blastn_nonpsr_refseqrna_consensus where symbol_consensus is not null;
-- Number of non-annotated probes;
select count(*) - 123 as no_non_annotated_probes from blastn_nonpsr_refseqrna_consensus where refseq_id is null; -- 123 is the number of control probes
select * from blastn_nonpsr_refseqrna_consensus where refseq_id is not null and name_consensus is null and probe_set_id not like 'AFFX-%'; -- Probes with refseq id but no gene symbol

select * from blastn_nonpsr_refseqrna_consensus where probe_set_id not like 'AFFX-%'; -- 47722 Non-control probes  
select * from blastn_nonpsr_refseqrna_consensus where probe_set_id not like 'AFFX-%' and refseq_id is not null order by symbol_consensus desc; -- 40072 Non-control probes with refseq_id (149 w/o gene symbol)
select * from blastn_nonpsr_refseqrna_consensus where probe_set_id not like 'AFFX-%' and refseq_id is null; -- 7650 non-control probes w/o refseq_id
select * from blastn_nonpsr_refseqrna_consensus where probe_set_id not like 'AFFX-%' and symbol_consensus is not null order by refseq_id desc; -- 40022 Non-control probes with gene symbol (99 w/o refseq_id)
select * from blastn_nonpsr_refseqrna_consensus where probe_set_id not like 'AFFX-%' and symbol_consensus is null and refseq_id is null; -- 7551 non-control probes w/o refseq_id and w/o gene symbol
select * from blastn_nonpsr_refseqrna_consensus where probe_set_id not like 'AFFX-%' and symbol_consensus is null and refseq_id is not null; -- 149 non-control probes w/ refseq_id and w/o gene symbol
select * from blastn_nonpsr_refseqrna_consensus where probe_set_id not like 'AFFX-%' and symbol_consensus is not null and refseq_id is null; -- 99 non-control probes w/o refseq_id and w/ gene symbol
select * from blastn_nonpsr_refseqrna_consensus where probe_set_id not like 'AFFX-%' and symbol_consensus is not null or refseq_id is not null; -- 40171 non-control probes w/ refseq_id OR w/ gene symbol
select * from blastn_nonpsr_refseqrna_consensus where probe_set_id not like 'AFFX-%' and source_consensus = 'HGNC'; -- 30640 probes with gene from HGNC 
select distinct(symbol_consensus) from blastn_nonpsr_refseqrna_consensus where probe_set_id not like 'AFFX-%' and source_consensus = 'HGNC'; -- 15772 genes from HGNC

select * from blastn_nonpsr_refseqrna_consensus where probe_set_id not like 'AFFX-%' and symbol_consensus is null and refseq_id is not null
and probe_set_id not in (select probe_set_id from blastn_nonpsr_refseqrna_consensus where probe_set_id not like 'AFFX-%' and symbol_consensus is not null and refseq_id is null);

select probe_set_id from blastn_nonpsr_refseqrna_consensus where probe_set_id not like 'AFFX-%' and symbol_consensus is not null and refseq_id is null
and probe_set_id not in (select probe_set_id from blastn_nonpsr_refseqrna_consensus where probe_set_id not like 'AFFX-%' and symbol_consensus is null and refseq_id is not null);


select * from blastn_nonpsr_refseqrna_consensus where refseq_id is not null and probe_set_id not like 'AFFX-%' and symbol_consensus is null;

select count(source_consensus) from blastn_nonpsr_refseqrna_consensus where symbol_consensus is not null;
select * from blastn_nonpsr_refseqrna_consensus where symbol_consensus is not null and probe_set_id not like 'AFFX-%';

-------------------------------------------------------------------------------
-- TRITUME
-------------------------------------------------------------------------------