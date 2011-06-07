/*
   20110519_snowball_csf1r.sql
   Why is CSF1R not in the snowball annotation?
*/

select names($$ t:'blastn_nonpsr_refseqrna', quote: "'"  $$)
select read_table($$ file:'F:/data/20110504_snowball_annotation/ENSSSCG00000014441_CSF1R.blastnout', table:'blastn_csf1r_nonpsr', temp:True, overwrite:True,
    header: ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'] $$)
comment on table blastn_csf1r_nonpsr is 'BLAST output from CSF1R cDNA from pig ENSSSCG00000014441 against non_PSR_fastas. See labbook 15/08/2011 and 20110519_snowball_csf1r.sql';
select distinct sseqid from blastn_csf1r_nonpsr where evalue < 1e-9;

select distinct blastn_nonpsr_refseqrna_species.* from blastn_csf1r_nonpsr inner join blastn_nonpsr_refseqrna_species on blastn_csf1r_nonpsr.sseqid = blastn_nonpsr_refseqrna_species.qseqid;
copy (
    select distinct blastn_nonpsr_refseqrna_best.* from blastn_csf1r_nonpsr inner join blastn_nonpsr_refseqrna_best on blastn_csf1r_nonpsr.sseqid = blastn_nonpsr_refseqrna_best.qseqid
) to 'D:/Tritume/csf1r_putative_hits.txt' with csv header delimiter E'\t';

select 

-- Blast hits of consensus annotation;
create temp table putative_csf1r AS(
select * from (
select distinct 
    blastn_nonpsr_refseqrna_consensus.refseq_id,
    blastn_nonpsr_refseqrna_consensus.species_name,
    blastn_nonpsr_refseqrna_consensus.symbol_consensus,
    blastn_nonpsr_refseqrna_consensus.name_consensus,
    qseqid,
    pident,
    length,
    evalue,
    bitscore,
    row_number() OVER (PARTITION BY qseqid order by bitscore desc, pident desc) as hit_rank
from blastn_nonpsr_refseqrna_consensus, blastn_nonpsr_refseqrna_species
where blastn_nonpsr_refseqrna_consensus.sequence_name like blastn_nonpsr_refseqrna_species.qseqid || ' %' and
 blastn_nonpsr_refseqrna_consensus.refseq_id = blastn_nonpsr_refseqrna_species.refseq_id and 
 blastn_nonpsr_refseqrna_species.qseqid in (select distinct sseqid from blastn_csf1r_nonpsr) ) as t1
where hit_rank = 1
order by qseqid
);
select * from putative_csf1r;

create temp table csf1r_hits_best AS(
select * from 
(select *, row_number() OVER (PARTITION BY sseqid order by bitscore desc, pident desc) as hit_rank from blastn_csf1r_nonpsr) as t1 
where hit_rank = 1
order by sseqid
);
select * from csf1r_hits_best;

copy (
select 
    putative_csf1r.*, 
    csf1r_hits_best.pident AS pident_csf1r, 
    csf1r_hits_best.length AS length_csf1r, 
    csf1r_hits_best.evalue AS evalue_csf1r, 
    csf1r_hits_best.bitscore AS bitscore_csf1r,
    putative_csf1r.bitscore - csf1r_hits_best.bitscore AS bitscore_ncbi_vs_csf1r
from csf1r_hits_best inner join putative_csf1r on putative_csf1r.qseqid = csf1r_hits_best.sseqid
) TO 'F:/data/20110504_snowball_annotation/csf1r_putative_hits.txt' with csv header delimiter E'\t';

