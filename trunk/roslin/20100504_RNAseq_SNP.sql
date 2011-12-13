
/*
  Import SNP data.
*/
-- truncate pileup_genotype;
-- select read_table($$ file:'C:/Tritume/20100409_RNAseq_CTRL_sscrofa9.56.pileup.snp', dest_table: 'pileup_genotype', append:True,  limit: -1, header:True $$);
-- update pileup_genotype SET dataset_id = '20100409_RNAseq_CTRL_snp' where dataset_id like '20100409_RNAseq_CTRL_tophat_gff';

-- select read_table($$ file:'C:/Tritume/20100409_RNAseq_LPS_sscrofa9.56.pileup.snp', dest_table: 'pileup_genotype', append:True,  limit: -1, header:True $$);
-- update pileup_genotype SET dataset_id = '20100409_RNAseq_LPS_snp' where dataset_id like '20100409_RNAseq_LPS_tophat_gff';

select dataset_id, count(dataset_id) from pileup_genotype group by dataset_id;

/*
  Select SNPs common to both CTRL and LPS in terms of location AND genotype (i.e. exclude SNPs where, for example, ctrl is A/C and lps is A/T).
  Exclude SNPs where the total count of alleles (ctrl + lps, counting only the allels NOT the total nreads) is > n
*/
-- drop table tmp_snp;
create table tmp_snp AS(
select
    ctrl.rname, ctrl.pos, ctrl.rbase, ctrl.allele_1, ctrl.allele_2,   -- SNP location and alleles
    (ctrl.count_allele_1 + ctrl.count_allele_2) AS ctrl_tot,          -- Sum of count allele 1 + allele 2 for ctrl
    (lps.count_allele_1 + lps.count_allele_2) AS lps_tot,             -- Same for LPS
    (ctrl.count_allele_1 + ctrl.count_allele_2 + lps.count_allele_1 + lps.count_allele_2) as tot, -- Allele 1 + Allele 2 in both LPS and CTRL 
    -- Calculate the Major Allele Frequency
    CASE WHEN ctrl.count_allele_1 > ctrl.count_allele_2 THEN ctrl.count_allele_1 / (ctrl.count_allele_1 + ctrl.count_allele_2)::numeric 
                                                        ELSE ctrl.count_allele_2 / (ctrl.count_allele_1 + ctrl.count_allele_2)::numeric END AS ctrl_majf,
    CASE WHEN lps.count_allele_1 > lps.count_allele_2 THEN lps.count_allele_1 / (lps.count_allele_1 + lps.count_allele_2)::numeric 
                                                        ELSE lps.count_allele_2 / (lps.count_allele_1 + lps.count_allele_2)::numeric END AS lps_majf,
    -- Calculate 'Allele 1' frequency. Allele 1 is the allele coming first in alphabetical order
    ctrl.count_allele_1 AS ctrl_a1,
    lps.count_allele_1 AS lps_a1,
    ctrl.count_allele_2 AS ctrl_a2,
    lps.count_allele_2 AS lps_a2,
    ctrl.count_allele_1 / (ctrl.count_allele_1 + ctrl.count_allele_2)::numeric AS ctrl_af,
    lps.count_allele_1 / (lps.count_allele_1 + lps.count_allele_2)::numeric AS lps_af
from 
    (select * from pileup_genotype where dataset_id like '20100409_RNAseq_CTRL_snp') as ctrl 
inner join 
    -- Select SNPs with location common to both libs
    (select * from pileup_genotype where dataset_id like '20100409_RNAseq_LPS_snp') as lps
        on ctrl.rname = lps.rname and ctrl.pos = lps.pos
where
    -- Filter for SNPs where scored alleles are the same in LPS and CTRL
    ctrl.allele_1 = lps.allele_1 AND ctrl.allele_2 = lps.allele_2 
    -- Accept SNPs with at least this many total (lps + ctrl) counts
      AND (ctrl.count_allele_1 + ctrl.count_allele_2 + lps.count_allele_1 + lps.count_allele_2) >= 20
order by rname, pos
);

select * from tmp_snp limit 20;
select * from tmp_snp where rbase not in (allele_1, allele_2);
select count(*) from tmp_snp;


-- SNPs with putative differential allelic expression and with transcription difference
-- Note table: tmp_edger_transcript is equivalent to:
-- select edger_toptags.* from edger_toptags where dataset_id like '20100408_LPSvsCTRL_GTF';
select distinct t1.rname, tmp_edger_transcript.transcript_id, count(t1.transcript_id) AS no_snp, avg(tmp_edger_transcript."logFC") AS "logFC", avg("FDR") AS "FDR", avg(ctrl_af) AS ctrl_afreq, avg(lps_af) AS lps_afreq, t1.gene_name
from cufflinks_transcript_gtf inner join 
    (
    select gtf_attribute(attributes, 'transcript_id') AS transcript_id, 
           gtf_attribute(attributes, 'exon_number') AS exon_number, 
           gtf_attribute(attributes, 'gene_name') AS gene_name, 
           snp_dae.*
    from sus_scrofa_sscrofa9_56_gtf 
    inner join snp_dae on
         sus_scrofa_sscrofa9_56_gtf.rname = snp_dae.rname and
         sus_scrofa_sscrofa9_56_gtf.f_start <= snp_dae.pos and
         sus_scrofa_sscrofa9_56_gtf.f_end >= snp_dae.pos
    ) as t1 
on cufflinks_transcript_gtf.transcript_id = t1.transcript_id
    inner join tmp_edger_transcript on
         tmp_edger_transcript.transcript_id = t1.transcript_id
where feature like 'transcript' and "FDR"< 0.01
group by t1.rname, tmp_edger_transcript.transcript_id, t1.gene_name;

select * from snp_dae;

--------------------------------------------[ Tritume ]--------------------------------------------

select * from tmp_snp where pos = 16429587;
select * from sus_scrofa_sscrofa9_56_gtf where rname like '5' and f_start < 16429587 and f_end > 16429587;

select * from cufflinks_transcript_gtf where rname like '5' and "start" < 16429587 and "end" > 16429587 and feature like 'transcript';
select * from cufflinks_transcript_gtf where transcript_id like 'ENSSSCT00000008324';
select * from sus_scrofa_sscrofa9_56_gtf where rname like '12' and f_start < 50062122 and f_end > 50062122;
select * from cufflinks_transcript_gtf where rname like '12' and "start" < 50062122 and "end" > 50062122 and feature like 'transcript';
50062122


dataset_id, rname, pos, rbase, consensus, consensus_phred, snp_phred, rms, nreads, allele_1, count_allele_1, allele_2, count_allele_2

-- SELECT WHERE nreads > n

---------------------------------------[ Tritume ]------------------------------------
select names($$ t:'pileup_genotype', quote: "" $$)



select * from pileup_genotype limit 10;

select * from pg_settings;

20100409_RNAseq_CTRL_SNP


