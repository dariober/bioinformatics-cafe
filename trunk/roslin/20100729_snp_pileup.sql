-- select read_table($$ file:'D:/Tritume/20100728_LPS_snp.pileup', dataset_id: '20100728_LPS', dest_table: 'pileup_genotype', append:True  $$);
-- select read_table($$ file:'D:/Tritume/20100728_CTRL_snp.pileup', dataset_id: '20100728_CTRL', dest_table: 'pileup_genotype', append:True  $$);

select * from pileup_genotype limit 10;

create temp table snp_ctrl AS(
    -- SNPs in CTRL library
    select rname || '_' || pos || '_' || consensus AS snp_id, * from pileup_genotype where dataset_id like '20100728_CTRL'
    );

create temp table snp_lps AS(
    -- SNPs in LPS library
    select rname || '_' || pos || '_' || consensus AS snp_id, * from pileup_genotype where dataset_id like '20100728_LPS'
    );

create temp table snp_id AS(
    -- All SNPs in CTRL and LPS
    select distinct rname || '_' || pos || '_' || consensus AS snp_id from pileup_genotype where dataset_id in ('20100728_LPS', '20100728_CTRL')
    );


-- SNPs in both libraries
create 
select snp_ctrl.snp_id,
       snp_ctrl.allele_1,
       snp_ctrl.allele_2,
       snp_ctrl.rbase,
       snp_ctrl.count_allele_1 AS ctrl_a1,
       snp_ctrl.count_allele_2 AS ctrl_a2,
       -- Frequency of reference allele in CTRL
       CASE WHEN snp_ctrl.rbase = snp_ctrl.allele_1 THEN snp_ctrl.count_allele_1 / (snp_ctrl.count_allele_1 + snp_ctrl.count_allele_2)::numeric 
            WHEN snp_ctrl.rbase = snp_ctrl.allele_2 THEN snp_ctrl.count_allele_2 / (snp_ctrl.count_allele_1 + snp_ctrl.count_allele_2)::numeric END AS ctrl_ref_freq,
       snp_lps.count_allele_1 AS lps_a1,
       snp_lps.count_allele_2 AS lps_a2,
       -- Frequency of reference allele in LPS
       CASE WHEN snp_lps.rbase = snp_lps.allele_1 THEN snp_lps.count_allele_1 / (snp_lps.count_allele_1 + snp_lps.count_allele_2)::numeric 
            WHEN snp_lps.rbase = snp_lps.allele_2 THEN snp_lps.count_allele_2 / (snp_lps.count_allele_1 + snp_lps.count_allele_2)::numeric END AS lps_ref_freq
/*
       -- Total count of reference base (CTRL + LPS)
       CASE WHEN (snp_ctrl.rbase = snp_ctrl.allele_1) THEN snp_ctrl.count_allele_1 + snp_lps.count_allele_1
            WHEN (snp_ctrl.rbase = snp_ctrl.allele_2) THEN snp_ctrl.count_allele_2 + snp_lps.count_allele_2 END AS ref_count,
       -- Total count of alternative base (CTRL + LPS)
       CASE WHEN (snp_ctrl.rbase != snp_ctrl.allele_1) THEN snp_ctrl.count_allele_1 + snp_lps.count_allele_1
            WHEN (snp_ctrl.rbase != snp_ctrl.allele_2) THEN snp_ctrl.count_allele_2 + snp_lps.count_allele_2 END AS alt_count
*/
from snp_ctrl left join snp_lps on snp_ctrl.snp_id = snp_lps.snp_id
-- Remove SNPs counted <20 times in  either library
where snp_lps.count_allele_1 + snp_lps.count_allele_2 >= 20 and
      snp_ctrl.count_allele_1 + snp_ctrl.count_allele_2 >= 20 and
      (snp_ctrl.rbase = snp_ctrl.allele_1 or snp_ctrl.rbase = snp_ctrl.allele_2)
order by snp_id;

-----------------------------------[ Tritume ]---------------------------------

select snp_id.snp_id,
       snp_ctrl.allele_1,
       snp_ctrl.allele_2,
       snp_ctrl.count_allele_1 AS ctrl_a1,
       snp_ctrl.count_allele_2 AS ctrl_a2,
       snp_lps.count_allele_1 AS lps_a1,
       snp_lps.count_allele_2 AS lps_a2
from snp_id left join snp_ctrl on snp_id.snp_id = snp_ctrl.snp_id
            left join snp_lps on snp_id.snp_id = snp_lps.snp_id
-- Remove SNPs counted <20 times in  either library
where snp_lps.count_allele_1 + snp_lps.count_allele_2 >= 20 and
      snp_ctrl.count_allele_1 + snp_ctrl.count_allele_2 >= 20
order by snp_id
limit 100;

select names('pileup_genotype');
"dataset_id", "rname", "pos", "rbase", "consensus", "consensus_phred", "snp_phred", "rms", "nreads", "allele_1", "count_allele_1", "allele_2", "count_allele_2"
