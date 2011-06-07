/* Import cuffcompare/cuffdiff output */ 

-- GTF file from cuffcompare:
select read_table($$
    file: 'F:/data/20110202_cufflinks/20110202_RNAseq_am_bmdm/am_bmdm_gtf/am_bmdm_gtf.combined.gtf',
    table: 'cufflinks.combined_gtf_am_bmdm_gtf',
    header: ['rname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'gene_id', 'transcript_id', 'exon_number', 'gene_name', 'oid', 'nearest_ref', 'class_code', 'tss_id', 'p_id'],
    overwrite: False,
    apply:""" ##
def apply(line):
        source= 'am_bmdm_gtf'
        gtf_line= ['']*17
        line= line.split('	')
        gtf_line[0:8]= line[0:8]
        attribute= line[8].split(';')
        for attr in attribute:
            if attr.startswith('gene_id "'):
                gtf_line[8]= attr[9:-1]
            elif attr.startswith(' transcript_id "'):
                gtf_line[9]= attr[16:-1]
            elif attr.startswith(' exon_number "'):
                gtf_line[10]= attr[14:-1]
            elif attr.startswith(' gene_name "'):
                gtf_line[11]= attr[12:-1]
            elif attr.startswith(' oId "'):
                gtf_line[12]= attr[6:-1]
            elif attr.startswith(' nearest_ref "'):
                gtf_line[13]= attr[14:-1]
            elif attr.startswith(' class_code "'):
                gtf_line[14]= attr[13:-1]
            elif attr.startswith(' tss_id "'):
                gtf_line[15]= attr[9:-1]            
            elif attr.startswith(' p_id "'):
                gtf_line[16]= attr[7:-1]
            elif attr == '':
                pass
            else:
                sys.exit('Unrecognized GTF attribute: ' + attr)
        return('	'.join(gtf_line))
"""
$$); --'
COMMENT ON TABLE combined_gtf_am_bmdm_gtf IS 'Output of cuffcompare/cuffdiff pipeline from transcripts assembled on the basis of the ensembl GTF annotation. See /exports/.../cufflinks/output/20110202_RNAseq_am_bmdm/am_bmdm_gtf and 20110501_import_cuffcompare.sql';
select * from combined_gtf_am_bmdm_gtf limit 10;
create index ind_gene_id on combined_gtf_am_bmdm_gtf (gene_id);
vacuum analyze combined_gtf_am_bmdm_gtf;

-- Gene tracking
select read_table($$
    file: 'F:/data/20110202_cufflinks/20110202_RNAseq_am_bmdm/am_bmdm_gtf/genes.fpkm_tracking',
    table: 'cufflinks.genes_fpkm_tracking_am_bmdm_gtf',
    header: True,
    overwrite: True
$$); --'
COMMENT ON TABLE genes_fpkm_tracking_am_bmdm_gtf IS 'Output of cuffcompare/cuffdiff pipeline from transcripts assembled on the basis of the ensembl GTF annotation. See /exports/.../cufflinks/output/20110202_RNAseq_am_bmdm/am_bmdm_gtf and 20110501_import_cuffcompare.sql';
select * from genes_fpkm_tracking_am_bmdm_gtf limit 10;

-- Gene differential expression
select read_table($$
    file: 'F:/data/20110202_cufflinks/20110202_RNAseq_am_bmdm/am_bmdm_gtf/gene_exp.diff',
    table: 'cufflinks.gene_exp_diff_am_bmdm_gtf',
    header: True,
    overwrite: True
$$);
COMMENT ON TABLE gene_exp_diff_am_bmdm_gtf IS 'Output of cuffcompare/cuffdiff pipeline from transcripts assembled on the basis of the ensembl GTF annotation. See /exports/.../cufflinks/output/20110202_RNAseq_am_bmdm/am_bmdm_gtf and 20110501_import_cuffcompare.sql';
select * from gene_exp_diff_am_bmdm_gtf limit 10;
create index ind_genediff_test_id on gene_exp_diff_am_bmdm_gtf (test_id);
vacuum analyze gene_exp_diff_am_bmdm_gtf;

-- Promoter differential expression
select read_table($$
    file: 'F:/data/20110202_cufflinks/20110202_RNAseq_am_bmdm/am_bmdm_gtf/promoters.diff',
    table: 'cufflinks.promoters_diff_am_bmdm_gtf',
    header: True,
    overwrite: True
$$);
COMMENT ON TABLE promoters_diff_am_bmdm_gtf IS 'Output of cuffcompare/cuffdiff pipeline from transcripts assembled on the basis of the ensembl GTF annotation. See /exports/.../cufflinks/output/20110202_RNAseq_am_bmdm/am_bmdm_gtf and 20110501_import_cuffcompare.sql';
select * from promoters_diff_am_bmdm_gtf limit 10;
vacuum analyze promoters_diff_am_bmdm_gtf;

-- Differential Splicing
select read_table($$
    file: 'F:/data/20110202_cufflinks/20110202_RNAseq_am_bmdm/am_bmdm_gtf/splicing.diff',
    table: 'cufflinks.splicing_diff_am_bmdm_gtf',
    header: True,
    overwrite: True
$$);
COMMENT ON TABLE splicing_diff_am_bmdm_gtf IS 'Output of cuffcompare/cuffdiff pipeline from transcripts assembled on the basis of the ensembl GTF annotation. See /exports/.../cufflinks/output/20110202_RNAseq_am_bmdm/am_bmdm_gtf and 20110501_import_cuffcompare.sql';
select * from splicing_diff_am_bmdm_gtf where significant = 'yes';


-- Isoform tracking
select read_table($$
    file: 'F:/data/20110202_cufflinks/20110202_RNAseq_am_bmdm/am_bmdm_gtf/isoforms.fpkm_tracking',
    table: 'cufflinks.isoforms_fpkm_tracking_am_bmdm_gtf',
    header: True,
    overwrite: True
$$); --'
COMMENT ON TABLE isoforms_fpkm_tracking_am_bmdm_gtf IS 'Output of cuffcompare/cuffdiff pipeline from transcripts assembled on the basis of the ensembl GTF annotation. See /exports/.../cufflinks/output/20110202_RNAseq_am_bmdm/am_bmdm_gtf and 20110501_import_cuffcompare.sql';
select * from isoforms_fpkm_tracking_am_bmdm_gtf limit 10;

-- Isoform differenatial expression
select read_table($$
    file: 'F:/data/20110202_cufflinks/20110202_RNAseq_am_bmdm/am_bmdm_gtf/isoform_exp.diff',
    table: 'cufflinks.isoform_exp_diff_am_bmdm_gtf',
    header: True,
    overwrite: True
$$); --'
COMMENT ON TABLE isoform_exp_diff_am_bmdm_gtf IS 'Output of cuffcompare/cuffdiff pipeline from transcripts assembled on the basis of the ensembl GTF annotation. See /exports/.../cufflinks/output/20110202_RNAseq_am_bmdm/am_bmdm_gtf and 20110501_import_cuffcompare.sql';
select * from isoform_exp_diff_am_bmdm_gtf limit 10;


-- TSS tracking
select read_table($$
    file: 'F:/data/20110202_cufflinks/20110202_RNAseq_am_bmdm/am_bmdm_gtf/tss_groups.fpkm_tracking',
    table: 'cufflinks.tss_groups_fpkm_tracking_am_bmdm_gtf',
    header: True,
    overwrite: True
$$); --
COMMENT ON TABLE tss_groups_fpkm_tracking_am_bmdm_gtf IS 'Output of cuffcompare/cuffdiff pipeline from transcripts assembled on the basis of the ensembl GTF annotation. See /exports/.../cufflinks/output/20110202_RNAseq_am_bmdm/am_bmdm_gtf and 20110501_import_cuffcompare.sql';
select * from tss_groups_fpkm_tracking_am_bmdm_gtf limit 10;

-- TSS differenatial expression
select read_table($$
    file: 'F:/data/20110202_cufflinks/20110202_RNAseq_am_bmdm/am_bmdm_gtf/tss_group_exp.diff',
    table: 'cufflinks.tss_group_exp_diff_am_bmdm_gtf',
    header: True,
    overwrite: True
$$); --'
COMMENT ON TABLE tss_group_exp_diff_am_bmdm_gtf IS 'Output of cuffcompare/cuffdiff pipeline from transcripts assembled on the basis of the ensembl GTF annotation. See /exports/.../cufflinks/output/20110202_RNAseq_am_bmdm/am_bmdm_gtf and 20110501_import_cuffcompare.sql';
select * from tss_group_exp_diff_am_bmdm_gtf limit 10;

-------------------------------------------------------------------------------
-- END OF IMPORTING
-------------------------------------------------------------------------------

-- Differential splicing
select 
    * -- sample_1, sample_2, count(*) 
from 
    splicing_diff_am_bmdm_gtf 
where  
    (significant = 'yes' and sample_1 = 'am_ctrl' and sample_2 = 'am_lps') OR
    (significant = 'yes' and sample_1 = 'bmdm_ctrl' and sample_2 = 'bmdm_lps')
/*group by 
    sample_1, sample_2*/
order by 
    gene;

-- Differential promoter usage
select 
    sample_1, sample_2, *
from 
    promoters_diff_am_bmdm_gtf
where  
    (significant = 'yes' and sample_1 = 'am_ctrl' and sample_2 = 'am_lps') OR
    (significant = 'yes' and sample_1 = 'bmdm_ctrl' and sample_2 = 'bmdm_lps')
group by 
    sample_1, sample_2;

-- Differential TSS usage
select 
    sample_1, sample_2, count(*)
from 
    tss_group_exp_diff_am_bmdm_gtf
where  
    (significant = 'yes' and sample_1 = 'am_ctrl' and sample_2 = 'am_lps') OR
    (significant = 'yes' and sample_1 = 'bmdm_ctrl' and sample_2 = 'bmdm_lps')
group by 
    sample_1, sample_2;


-- Genes differentially expressed (FDR < 0.05)
drop table if exists genes_de;
create temp table genes_de AS(
    select distinct *, log2((value_2 + 1e-6) / (value_1 + 1e-6)) AS log_fc,
        CASE WHEN sample_1 = 'am_ctrl' and significant = 'yes' THEN '1' ELSE 0 END AS am_sign,
        CASE WHEN sample_1 = 'bmdm_ctrl' and significant = 'yes' THEN '1' ELSE 0 END AS bmdm_sign 
    from gene_exp_diff_am_bmdm_gtf 
    where significant = 'yes' and gene not like '-' and
        ((sample_1 = 'am_ctrl' and sample_2 = 'am_lps') OR (sample_1 = 'bmdm_ctrl' and sample_2 = 'bmdm_lps'))
    order by gene
);
select * from genes_de order by gene limit 10;

/* ----------------------------------------------------------------------------
   Comparison cuffdiff / affymetrix

   How many genes in common/different between the two platforms?
 --------------------------------------------------------------------------- */

-- Consensus table: ensembl transcript IDs in bith affymetrix and RNAseq platforms:
drop table if exists cons_trans;
create temp table cons_trans AS(
    select distinct combined_gtf_am_bmdm_gtf.nearest_ref as ensembl_transcript_id
    -- Transcripts tested for D.E. by cuffdiff
    from isoform_exp_diff_am_bmdm_gtf inner join combined_gtf_am_bmdm_gtf on 
        isoform_exp_diff_am_bmdm_gtf.test_id= combined_gtf_am_bmdm_gtf.transcript_id
    where isoform_exp_diff_am_bmdm_gtf.sample_1 = 'bmdm_ctrl' and isoform_exp_diff_am_bmdm_gtf.sample_2 = 'bmdm_lps' -- 21980
    intersect 
       -- transcripts available in affymetrix annotation
       (select distinct ensembl_transcript_id from annotation_ctuggle_ensembl inner join annotation_ctuggle on annotation_ctuggle_ensembl.gene_name = annotation_ctuggle.gene_name)-- 7792
);-- intersection: 7825
select count(*) from cons_trans limit 10;

-- Alternative link: Connect by gene_name instead of ensembl_transcript_id. However, less genes in common are returned.
/*
create temp table cons_trans AS(
    select distinct upper(gene_exp_diff_am_bmdm_gtf.gene) -- Remove possible differences due to case.
    -- Transcripts tested for D.E. by cuffdiff
    from gene_exp_diff_am_bmdm_gtf 
    where gene_exp_diff_am_bmdm_gtf.sample_1 = 'bmdm_ctrl' and gene_exp_diff_am_bmdm_gtf.sample_2 = 'bmdm_lps' -- 12019
    intersect 
       -- transcripts available in affymetrix annotation
       select distinct upper(gene_name) from annotation_ctuggle where gene_name is not null -- 11373
);-- intersection: 7177
*/

-- Probes and genes at 0 vs 7 h:
drop table affy_de;
create temp table affy_de AS(
  select ensembl_transcript_id, annotation_ctuggle_ensembl.gene_name, avg(limma_toptable.logfc) as avg_log_fc, avg(adjpval) as avg_fdr
  from limma_toptable inner join annotation_ctuggle on 
      annotation_ctuggle.probe_set_id = limma_toptable.id
                    inner join annotation_ctuggle_ensembl on 
      annotation_ctuggle_ensembl.gene_name = annotation_ctuggle.gene_name
  where coef = '7vs0'
  group by ensembl_transcript_id, annotation_ctuggle_ensembl.gene_name
);
select * from affy_de;
select count(distinct ensembl_transcript_id) from affy_de; -- 7792 (annotated) transcripts.

-- transcripts from rnaseq:
drop table if exists rnaseq_de;
create temp table rnaseq_de AS(
    select distinct combined_gtf_am_bmdm_gtf.nearest_ref, isoform_exp_diff_am_bmdm_gtf.*, log2((value_2 + 1e-3)/(value_1 + 1e-3)) as log_fc 
    from isoform_exp_diff_am_bmdm_gtf inner join combined_gtf_am_bmdm_gtf on 
        isoform_exp_diff_am_bmdm_gtf.test_id = combined_gtf_am_bmdm_gtf.transcript_id
    where sample_1 = 'bmdm_ctrl' and sample_2 = 'bmdm_lps'
); -- 21980
select * from rnaseq_de limit 10;

-- For transcripts found in both platforms assign differential expression from limma and from cuffdiff:
create table public.platforms AS(
select cons_trans.ensembl_transcript_id, 
       affy_de.avg_log_fc as affy_log_fc, 
       affy_de.avg_fdr,
       value_1 as rnaseq_ctrl, 
       value_2 as rnaseq_lps, 
       log_fc as rnaseq_log_fc, 
       rnaseq_de.p_value, rnaseq_de.significant
from cons_trans inner join affy_de   on cons_trans.ensembl_transcript_id = affy_de.ensembl_transcript_id
                inner join rnaseq_de on cons_trans.ensembl_transcript_id = rnaseq_de.nearest_ref
order by ensembl_transcript_id
);
comment on table public.platforms is 'Comparison in log fold change and p-value between RNAseq and affymetrix. Using output of limma_toptable and output of cuffdiff (isoforms level). See 20110501_import_cuffcompare.sql'

-------------------------------------------------------------------------------
drop table if exists public.cuffcompare_genes_de;
create table public.cuffcompare_genes_de AS(
  select distinct genes_fpkm_tracking_am_bmdm_gtf.*, am_sign, bmdm_sign,
      log2(("am_lps_FPKM"+1e-3)/("am_ctrl_FPKM"+1e-3)) as log_fc_am, log2(("bmdm_lps_FPKM"+1e-3)/("bmdm_ctrl_FPKM"+1e-3)) as log_fc_bmdm
  from (select distinct test_id, sum(am_sign) as am_sign, sum(bmdm_sign) as bmdm_sign from genes_de group by test_id) AS genes_de inner join genes_fpkm_tracking_am_bmdm_gtf on 
      genes_fpkm_tracking_am_bmdm_gtf.tracking_id = genes_de.test_id
  order by log_fc_am
);
comment on table public.cuffcompare_genes_de is 'Temporary table prepared with 20110501_import_cuffcompare.sql. Genes identified as D.E. by cuffdiff in am ctrl vs am lps or bmdm ctrl vs bmdm lps';

-- Differential splicing
select distinct nearest_ref, splicing_diff_am_bmdm_gtf.*, isoforms_fpkm_tracking_am_bmdm_gtf.* 
from splicing_diff_am_bmdm_gtf 
    inner join combined_gtf_am_bmdm_gtf on splicing_diff_am_bmdm_gtf.test_id = combined_gtf_am_bmdm_gtf.tss_id
    inner join isoforms_fpkm_tracking_am_bmdm_gtf on isoforms_fpkm_tracking_am_bmdm_gtf.tss_id = splicing_diff_am_bmdm_gtf.test_id
where (sample_1 = 'am_ctrl' and sample_2 = 'am_lps' and significant = 'yes') or
      (sample_1 = 'bmdm_ctrl' and sample_2 = 'bmdm_lps' and significant = 'yes')
order by nearest_ref;
-------------------------------------------------------------------------------
-- TRITUME
-------------------------------------------------------------------------------

select distinct tss_id, test_id, gene_name, gene from combined_gtf_am_bmdm_gtf inner join splicing_diff_am_bmdm_gtf on tss_id = test_id



-- Transcripts
create temp table rnaseq_de AS(
  select distinct genes_de.*, combined_gtf_am_bmdm_gtf.nearest_ref
  from genes_de inner join combined_gtf_am_bmdm_gtf on 
      genes_de.test_id= combined_gtf_am_bmdm_gtf.gene_id
  where genes_de.sample_1 = 'bmdm_ctrl' and genes_de.sample_2 = 'bmdm_lps'
);


-- Assign to cuffdiff transcripts the ensembl_transcript_id:
create temp table rnaseq_de AS(
  select distinct genes_de.*, combined_gtf_am_bmdm_gtf.nearest_ref
  from genes_de inner join combined_gtf_am_bmdm_gtf on 
      genes_de.test_id= combined_gtf_am_bmdm_gtf.gene_id
  where genes_de.sample_1 = 'bmdm_ctrl' and genes_de.sample_2 = 'bmdm_lps'
);
select count(distinct nearest_ref) from rnaseq_de; -- 2827

-- How many of the RNAseq transcripts (3110) are annotated in affymetrix chip?
select count(distinct rnaseq_de.nearest_ref)
from rnaseq_de inner join annotation_ctuggle_ensembl on 
    annotation_ctuggle_ensembl.ensembl_transcript_id = rnaseq_de.nearest_ref; -- 2162

-- Transcripts found both by affymetrix and RNAseq:
select affy_de.ensembl_transcript_id 
from rnaseq_de inner join affy_de on 
    affy_de.ensembl_transcript_id = rnaseq_de.nearest_ref; -- 1249

-- Affymetrix DE, not :
select affy_de.ensembl_transcript_id 
from rnaseq_de inner join affy_de on 
    affy_de.ensembl_transcript_id = rnaseq_de.nearest_ref; -- 1249

-- GTF file from cuffcompare:
select * from combined_gtf_am_bmdm_gtf where gene_name = 'U4' limit 10;

-- Gene tracking
select * from genes_fpkm_tracking_am_bmdm_gtf order by "am_ctrl_FPKM" where tracking_id = 'XLOC_011634' limit 1000;

-- Gene differential expression
select * from gene_exp_diff_am_bmdm_gtf where test_id = 'XLOC_011634' limit 10;

-- Differential Splicing
select * from splicing_diff_am_bmdm_gtf limit 10;



select * from combined_gtf_am_bmdm_gtf where tss_id = 'TSS100';

select * from limma_toptable where id like 'Ssc.3706.1.%'

select nullif()