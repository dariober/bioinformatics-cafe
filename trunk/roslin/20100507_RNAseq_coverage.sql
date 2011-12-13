create temp table cufflinks_nogtf AS(
    select * from cufflinks_transcript_gtf where source in ('20100317_RNAseq_CTRL_noGTF', '20100317_RNAseq_LPS_noGTF') 
    )
select * from cufflinks_nogtf limit 10;
select source, count(*) from cufflinks_nogtf group by source;


copy (
select distinct 
    ens_gtf.rname, 
    cufflinks_nogtf.source,
    cufflinks_nogtf.transcript_id AS nogtf_transcript_id,
    cufflinks_nogtf.exon_number AS nogtf_exon,
    cufflinks_nogtf.start, 
    cufflinks_nogtf.end, 
    cufflinks_nogtf.strand AS nogft_transcript_id,
    ens_gtf.exon_number AS ens_exon,
    ens_gtf.f_start AS ens_start, 
    ens_gtf.f_end AS ens_end, 
    ens_gtf.strand AS ens_strand
from cufflinks_nogtf inner join (
    select gtf_attribute(attributes, 'exon_number') AS exon_number, * from sus_scrofa_sscrofa9_56_gtf where gtf_attribute(attributes, 'gene_name') like 'ACTB_PIG'
    ) AS ens_gtf on 
               -- Link chromosome to chromosome
               ens_gtf.rname = cufflinks_nogtf.rname AND
            -- de novo transcripts fully contained in Sscrofa GTF
            --   (ens_gtf.f_start <= cufflinks_nogtf.start AND
            --   ens_gtf.f_end >= cufflinks_nogtf.end) OR
            /* 
            De novo transcripts partially contained: denovo-start BEFORE the reference-end and 'denovo-end' AFTER reference-end
                             dn_start |----------------| dn_end
                ens_start |-----------------|ens_end
            */
                (cufflinks_nogtf.start <= ens_gtf.f_end AND
                 cufflinks_nogtf.end > ens_gtf.f_end)   OR
            /* 
            De novo transcripts partially contained: denovo-start BEFORE the reference-end and 'denovo-end' AFTER reference-end
                dn_start |----------------| dn_end
                            ens_start |-----------------|ens_end
            */
                (cufflinks_nogtf.start < ens_gtf.f_start AND
                 cufflinks_nogtf.end >= ens_gtf.f_end)
where ens_gtf.feature like 'exon' and cufflinks_nogtf.feature like 'exon'
) to 'C:/Tritume/denovo.txt' with csv header delimiter E'\t';

where "start"
ENSSSCT00000001533

select source, strand, count(strand) from cufflinks_nogtf group by strand, source;

select * from sus_scrofa_sscrofa9_56_gtf where attributes like '%ACTB_PIG%';
" gene_id "ENSSSCG00000007585"; transcript_id "ENSSSCT00000008324"; exon_number "6"; gene_name "ACTB_PIG"; transcript_name "ACTB_PIG"; protein_id "ENSSSCP00000008105";"

select gtf_attribute($$" gene_id "ENSSSCG00000001404"; transcript_id "ENSSSCT00000001533"; exon_number "1"; gene_name "TNFA_PIG"; transcript_name "TNFA_PIG";"$$, 'gene_name');



select * from cuffcompare_combined_gtf where v9 like '%ENSSSCT00000008324%';