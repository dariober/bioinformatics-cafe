/* Importing cufflinks output:
- Alveolar macrophages CTRL and LPS
- BMDM CTRL and LPS
From: 
/exports/work/vet_roslin_nextgen/dario/tophat/output/20110202_rnaseq_am_ctrl
/exports/work/vet_roslin_nextgen/dario/tophat/output/20110202_rnaseq_am_lps
/exports/work/vet_roslin_nextgen/dario/tophat/output/20110202_rnaseq_bmdm_lps
/exports/work/vet_roslin_nextgen/dario/tophat/output/20110202_rnaseq_bmdm_lps
*/
create table cufflinks_transcript_gtf AS (select * from cufflinks_transcript_gtf_old where 1=2);
-- Transfer files
select read_table($$ file: 'D:/Tritume/transcripts_gtf_am_ctrl.gtf', table: 'cufflinks_transcript_gtf', append:True $$)
update