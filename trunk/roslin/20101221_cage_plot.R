library(RODBC)
library(moments)
conn<- odbcConnect(dsn= 'pgCAGE')

tss<- sqlFetch(conn, 'tss_to_transcripts')
hist(c(  tss[tss$strand == '+', 'enstss2cluster_end'],
       - tss[tss$strand == '-', 'enstss2cluster_start']), 
     breaks= 40,
     main= "Distance between transcript 5'end (ENSEMBL 59)\nand trascription start clusters (CAGE)",
     xlab= "Distance (bp)\n(Negative dist.: The TSS is upstream to the 5'ends)"
     )
savePlot('M:/Documents/LabBook/LabBook_Figures/20101220_cage_to_transscripts.emf', 'emf')

tc<- sqlQuery(conn, 
  "select * from cage050810 where tss_cluster_id = 'SSTSC_1_R_135267798'")
plot(c(0, tc$read_count, 0), type= 'h')


tc2<- sqlQuery(conn, 
  "select * from cage050810 where tss_cluster_id = 'SSTSC_1_R_294403751'")
plot(c(0, tc2$read_count, 0), type= 'h')




kurtosis(c(0, tc$read_count, 0))


xv<- sqlQuery(conn, "
select distinct rnaseq_ctrl_gtf.transcript_id, 
                promoter_id, 
                fpkm AS transcript_fpkm, 
                promoters.promoter_rpkm 
from rnaseq_ctrl_gtf inner join promoters on 
    rnaseq_ctrl_gtf.transcript_id = promoters.transcript_id
where feature like 'transcript'; 
")

odbcClose(conn)