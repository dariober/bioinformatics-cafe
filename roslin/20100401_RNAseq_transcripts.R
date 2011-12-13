library(RODBC)
conn_db<- odbcConnect(dsn= "pgCAGE", case= "nochange", uid= "dberaldi", pwd= "mypassword")
gtf<- sqlQuery(conn_db, "select * from cufflinks_transcript_gtf where feature like 'transcript'")

gtf_hist<- hist(log(gtf$fpkm[gtf$source == '20100317_RNAseq_CTRL_noGTF']), main= "Distribution of FPKM", xlab= 'Log(FPKM)', border= 'blue')
gtf_hist<- hist(log(gtf$fpkm[gtf$source == '20100317_RNAseq_LPS_noGTF']), main= "Distribution of FPKM", xlab= 'Log(FPKM)', border= 'red', add=TRUE)
abline(v= c(log(2), log(5)), col= "lightgrey", lty= "dashed")
legend('topright', legend= c('CTRL', 'LPS', 'FPKM 2 and 5'), col= c('blue', 'red', 'lightgrey'), lty= 'solid', lwd= 2)

gtf_hist<- hist(log(gtf$fpkm[gtf$source == '20100317_RNAseq_CTRL']), 
  ylim= c(0, 18000), border= "blue", breaks= 40, main= "Distribution of FPKM", xlab= 'Log(FPKM)')
nogtf_hist<- hist(log(gtf$fpkm[gtf$source == '20100317_RNAseq_CTRL_noGTF']), add= T, border= "red")
abline(v= c(log(3), log(5)), col= "lightgrey", lty= "dashed")
legend(legend= c())

par(mfrow= c(4,1), mar= c(5, 4.1, 2.5, 1))
# -------------------[ ACTB ]-------------------
ref<- sqlQuery(conn_db, "select * from sus_scrofa_sscrofa9_56_gtf where attributes like '%ENSSSCT00000008324%' order by f_start")
lps<- sqlQuery(conn_db, "select * from cufflinks_transcript_gtf where transcript_id like 'LPS.115774.0' order by start")
ctrl<-sqlQuery(conn_db, "select * from cufflinks_transcript_gtf where transcript_id like 'CTRL.122412.0' order by start")
ctrl2<-sqlQuery(conn_db, "select * from cufflinks_transcript_gtf where transcript_id like 'CTRL.122415.0' order by start")
lim1<- min(c(ref$f_start, lps$start, ctrl$start))
lim2<- max(c(ref$f_end, lps$end, ctrl$end))

main= paste('ACTB (', round(mean(c(lps$fpkm, ctrl$fpkm)),0), ' FPKM)', sep= '')
plot(x= 1, y=1, xlim= c(lim1, lim2), type= 'n', xlab= '', ylab= '', axes= F, ylim=c(0,3), main= main)
axis(side=1)
ref_ex<- ref[ref$feature == 'exon',]
segments(x0= ref_ex$f_start, x1= ref_ex$f_end, y0= 1, y1= 1, lwd= 5, col= 'red')
lps_ex<- lps[lps$feature == 'exon',]
segments(x0= lps_ex$start, x1= lps_ex$end, y0= 2, y1= 2, lwd= 5, col= 'lightblue')
ctrl_ex<- ctrl[ctrl$feature == 'exon',]
segments(x0= ctrl_ex$start, x1= ctrl_ex$end, y0= 3, y1= 3, lwd= 5, col= 'lightblue')
ctrl2_ex<- ctrl2[ctrl2$feature == 'exon',]
segments(x0= ctrl2_ex$start, x1= ctrl2_ex$end, y0= 3, y1= 3, lwd= 5, col= 'darkblue')
mtext(side= 2, at= 1:3, text= c('ensembl', 'LPS', 'CTRL'), adj= 0.5, las= 1, cex= 0.75)

# --------------------[ HPRT ]----------------------
ref<- sqlQuery(conn_db, "select * from sus_scrofa_sscrofa9_56_gtf where attributes like '%ENSSSCT00000013865%' order by f_start")
lps<- sqlQuery(conn_db, "select * from cufflinks_transcript_gtf where transcript_id like 'LPS.211433.0' order by start")
ctrl<-sqlQuery(conn_db, "select * from cufflinks_transcript_gtf where transcript_id like 'CTRL.222634.0' order by start")
lim1<- min(c(ref$f_start, lps$start, ctrl$start))
lim2<- max(c(ref$f_end, lps$end, ctrl$end))

main= paste('HPRT (', round(mean(c(lps$fpkm, ctrl$fpkm)),0), ' FPKM)', sep= '')
plot(x= 1, y=1, xlim= c(lim1, lim2), type= 'n', xlab= '', ylab= '', axes= F, ylim=c(0,3), main=main)
axis(side=1)
ref_ex<- ref[ref$feature == 'exon',]
segments(x0= ref_ex$f_start, x1= ref_ex$f_end, y0= 1, y1= 1, lwd= 5, col= 'red')
lps_ex<- lps[lps$feature == 'exon',]
segments(x0= lps_ex$start, x1= lps_ex$end, y0= 2, y1= 2, lwd= 5, col= 'lightblue')
ctrl_ex<- ctrl[ctrl$feature == 'exon',]
segments(x0= ctrl_ex$start, x1= ctrl_ex$end, y0= 3, y1= 3, lwd= 5, col= 'lightblue')
ctrl2_ex<- ctrl2[ctrl2$feature == 'exon',]
segments(x0= ctrl2_ex$start, x1= ctrl2_ex$end, y0= 3, y1= 3, lwd= 5, col= 'darkblue')

mtext(side= 2, at= 1:3, text= c('ensembl', 'LPS', 'CTRL'), adj= 0.5, las= 1, cex= 0.75)

# -----------------------[ GAPDH ]-----------------------
ref<- sqlQuery(conn_db, "select * from sus_scrofa_sscrofa9_56_gtf where attributes like '%ENSSSCT00000000756%' order by f_start")
lps<- sqlQuery(conn_db, "select * from cufflinks_transcript_gtf where transcript_id like 'LPS.152114.0' order by start")
ctrl<-sqlQuery(conn_db, "select * from cufflinks_transcript_gtf where transcript_id like 'CTRL.160181.0' order by start")

lim1<- min(c(ref$f_start, lps$start, ctrl$start))
lim2<- max(c(ref$f_end, lps$end, ctrl$end))

main= paste('GAPDH (', round(mean(c(lps$fpkm, ctrl$fpkm)),0), ' FPKM)', sep= '')
plot(x= 1, y=1, xlim= c(lim1, lim2), type= 'n', xlab= '', ylab= '', axes= F, ylim=c(0,3), main=main)
axis(side=1)
ref_ex<- ref[ref$feature == 'exon',]
segments(x0= ref_ex$f_start, x1= ref_ex$f_end, y0= 1, y1= 1, lwd= 5, col= 'red')
lps_ex<- lps[lps$feature == 'exon',]
segments(x0= lps_ex$start, x1= lps_ex$end, y0= 2, y1= 2, lwd= 5, col= 'lightblue')
ctrl_ex<- ctrl[ctrl$feature == 'exon',]
segments(x0= ctrl_ex$start, x1= ctrl_ex$end, y0= 3, y1= 3, lwd= 5, col= 'lightblue')
ctrl2_ex<- ctrl2[ctrl2$feature == 'exon',]
segments(x0= ctrl2_ex$start, x1= ctrl2_ex$end, y0= 3, y1= 3, lwd= 5, col= 'darkblue')

mtext(side= 2, at= 1:3, text= c('ensembl', 'LPS', 'CTRL'), adj= 0.5, las= 1, cex= 0.75)

# --------------------------[ TNFA_PIG ]----------------------------
ref<- sqlQuery(conn_db, "select * from sus_scrofa_sscrofa9_56_gtf where attributes like '%ENSSSCT00000001533%' order by f_start")
lps<- sqlQuery(conn_db, "select * from cufflinks_transcript_gtf where transcript_id like 'LPS.175619.0' order by start")
ctrl<-sqlQuery(conn_db, "select * from cufflinks_transcript_gtf where transcript_id like 'CTRL.184978.0' order by start")

lim1<- min(c(ref$f_start, lps$start, ctrl$start))
lim2<- max(c(ref$f_end, lps$end, ctrl$end))

main= paste('TNFA (', round(mean(c(lps$fpkm, ctrl$fpkm)),0), ' FPKM)', sep= '')
plot(x= 1, y=1, xlim= c(lim1, lim2), type= 'n', xlab= "Chromosome position", ylab= '', axes= F, ylim=c(0,3), main=main)
axis(side=1)
ref_ex<- ref[ref$feature == 'exon',]
segments(x0= ref_ex$f_start, x1= ref_ex$f_end, y0= 1, y1= 1, lwd= 5, col= 'red')
lps_ex<- lps[lps$feature == 'exon',]
segments(x0= lps_ex$start, x1= lps_ex$end, y0= 2, y1= 2, lwd= 5, col= 'lightblue')
ctrl_ex<- ctrl[ctrl$feature == 'exon',]
segments(x0= ctrl_ex$start, x1= ctrl_ex$end, y0= 3, y1= 3, lwd= 5, col= 'lightblue')
ctrl2_ex<- ctrl2[ctrl2$feature == 'exon',]
segments(x0= ctrl2_ex$start, x1= ctrl2_ex$end, y0= 3, y1= 3, lwd= 5, col= 'darkblue')

mtext(side= 2, at= 1:3, text= c('ensembl', 'LPS', 'CTRL'), adj= 0.5, las= 1, cex= 0.75)

# ------------------------------------------------------------------------------------------------
dim(gtf)

(34934 - 61323)/61323 

(8105935 - 7801707) / 8105935