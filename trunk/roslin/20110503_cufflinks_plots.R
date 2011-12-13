# ------------------------------------------------------------------------------
# Number of transcripts identified by RNAseq by FPKM
# ------------------------------------------------------------------------------

library(RODBC)
library(sqldf)
conn<- odbcConnect(dsn= 'pgVitelleschi')

cuff_gtf<- sqlQuery(conn, " select * from cufflinks_transcript_gtf where feature = 'transcript' and fpkm > 0 order by fpkm")
cuff_gtf[1:10,]

cuff_am_ctrl<- sqldf("select * from cuff_gtf where source = '20110202_am_ctrl' ")
cuff_am_ctrl$ord<- ntrans<- 1:nrow(cuff_am_ctrl)
cuff_am_lps<- sqldf("select * from cuff_gtf where source = '20110202_am_lps' ")
cuff_am_lps$ord<- ntrans<- 1:nrow(cuff_am_lps)
cuff_am_ctrl[1:10,]

transcr<-  data.frame(ensembl_transcript_id= c('ENSSSCT00000008324', 'ENSSSCT00000001533', 'ENSSSCT00000008863'), gene_name= c('ACTB', 'TNFA', 'IL1A'))
ntrans_ctrl<- sqldf("
    SELECT cuff_am_ctrl.*, gene_name
    FROM cuff_am_ctrl, transcr
    WHERE transcript_id IN (transcr.ensembl_transcript_id)
")
ntrans_lps<- sqldf("
    SELECT cuff_am_lps.*, gene_name
    FROM cuff_am_lps, transcr
    WHERE transcript_id IN (transcr.ensembl_transcript_id)
")


graphics.off()
windows(height= 13/2.54, width= 12/2.54)
par(mgp= c(2, 0.75, 0), cex.lab= 0.85, cex.axis= 0.85, las= 1, mfrow= c(2,1), mar= c(0.5, 3, 0.5, 0.1), oma= c(3,0,1,0), xpd= FALSE)
    plot(cuff_am_ctrl$fpkm, cuff_am_ctrl$ord/1000, type= 'l', log= 'x', xaxt= 'n',
         xlab= '', ylab= 'No. transcripts x1000', lwd= 2, col= 'grey25', xlim= c(0.1, max(cuff_gtf$fpkm)))
    points(x= ntrans_ctrl$fpkm, y= ntrans_ctrl$ord/1000, pch= 19, col= c('red', 'blue', 'grey25'))
    segments(x0= ntrans_ctrl$fpkm, x1= ntrans_ctrl$fpkm,
             y0= 0, y1= ntrans_ctrl$ord/1000, col= c('red', 'blue', 'grey25'), lty= 'dotted', xpd= FALSE)
    text(x= ntrans_ctrl$fpkm, y= ntrans_ctrl$ord/1000, labels= ntrans_ctrl$gene_name, adj= c(0,1), cex= 0.8)
    grid()
    legend('topleft', legend= 'AM CTRL', bg= 'white')
    plot(cuff_am_lps$fpkm, cuff_am_lps$ord/1000, type= 'l', log= 'x', xaxt= 'n',
         xlab= 'FPKM', ylab= 'No. transcripts x1000', lwd= 2, col= 'grey25', xlim= c(0.1, max(cuff_gtf$fpkm)))
    points(x= ntrans_lps$fpkm, y= ntrans_lps$ord/1000, pch= 19, col= c('red', 'blue', 'grey25'))
    segments(x0= ntrans_lps$fpkm, x1= ntrans_lps$fpkm,
             y0= 0, y1= ntrans_lps$ord/1000, col= c('red', 'blue', 'grey25'), lty= 'dotted', xpd= FALSE)
    axis(side= 1, at= c(0.1, 10^(seq(0, 5) )), labels= formatC(c(0.1, 10^(seq(0, 5) )), format= 'd', big.mark= ','))
    grid()
    legend('topleft', legend= 'AM LPS', bg= 'white')
mtext(side= 1, xpd= TRUE, text= 'Expression level (FPKM)', cex= 0.85, line= 2)
title(main= 'No. of transcripts and expression', outer= TRUE, cex.main= 1)
savePlot('M:/Documents/LabBook/LabBook_Figures/20110503_cufflinks_plot_am.emf', 'emf')

# -----------------------------------------------------------------------------
#  TRITUME
# -----------------------------------------------------------------------------

#plot(cuff_am_ctrl$fpkm, ntrans, type= 'l', log= 'x', xaxt= 'n')
#axis(side= 1, at= 10^(-5:5), labels= formatC(10^(-5:5), format= 'd', big.mark = ","))

points(cuff_am_ctrl$fpkm[cuff_am_ctrl$transcript_id %in% transcr], ntrans)
hist(cuff_am_ctrl$fpkm[cuff_am_ctrl$fpkm < 100])

length(cuff_am_ctrl$fpkm[cuff_am_ctrl$fpkm < 10])

options()
tfpkm<- sqldf("select * from cuff_am_lps, transcr where transcript_id in (transcr.ensembl_transcript_id) ")
ntrans<- 1:nrow(cuff_am_ctrl)
