#
#  Explore expression of IDO1 in LPS and CTRL libraries
#


library(RODBC)

conn<- odbcConnect(dsn= 'pgCAGE')

cuff<- sqlQuery(conn, 
  "SELECT * FROM cufflinks_transcript_gtf 
       WHERE feature LIKE 'transcript' AND
       source IN('20100317_RNAseq_LPS', '20100317_RNAseq_CTRL')"
  )
dim(cuff)

hist(log(cuff$fpkm), xlab= 'log(FPKM)', 
    main= 'Distibution of gene expression in CTRL and LPS', breaks= 30)
abline(v= log(cuff$fpkm[cuff$transcript_id == 'ENSSSCT00000007675']), 
    col= c('dodgerblue', 'firebrick'), lty= 'dotted', lwd= 2)
legend('topleft', legend= c('IDO1 in CTRL', 'IDO1 in 7h LPS'), 
    col= c('dodgerblue', 'firebrick'), lty= 'dashed', lwd= 2)
savePlot('F:/data/20101104_IDO1/rnaseq_expr.emf', 'emf')

pile_ctrl<- read.table('F:/data/20100630_RNAseq_pileup/20100630_RNAseq_CTRL_contig.pileup/20100630_RNAseq_CTRL_17.pileup', 
    header=T, sep='\t')
pile_lps<- read.table('F:/data/20100630_RNAseq_pileup/20100630_RNAseq_LPS_contig.pileup/20100630_RNAseq_LPS_17.pileup', 
    header=T, sep='\t')
pile_lps[1:10,]

left_pos<- 9666000 ## 9666698
right_pos<- 9682000 ## 9681681

ensembl_gtf<- sqlQuery(conn, paste("
select *, gtf_attribute(attributes, 'transcript_id') AS transcript_id from sus_scrofa_sscrofa9_59_gtf 
    where rname = '17' and 
          f_start >= ", left_pos,  " and 
          f_end <= ", right_pos, " and
          feature like 'exon'
    order by rname, f_start;
  "))

par(mar=c(5,4,6,2), cex= 0.8)
plot(x= pile_lps[pile_lps$pos>= left_pos & pile_lps$pos<= right_pos, 'pos'],
     y= pile_lps[pile_lps$pos>= left_pos & pile_lps$pos<= right_pos, 'no_reads'],
     type= 'h', col= 'firebrick', 
     xlab= 'chr 17: 9,666,000 - 9,682,000', xaxt='n',
     ylab= 'No. reads', main= 'Raw expression of IDO1 in CTRL and LPS alveolar macrophages')
axis(side= 1, at= seq(left_pos, right_pos, by= 2000), 
    labels= format(seq(left_pos, right_pos, by= 2000), big.mark=','))

points(x= pile_ctrl[pile_ctrl$pos>= left_pos & pile_ctrl$pos<= right_pos, 'pos'],
     y= pile_ctrl[pile_ctrl$pos>= left_pos & pile_ctrl$pos<= right_pos, 'no_reads'],
     type= 'h', col= 'dodgerblue')

legend(x= left_pos, y= 3600, , legend= c('CTRL', 'LPS'), xpd= T, horiz= T,
    col= c('dodgerblue', 'firebrick'), lty= 'solid', lwd= 2, bty= 'n')

rect(xleft= ensembl_gtf$f_start, xright= ensembl_gtf$f_end, ybottom= 3000, ytop= 3200, col= 'olivedrab4', border= NA)
segments(x0= min(ensembl_gtf$f_start), x1= max(ensembl_gtf$f_end), y0= 3100, col= 'olivedrab4', lwd= 1.5)

savePlot('F:/data/20101104_IDO1/ido1_expr.jpg', 'jpg')
write.table(x= pile_ctrl[pile_ctrl$pos>= left_pos & pile_ctrl$pos<= right_pos, ],
            file= 'F:/data/20101104_IDO1/ido1_ctrl.csv', sep=',', row.names=F)
write.table(x= pile_lps[pile_lps$pos>= left_pos & pile_lps$pos<= right_pos, ],
            file= 'F:/data/20101104_IDO1/ido1_lps.csv', sep=',', row.names=F)

