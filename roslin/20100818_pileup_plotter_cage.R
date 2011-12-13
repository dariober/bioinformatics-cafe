## For these files see eddie on /exports/work/vet_roslin_nextgen/dario/miscellanea/20100701
pileup_chr4_ctrl<- read.table('D:/Tritume/20100630_RNAseq_CTRL_chr4.pileup', header= T, sep= '\t')
pileup_chr4_lps<- read.table('D:/Tritume/20100630_RNAseq_LPS_chr4.pileup', header= T, sep= '\t')

names(pileup_chr4_ctrl)
pileup_chr4_ctrl<- pileup_chr4_ctrl[pileup_chr4_ctrl$pos >= 94700000 & pileup_chr4_ctrl$pos <= 95050000, ]
pileup_chr4_lps<- pileup_chr4_lps[pileup_chr4_lps$pos >= 94700000 & pileup_chr4_lps$pos <= 95050000, ]

dim(pileup_chr4_ctrl); dim(pileup_chr4_lps)
write.table(pileup_chr4_ctrl, file= 'D:/Tritume/chr4_ctrl_slice.pileup', sep= '\t', row.names=F)
write.table(pileup_chr4_lps, file= 'D:/Tritume/chr4_lps_slice.pileup', sep= '\t', row.names=F)

library(RODBC)
conn_db<- odbcConnect(dsn= 'pgCAGE', case= 'nochange')

ensembl_gtf<- sqlQuery(conn_db, "
select *, gtf_attribute(attributes, 'transcript_id') AS transcript_id from sus_scrofa_sscrofa9_56_gtf 
    where rname = '4' and 
          f_start >= 94700000 and 
          f_end <= 95050000 and
          feature like 'exon'
    order by rname, f_start;
  ")

## Import CAGE tags (see LabBook 16/08/2010)
cage_tags<- sqlQuery(conn_db, "
    SELECT * FROM sam_bowtie_human_vs_pig 
    WHERE rname like '4' AND 
          pos BETWEEN 94928458-1000 AND 94956044+1000 
    ORDER BY pos")

## Region to zoom in. Plot large window first than use locator()
## x<- locator() 
## Or use the following coordinates:
x<- list(); x$x<- c(94924826, 94961651)

par(mfrow= c(2,1), mar=c(3, 4, 1, 1), oma= c(1, 0, 0, 0))
  # ----------------------[ Large window ]------------------------------
  largew<- c(94700000, 95050000)
  # Plot region
  plot(x= pileup_chr4_lps$pos, y= pileup_chr4_lps$no_reads, type= 'n',
    xlab= '', ylab= 'Coverage', main= '', xlim= largew, ylim= c(0, 80), xaxt= 'n');
  mtext(text=paste('Chr 4:', largew[1], '-', largew[2]), side= 3, line= -2, adj= 0.05, font= 2)

  # Annotate x axis
  px<- axis(1, labels=F, tick= F)
  axis(1, at= px, labels= formatC(x= px/1000, big.mark = ",", format = "d"), las= 1)

  # Highlight zoom'd region.
  rect(xleft= x$x[1], xright= x$x[2], ybottom= -10, ytop= 100, col= 'azure2', border = NA )
  box() ## Redraw box

  # Draw actual data: LPS
  points(x= pileup_chr4_lps$pos, y= pileup_chr4_lps$no_reads, type= 'h', lwd= 2, col= 'red')

  # Draw actual data: CTRL. Note that the bars are shifted right using +1500
  points(x= pileup_chr4_ctrl$pos + 1500, y= pileup_chr4_ctrl$no_reads, type= 'h', lwd= 2, col= 'blue')
 
  # Show where transcripts are. Different colour = different transcript
  segments(x0= ensembl_gtf$f_start, x1= ensembl_gtf$f_end, y0= 20, col= as.numeric(ensembl_gtf$transcript_id), lwd= 5)


  #-----------------------[ Zoom into region of interest ]-------------------------------
  # LPS
  plot(x= pileup_chr4_lps$pos, y= pileup_chr4_lps$no_reads, type= 'h', ylab= 'Coverage', main= '', lwd= 2, col= 'red',
    xlim= x$x, ylim= c(0, 20), xaxt='n');
  mtext(side= 1, text= 'Genome coordinate (kb)', line= 2.5)
  mtext(text= 'ENSSSCT00000007031', font= 2, side= 3, adj= 0.05, line= 0.1)
  p<- axis(1, at= seq(round(x$x[1]/1000)*1000, round(x$x[2]/1000)*1000, length.out= 11), labels= FALSE, tick= F)
  axis(1, at= p, labels= formatC(x= p/1000, big.mark = ",", format = "d"), las= 1) 

  abline(v= c(ensembl_gtf$f_start, ensembl_gtf$f_end), col= 'blue', lty= 'dotted')
  segments(x0= ensembl_gtf$f_start, x1= ensembl_gtf$f_end, y0= 20, col= 'green', lwd= 3)
  
  # CTRL
  points(x= pileup_chr4_ctrl$pos, y= pileup_chr4_ctrl$no_reads, type= 'h', col= 'blue', lwd= 2); 

  # Human CAGE tags
  with(cage_tags[cage_tags$flag == 0, ],
         arrows(x0= pos, y0= rep(19, length(pos)),
           x1= pos + 20, y1= rep(19, length(pos)), 
           code= 2, length= 0.085)
       )
   with(cage_tags[cage_tags$flag == 16, ],
         arrows(x0= pos, y0= rep(18, length(pos)),
           x1= pos + 20, y1= rep(18, length(pos)), 
           code= 1, length= 0.085)
       )
 
  # Legend
   legend('topright', legend= c("LPS", "CTRL", 'Exons'), col= c("red", "blue", "green"), lwd= 2, cex= 0.85, bty='n', xpd= TRUE)

cage_tags_chr4<- sqlQuery(conn_db, "
    SELECT * FROM sam_bowtie_human_vs_pig 
    WHERE rname like '4'
    ORDER BY pos")
dim(cage_tags_chr4)
h1<- hist(cage_tags_chr4$pos, breaks= seq(0, max(cage_tags_chr4$pos)+30000, by= 30000), plot=T)
mean(h1$counts)