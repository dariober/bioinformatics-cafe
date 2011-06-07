#
# Plot a slice of a pileup file
# Slice produced typically by (20100629)_pileup_slice.py
# and should look like:
# -------------------
 rname	pos	no_reads
 4	94788700	0
 4	94788701	0
 4	94788702	0
 4	94788703	0
 4	94788704	0
 4	94788705	0
# -------------------

# ---------------------[ Import data ]---------------------------

## For these files see also eddie on /exports/work/vet_roslin_nextgen/dario/miscellanea/20100701
##
pileup_chr4_ctrl<- read.table('F:/data/20100630_RNAseq_pileup/20100630_RNAseq_CTRL_contig.pileup/20100630_RNAseq_CTRL_4.pileup', 
    header= T, sep= '\t', colClasses= c('character', 'numeric', 'character', 'character', 'numeric', 'numeric', 'numeric', 'numeric'))
## Remove unnecessary columns
pileup_chr4_ctrl<- pileup_chr4_ctrl[, c(1,2,8)]
pileup_chr4_ctrl[1:10,]

pileup_chr4_lps<- read.table('F:/data/20100630_RNAseq_pileup/20100630_RNAseq_LPS_contig.pileup/20100630_RNAseq_LPS_4.pileup',
    header= T, sep= '\t', colClasses= c('character', 'numeric', 'character', 'character', 'numeric', 'numeric', 'numeric', 'numeric'))
## Remove unnecessary columns
pileup_chr4_lps<- pileup_chr4_lps[, c(1,2,8)]
pileup_chr4_lps[1:10,]


# ---------------------------[ Interval to plot ]--------------------
#
# ENSSSCG00000010716 14:138735210-138752998
# ENSSSCG00000010787 14:148349161-148359147
# ENSSSCG00000006350 4:92582670-92601887
# ENSSSCG00000006678 4:103242867-103251849
left_pos<-  103241000
right_pos<- 103253000

pileup_chr4_ctrl<- pileup_chr4_ctrl[pileup_chr4_ctrl$pos >= left_pos & pileup_chr4_ctrl$pos <= right_pos, ]
pileup_chr4_lps<- pileup_chr4_lps[pileup_chr4_lps$pos >= left_pos & pileup_chr4_lps$pos <= right_pos, ]

dim(pileup_chr4_ctrl); dim(pileup_chr4_lps)

# write.table(pileup_chr4_ctrl, file= 'D:/Tritume/chr4_ctrl_slice.pileup', sep= '\t', row.names=F)
# write.table(pileup_chr4_lps, file= 'D:/Tritume/chr4_lps_slice.pileup', sep= '\t', row.names=F)

library(RODBC)
conn_db<- odbcConnect(dsn= 'pgCAGE', case= 'nochange')

ensembl_gtf<- sqlQuery(conn_db, paste("
select *, gtf_attribute(attributes, 'transcript_id') AS transcript_id from sus_scrofa_sscrofa9_56_gtf 
    where rname = '4' and 
          f_start >= ", left_pos,  " and 
          f_end <= ", right_pos, " and
          feature like 'exon'
    order by rname, f_start;
  "))

## Region to zoom in. Plot large window first than use locator()
## x<- locator() 
## Or use the following coordinates:
x<- list(); x$x<- c(94924826, 94961651)

par(mfrow= c(2,1), mar=c(3, 4, 1, 1), oma= c(1, 0, 0, 0))
  
# ----------------------[ Large window ]------------------------------

  largew<- c(left_pos, right_pos)
  # Plot region
  plot(x= pileup_chr4_lps$pos, y= pileup_chr4_lps$no_reads, type= 'n',
    xlab= '', ylab= 'Coverage', main= '', xlim= largew, ylim= c(0, 60), xaxt='n');
  mtext(text=paste('Chr 4:', largew[1], '-', largew[2]), side= 3, line= -3, adj= 0.05, font= 2)

  # Annotate x axis
  px<- axis(1, labels=F, tick= F)
  axis(1, at= px, labels= formatC(x= px/1000, big.mark = ",", format = "d"), las= 1)

  # Draw actual data: LPS
  points(x= pileup_chr4_lps$pos, y= pileup_chr4_lps$no_reads, type= 'h', lwd= 1, col= 'red')

  # Draw actual data: CTRL. Note that the bars are shifted right
  points(x= pileup_chr4_ctrl$pos+20, y= pileup_chr4_ctrl$no_reads, type= 'h', lwd= 1, col= 'blue')
 
  # Show where transcripts are. Different colour = different transcript
  rect(xleft= ensembl_gtf$f_start, xright= ensembl_gtf$f_end, ybottom= 56, ytop= 60, col= 'firebrick4', border= NA)
  segments(x0= min(ensembl_gtf$f_start), x1= max(ensembl_gtf$f_end), y0= 58, col= 'firebrick4', lwd= 1.5)

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
  
  # Legend
  legend('topright', legend= c("LPS", "CTRL", 'Exons'), col= c("red", "blue", "green"), lwd= 2, cex= 0.85, bty='n', xpd= TRUE)

savePlot('M:/Documents/LabBook/LabBook_Figures/20100705_ENSSSCT00000007031.emf', 'emf')
savePlot('M:/Documents/LabBook/LabBook_Figures/20100705_ENSSSCT00000007031.jpeg', 'jpeg')

# -----------------------------------[ End of graph ]--------------------------------------


#------------------------------------[ Tritume ]------------------------------------

plot(x= pileup_chr12_ctrl$pos, y= pileup_chr12_ctrl$no_reads, type='h', 
    xlab= 'Genome coordinate', ylab= 'Coverage', main= 'Coverage', lwd= 2, col= 'lightgrey',
    xlim= c(94928458, 94956044), ylim= c(0,15)); 



is.numeric(pileup_chr4_ctrl$cov)

library(RODBC)
conn_db<- odbcConnect(dsn= 'pgCAGE', case= 'nochange')
cufflinks_gtf<- sqlQuery(conn_db, "
    SELECT * FROM cufflinks_transcript_gtf 
        WHERE transcript_id LIKE 'ENSSSCT00000007031' AND
              source LIKE '20100317_RNAseq_LPS' AND
              feature LIKE 'exon'
    ")

ensembl_gtf<- sqlQuery(conn_db, "
select * from sus_scrofa_sscrofa9_56_gtf 
    where rname = '4' and 
          f_start >= 94700000 and 
          f_end <= 95050000 and
          feature like 'exon'
    order by rname, f_start;
  ")


slice<- read.table('D:/Tritume/20100409_RNAseq_LPS_sscrofa9.56.pileup.slice', header=T, sep='\t')
names(slice)


# Entire window 
# 
plot(x= slice$pos, y= slice$cov, type='l', xlab= 'Genome coordinate', ylab= 'Coverage', main= 'Coverage', lwd= 2); 
abline(v= c(ensembl_gtf$f_start, ensembl_gtf$f_end), col= 'blue', lty= 'dotted')
segments(x0= ensembl_gtf$f_start, x1= ensembl_gtf$f_end, y0= 15, col= 'red', lwd= 3)


# Zoom in using locator()
#
x<- locator()$x;
plot(x= slice$pos, y= slice$cov, type='l', xlim= x, xlab= 'Genome coordinate', ylab= 'Coverage', main= paste('Coverage between\n', round(x[1]), 'and', round(x[2])))
abline(v= c(ensembl_gtf$f_start, ensembl_gtf$f_end), col= 'blue', lty= 'dotted')

x<- c(94927000, 94957000)

savePlot('D:/Tritume/slice.emf', 'emf')

seq(round(x$x[1]/1000)*1000, round(x$x[2]/1000)*1000, by= 5000)







