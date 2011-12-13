#
#  Explore CAGE expression
#  Library CAGE050810, alignment and processing 16/12/2010
#  LabBook ref 16/12/2010


library(RODBC)

conn<- odbcConnect(dsn= 'pgCAGE')

#----------------------[ I/O ]------------------
# Coordinates to extract
## IDO; ENSSSCT00000007675; Chromosome 17: 9,666,698-9,681,681
## TNFA; ENSSSCT00000001533; 7: 27658098-27661703
## IL1; ENSSSCT00000008863; 3: 39059355-39072049; reverse ## IL1
## ACTB; ENSSSCG00000007585; 3:3,868,537-3,871,915 forward
rname<- 3
left_pos<-  3868537 - 1000
right_pos<- 3871915 + 1000


#pile_f<- read.table('F:/data/20101217_CAGE050810/cage_050810_bwa_20101216.clean.single.forward.slice.pileup', 
#    header=T, sep='\t')

#pile_r<- read.table('F:/data/20101217_CAGE050810/cage_050810_bwa_20101216.clean.single.reverse.slice.pileup', 
#    header=T, sep='\t')

#rnaseq<- read.table('F:/data/20100630_RNAseq_pileup/20100630_RNAseq_CTRL_contig.pileup/20100630_RNAseq_CTRL_3.pileup',
#    header=T, sep='\t')
#rnaseq_tnfa<- rnaseq[rnaseq$pos>= 27658098 & rnaseq$pos<= 27661703 & rnaseq$rname == 7,]
#rnaseq_actb<- rnaseq[rnaseq$pos>= (3868537 - 1000) & rnaseq$pos<= (3871915 + 1000) & rnaseq$rname == 3,]

rm(rnaseq)

ensembl_gtf<- sqlQuery(conn, paste("
    SELECT *, gtf_attribute(attributes, 'transcript_id') AS transcript_id FROM sus_scrofa_sscrofa9_59_gtf 
    WHERE rname = '", rname, "' AND f_start >= ", left_pos,  " AND f_end <= ", right_pos, " AND feature LIKE 'exon'
    ORDER BY rname, f_start;
    ", sep= ''))

slice<- pile_f[pile_f$pos>= left_pos & pile_f$pos<= right_pos & pile_f$rname == rname,]
slice<- pile_r[pile_r$pos>= left_pos & pile_r$pos<= right_pos & pile_r$rname == rname,]


## ACTB
rname<- 3
left_pos<-  3868537 - 500
right_pos<- 3871915 + 500

#rnaseq_actb<- rnaseq[rnaseq$pos>= (3868537 - 1000) & rnaseq$pos<= (3871915 + 1000) & rnaseq$rname == 3,]


slice<- pile_f[pile_f$pos>= left_pos & pile_f$pos<= right_pos & pile_f$rname == rname,]
rnaseq_actb<- rnaseq[rnaseq$pos>= (3868537 - 1000) & rnaseq$pos<= (3871915 + 1000) & rnaseq$rname == 3,]

par(cex=0.8)
with(rnaseq_actb, {
    plot(x= pos, y= no_reads, type='h', col= 'dodgerblue', ylab= 'Read count', xlab= 'Position', 
         ylim= c(0, 12000), xaxt='n', xlim= c(3868537 - 250, 3871915 + 250))
    })
with(slice, {
    points(x= pos, y= read_count, type='h', col= 'firebrick')
    })
axis(side= 1, at= seq(left_pos, right_pos, length.out= 10), 
    labels= format(seq(left_pos, right_pos, length.out= 10), big.mark=','))

rect(xleft= ensembl_gtf$f_start, xright= ensembl_gtf$f_end, ybottom= 11500, ytop= 12000, col= 'olivedrab4', border= NA)
segments(x0= min(ensembl_gtf$f_start), x1= max(ensembl_gtf$f_end), y0= 11750, col= 'olivedrab4', lwd= 1.5)

## CSF1R; Chromosome 2: 136,526,706-136,552,808 reverse
rname<- 2
left_pos<-  136526706 - 500
right_pos<- 136552808 + 500
rnaseq<- read.table('F:/data/20100630_RNAseq_pileup/20100630_RNAseq_CTRL_contig.pileup/20100630_RNAseq_CTRL_2.pileup',
    header=T, sep='\t')

slice<- pile_r[pile_r$pos>= left_pos & pile_r$pos<= right_pos & pile_r$rname == rname,]
rnaseq_csf1r<- rnaseq[rnaseq$pos>= left_pos & rnaseq$pos<= right_pos & rnaseq$rname == rname,]

par(cex=0.8)
with(rnaseq_csf1r, {
    plot(x= pos, y= no_reads, type='h', col= 'dodgerblue', ylab= 'Read count', xlab= 'Position', 
         ylim= c(0, 2500), xaxt='n', xlim=c(136552808-200, 136552808+200))
    })
with(slice, {
    points(x= pos, y= read_count, type='h', col= 'firebrick')
    })
# axis(side= 1, at= seq(left_pos, right_pos, length.out= 10), 
#     labels= format(seq(left_pos, right_pos, length.out= 10), big.mark=','))
axis(side= 1, at= seq(136552808-200, 136552808+200, length.out= 10), 
    labels= format(seq(136552808-200, 136552808+200, length.out= 10), big.mark=','))

rect(xleft= ensembl_gtf$f_start, xright= ensembl_gtf$f_end, ybottom= 2350, ytop= 2500, col= 'olivedrab4', border= NA)
segments(x0= min(ensembl_gtf$f_start), x1= max(ensembl_gtf$f_end), y0= 2425, col= 'olivedrab4', lwd= 1.5)

