#!/usr/bin/env Rscript

args<- commandArgs(trailingOnly = TRUE)

if(length(args) == 0){
	quit(save= 'no', status= 100)
}


dat<- read.table(args[1], header= TRUE, sep= "\t")
dat$tot<- rowSums(dat[, c('cnt_met', 'cnt_unmet', 'cnt_mismatch', 'cnt_noncytosine')])
dat$pct_met<- 100 * (dat$cnt_met / (dat$cnt_met + dat$cnt_unmet))
dat$pct_met<- ifelse(is.nan(dat$pct_met), 0, dat$pct_met)

dat$pct_mism<- 100 * (dat$cnt_mism / dat$tot)
dat$pct_mism<- ifelse(is.nan(dat$pct_mism), 0, dat$pct_mism)

dat$pct_noncytosine<- 100 * (dat$cnt_noncytosine / dat$tot)
dat$pct_noncytosine<- ifelse(is.nan(dat$pct_noncytosine), 0, dat$pct_noncytosine)

## Overwrite input file with now additional columns
write.table(x= dat, file= args[1], col.names= TRUE, sep= "\t", quote= FALSE, row.names= FALSE)

nreads<- length(unique(dat$read))

outpdf<- paste(sub("\\.txt$", "", args[1], perl= TRUE), ".pdf", sep= "")
pdf(outpdf, height= ifelse(nreads == 1, 9.5/2.54, 18/2.54), width= 16/2.54, 
	pointsize= ifelse(nreads == 1, 9, 10))
par(mfrow= c(nreads * 2, 1), las= 1, mgp= c(2, 0.5, 0), mar= c(2, 3, 1, 3), 
	oma= c(1, 0, 0, 0), bty= 'u', tcl= -0.25)
for(x in unique(dat$read) ){
	pdat<- dat[dat$read == x, ]
	
	## Pct methylated
	plot(x= pdat$position, y= pdat$pct_met, ylab= '% Unconverted', xlab= '', type= 'n', 
		ylim= c(0, max(dat$pct_met)))
	rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col= "gray95", border= NA)
	points(x= pdat$position, y= pdat$pct_met, pch= 19, cex= 0.5, col= 'grey30', type= 'o')
	
	## Tot reads
	par(new= TRUE)
	plot(x= pdat$position, y= pdat$tot, type= 'l', yaxt= 'n', xlab= '', ylab= '', col= 'blue', 
		ylim= c(0, max(pdat$tot)))
	axis(side= 4, col= 'blue', col.axis= 'blue', las= 0)
	mtext(side= 4, text= 'N reads', line= 2, col= 'blue', las= 0)
	mtext(paste("Read", x), side= 3, line= -0.1)
	grid(col= 'grey30')
	
	## Percent mismatch
	plot(x= pdat$position, y= pdat$pct_mism, ylab= '% Mismatches (A, G, N)', xlab= '', type= 'n', 
		ylim= c(0, max(dat$pct_mism))) #4682B4
	rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col= "grey85", border= NA) ## #F5DEB3  #FFE4B580
	points(x= pdat$position, y= pdat$pct_mism, pch= 19, cex= 0.5, col= 'grey30', type= 'o')
	grid(col= 'grey30')
	
	## Percent non cytosine
	par(new= TRUE)
	plot(x= pdat$position, y= pdat$pct_noncytosine, type= 'l', yaxt= 'n', xlab= '', ylab= '', col= 'blue', 
		ylim= c(0, max(pdat$pct_noncytosine)))
	axis(side= 4, col= 'blue', col.axis= 'blue', las= 0)
	mtext(side= 4, text= '% Non cytosine', line= 2, col= 'blue', las= 0)
}
mtext(text= "Position", side= 1, line= -0.2, outer= TRUE)
dev.off()

