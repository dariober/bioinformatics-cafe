#!/usr/bin/env Rscript

args<- commandArgs(trailingOnly = TRUE)

if(length(args) == 0){
	quit(save= 'no', status= 100)
}

dat<- read.table(args[1], header= TRUE, sep= "\t")
dat$tot<- rowSums(dat)
dat$pct_met<- 100 * (dat$cnt_met / (dat$cnt_met + dat$cnt_unmet))
dat$pct_met<- ifelse(is.nan(dat$pct_met), 0, dat$pct_met)

dat$pct_mism<- 100 * (dat$cnt_mism / dat$tot)
dat$pct_mism<- ifelse(is.nan(dat$pct_mism), 0, dat$pct_mism)

write.table(x= dat, file= args[1], col.names= TRUE, sep= "\t", quote= FALSE, row.names= FALSE)

nreads<- length(unique(dat$read)) 

outpdf<- paste(sub("\\.txt$", "", args[1], perl= TRUE), ".pdf", sep= "")
pdf(outpdf, height= ifelse(nreads == 1, 7.5/2.54, 14/2.54), width= 16/2.54, pointsize= 10)
par(mfrow= c(nreads, 1), las= 1, mgp= c(2, 0.7, 0), mar= c(2, 3, 1, 0.5), 
	oma= c(1, 0, 0, 0), bty= 'l')
for(x in unique(dat$read) ){
	pdat<- dat[dat$read == x, ]
	plot(x= pdat$position, y= pdat$pct_met, ylab= '% Unconverted', xlab= '', type= 'n', 
		ylim= c(0, max(dat$pct_met)))
	rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col= "gray95", border= NA)
	points(x= pdat$position, y= pdat$pct_met, pch= 19, cex= 0.5, col= 'grey30', type= 'o')
	mtext(paste("Read", x), side= 3, line= -0.3)
	grid(col= 'grey30')
}
mtext(text= "Position", side= 1, line= -0.2, outer= TRUE)
dev.off()

