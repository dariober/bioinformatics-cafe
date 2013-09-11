"
DESCRIPTION:
    Read matrix files for chrom and produce GLM output 

USAGE:
    -- Processing 'chrM'
    R CMD BATCH -chrM bsdata_glm.R chrM/bsdata_glm.Rout

SEE ALSO:
    For processing command line args from R BATCH

cd /lustre/sblab/berald01/projects/20120522-oxbsseq-mesc/slx-6278/genome/bsdata
"

# ------------------------- Setting up ----------------------------------------
source('/home/berald01/svn_checkout/bioinformatics-misc/glm.R')
source('/home/berald01/svn_checkout/bioinformatics-misc/BSdata.R')
source('/home/berald01/svn_checkout/bioinformatics-misc/makeTransparent.R')

NPROCS<- 4

args <- commandArgs(trailingOnly = F)
CHROM<- args[length(args)]
CHROM<- sub("^-", "", CHROM, perl= TRUE)

setwd(CHROM)
cat(getwd())
# -----------------------------------------------------------------------------

cpg<- read.bsdata('genome', mat= c('loci', 'cnt_met', 'tot_reads'), gzip= FALSE, nrows= -1)
cpg@design$bs<- c('oxBS', 'oxBS', 'oxBS', 'BS', 'BS', 'BS')
cpg@pct_met<- as.ffdf(100 * (cpg@cnt_met[,] / cpg@tot_reads[,]))

# save.BSdata(cpg, 'cpg')
# ffload('cpg')


# Read count
# ----------
pdf('tot_reads.hist.pdf', width= 15/2.54, height= 12/2.54, pointsize= 9)
par(mfrow= c(2, 3), mar= c(2.5, 3, 2, 0.2), oma= c(3,3,2,0))
for(x in colnames(cpg@tot_reads)){
    readQ<- quantile(cpg@tot_reads[, x], 0.95)
    y<- ifelse(cpg@tot_reads[, x] > readQ, readQ, cpg@tot_reads[, x] )
    hist(y, main= '', col= 'grey80', xlab= '', ylab= '', breaks= 10)
    mtext(side= 3, line= -1, text= x, cex= 0.8)
    print(x)
}
mtext(side= 1, text= 'Read depth', font= 2, outer= TRUE)
mtext(side= 2, text= 'Frequency', outer= TRUE, font= 2)
mtext(side= 3, text= sprintf('%s Depth of coverage at CpG sites', CHROM), outer= TRUE, font= 2, line= 0)
dev.off()
system(sprintf('scp tot_reads.hist.pdf $mac_office:$cri_public_projects/20120522-oxbsseq-mesc/Documents/slx-6278/genome/%s', CHROM))

# Methylation 
# -----------
## Select CpG with coverage n<x<m
m<- 5
Q<- 0.95
pdf('pct_met.hist.pdf', width= 15/2.54, height= 12/2.54, pointsize= 9)
par(mfrow= c(2, 3), mar= c(2.5, 3, 2, 0.2), oma= c(3,3,2,0))
    for(x in colnames(cpg@pct_met)){
        y<- cpg@pct_met[, x]
        tot<- cpg@tot_reads[, x]
        readQ<- quantile(cpg@tot_reads[, x], Q)
        y<- y[which(tot >= m & tot <= readQ)]  
        hist(y, main= '', col= 'grey80', xlab= '', ylab= '', breaks= 25)
        mtext(side= 3, line= -1, text= paste(x, '\nn= ', length(y), sep= ''), cex= 0.8, col= 'darkblue')
        print(x)
    }
mtext(side= 1, text= '% methylated', font= 2, outer= TRUE)
mtext(side= 2, text= 'Frequency', outer= TRUE, font= 2)
mtext(side= 3, text= sprintf('%s CpG sites with depth %sx to quantile(%s)', CHROM, m, Q), outer= TRUE, font= 2, line= 0)
dev.off()
system(sprintf('scp pct_met.hist.pdf $mac_office:$cri_public_projects/20120522-oxbsseq-mesc/Documents/slx-6278/genome/%s', CHROM))

# Filter positions with not-too-high coverage
# -------------------------------------------
filter<- which(rowSums(cpg@tot_reads[,] >= 5 & cpg@tot_reads[,] <= 1000) >= 4) ## Positions with coverage x in at least y libs 
cpg.2<- BSdataApply(cpg, FUN= function(x) x[filter,])
pdf('hist_depth_comb_libs.pdf', width= 18/2.54, height= 10/2.54, pointsize= 10)
par(mfrow= c(1,2))
hist(rowSums(cpg.2@tot_reads[,]),
    main= sprintf('Depth of coverage summed across %s libraries\nn= %s', ncol(cpg.2@tot_reads), nrow(cpg.2@tot_reads)), 
    xlab= 'Cumulated depth of coverage', ylab= 'No. CpG', breaks= 20)
hist(rowMeans(cpg.2@pct_met[,]),
    main= sprintf('%s Methylation averaged across %s libraries\nn= %s', CHROM, ncol(cpg.2@tot_reads), nrow(cpg.2@tot_reads)), 
    xlab= '% methylation', ylab= 'No. CpG', breaks= 20)
dev.off()
system(sprintf('scp hist_depth_comb_libs.pdf $mac_office:$cri_public_projects/20120522-oxbsseq-mesc/Documents/slx-6278/genome/%s', CHROM))

# GLM 
# ---
# bsobj<- BSdataApply(cpg.2, FUN= function(x) x[1:10000,])

## UNCOMMENT TO REPEAT GLM
# glmbs<- glmBS.2(bsobj= cpg.2, bs= cpg.2@design$bs, contrast= c('BS', 'oxBS'), family= 'binomial', nthreads= NPROCS)
# row.names(glmbs)<- row.names(cpg.2@cnt_met)
# write.table(glmbs, file= 'glmbs.txt', row.names= FALSE, col.names= TRUE, sep= '\t', quote= FALSE)

## UNCOMMENT TO GET PREVIOUS GLM
glmbs<- read.table('glmbs.txt', header= TRUE, sep= '\t', stringsAsFactors= FALSE)

## MEMO: Estimate is "oxBS - BS" (reverse alphanumeric order)
pneg<- glmbs$pvalue[which(glmbs$estimate > 0)] ## Pvalues of negative estimates.
ppos<- glmbs$pvalue[which(glmbs$estimate < 0)]

pdf('hist_pvals_glmbn.paired.pdf', width= 18/2.54, height= 10/2.54, pointsize= 8)
par(las= 0, mfrow= c(1, 2))
    h1<- hist(pneg, border= 'transparent', add= FALSE, col= 'grey70', freq= FALSE, breaks= 20, xlab= 'P-value', main= "P-values from GLM binom")
    h2<- hist(ppos, border= 'red', freq= FALSE, add= TRUE, breaks= 20)
    points(pch= 17, col= 'blue', x= 0.05, y= 0)
    
    h1<- hist(pneg, border= 'transparent', add= FALSE, col= 'grey70', freq= TRUE, breaks= 20, xlab= 'P-value', main= sprintf("%s P-values from GLM binom\nn= %s", CHROM, nrow(glmbs)))
    h2<- hist(ppos, border= 'red', freq= TRUE, add= TRUE, breaks= 20)
    points(pch= 17, col= 'blue', x= 0.05, y= 0)
    legend('topleft', legend=c('Pval +ve %', 'Pval -ve %'), col= c('red', 'grey70'), pch= 15, bg= 'white')
dev.off()
system(sprintf('scp hist_pvals_glmbn.paired.pdf $mac_office:$cri_public_projects/20120522-oxbsseq-mesc/Documents/slx-6278/genome/%s', CHROM))

quit(save= 'no')
