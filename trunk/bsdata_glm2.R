#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

description<- 'Read matrix files for chrom and produce GLM output. Example bsdata_glm2.R -d design.txt -c BS oxBS  -i chr1 -m 5 -M 1000 -l 4'

# create parser object
parser <- ArgumentParser(description= description)

parser$add_argument("-i", "--indir", required= TRUE,
    help="Directory with input files")

parser$add_argument("-b", "--basename", default= 'genome',
    help="Basename to find input files inside --indir. Default 'genome'")

parser$add_argument("-d", "--design", required= TRUE,
    help="File with the experimental design")

parser$add_argument("-c", "--contrast", required= TRUE,
    nargs= 2,
    help= "Contrasts to apply")

parser$add_argument("-m", "--minCov", default= 0,
    type= "integer",
    help= "Minimum coverage for a position to be retained")

parser$add_argument("-M", "--maxCov", default= 1000000,
    type= "integer",
    help= "Max coverage for a position to be retained")

parser$add_argument("-l", "--minLibs", default= 1,
    type= "integer",
    help= "Min number of libraries that must satisfy -m and -M. Default 1: It's enough that one lib has >= minCov and <= maxCov for the locus to be included.")

parser$add_argument("-n", "--nprocs", default= 1,
    type= "integer",
    help= "Number of processors to use")

parser$add_argument("--nrows", default= -1,
    type= "integer",
    help= "For debugging: Number of rows from input to process. Default -1 (all rows)")

args <- parser$parse_args()

# ------------------------- Setting up ----------------------------------------
source('/home/berald01/svn_checkout/bioinformatics-misc/glm.R')
source('/home/berald01/svn_checkout/bioinformatics-misc/BSdata.R')
source('/home/berald01/svn_checkout/bioinformatics-misc/makeTransparent.R')
# -----------------------------------------------------------------------------

## Read design
## -----------
design<- read.table(args$design, header= TRUE, sep= '\t')
if(('library_id' %in% names(design)) == FALSE){
    stop('design file must contain a column "library_id"')
}
if(('bs' %in% names(design)) == FALSE){
    stop('\n\ndesign file must contain a column "bs"')
}
if(!all(args$contrast %in% design$bs)){
    stop('\nInvalid contrasts')
}
if(length(unique(args$contrast)) != 2){
    stop('Invalid number of contrasts')
}
if(args$minLibs > nrow(design)){
    stop('More libraries in minLibs than available.')
}

## Read data
## ---------
cpg<- read.bsdata(file.path(args$indir, args$basename), mat= c('loci', 'cnt_met', 'tot_reads'), gzip= FALSE, nrows= args$nrows)

if(all(as.character(design$library_id) == as.character(cpg@design$library_id)) == FALSE){
    print(design)
    stop('\n\nLibraries IDs in design file do not match those in input files')
}
cpg@design$bs<- design$bs

cpg@pct_met<- as.ffdf(100 * (cpg@cnt_met[,] / cpg@tot_reads[,]))

# Read count
# ----------
cat('Plot of read counts...\n')
pdf(file.path(args$indir, 'tot_reads.hist.pdf'), width= 15/2.54, height= 12/2.54, pointsize= 9)
par(mfrow= c(2, 3), mar= c(2.5, 3, 2, 0.2), oma= c(3,3,2,0))
for(x in colnames(cpg@tot_reads)){
    readQ<- quantile(cpg@tot_reads[, x], 0.95)
    y<- ifelse(cpg@tot_reads[, x] > readQ, readQ, cpg@tot_reads[, x] )
    hist(y, main= '', col= 'grey80', xlab= '', ylab= '', breaks= 10)
    mtext(side= 3, line= -1, text= x, cex= 0.8)
}
mtext(side= 1, text= 'Read depth', font= 2, outer= TRUE)
mtext(side= 2, text= 'Frequency', outer= TRUE, font= 2)
mtext(side= 3, text= sprintf('%s Depth of coverage at CpG sites', args$indir), outer= TRUE, font= 2, line= 0)
dev.off()

# Methylation 
# -----------
cat('Plot of methylation...\n')
## Select CpG with coverage n<x<m
m<- 5
Q<- 0.95

pdf(file.path(args$indir, 'pct_met.hist.pdf'), width= 15/2.54, height= 12/2.54, pointsize= 9)
par(mfrow= c(2, 3), mar= c(2.5, 3, 2, 0.2), oma= c(3,3,2,0))
    for(x in colnames(cpg@pct_met)){
        y<- cpg@pct_met[, x]
        tot<- cpg@tot_reads[, x]
        readQ<- quantile(cpg@tot_reads[, x], Q)
        y<- y[which(tot >= m & tot <= readQ)]  
        hist(y, main= '', col= 'grey80', xlab= '', ylab= '', breaks= 25)
        mtext(side= 3, line= -1, text= paste(x, '\nn= ', length(y), sep= ''), cex= 0.8, col= 'darkblue')
    }
mtext(side= 1, text= '% methylated', font= 2, outer= TRUE)
mtext(side= 2, text= 'Frequency', outer= TRUE, font= 2)
mtext(side= 3, text= sprintf('%s CpG sites with depth %sx to quantile(%s)', args$indir, m, Q), outer= TRUE, font= 2, line= 0)
dev.off()

# Filter positions with not-too-high coverage
# -------------------------------------------

# Select correct columns
cpg.2<- BSdataApply(cpg, FUN= function(x) x[, which(cpg@design$bs %in% args$contrast)])
cpg.2@design<- cpg@design[which(cpg@design$bs %in% args$contrast),]

# Filter for coverage
filter<- which(rowSums(cpg.2@tot_reads[,] >= args$minCov & cpg.2@tot_reads[,] <= args$maxCov) >= args$minLibs) ## Positions with coverage x in at least y libs 
cat(sprintf('Applying GLM to %s loci\n', length(filter)))
cpg.2<- BSdataApply(cpg.2, FUN= function(x) x[filter,])

pdf(file.path(args$indir, 'hist_depth_comb_libs.pdf'), width= 18/2.54, height= 10/2.54, pointsize= 10)
par(mfrow= c(1,2))
hist(rowSums(cpg.2@tot_reads[,]),
    main= sprintf('Depth of coverage summed across %s libraries\nn= %s', ncol(cpg.2@tot_reads), nrow(cpg.2@tot_reads)), 
    xlab= 'Cumulated depth of coverage', ylab= 'No. CpG', breaks= 20)
hist(rowMeans(cpg.2@pct_met[,]),
    main= sprintf('%s Methylation averaged across %s libraries\nn= %s', args$indir, ncol(cpg.2@tot_reads), nrow(cpg.2@tot_reads)), 
    xlab= '% methylation', ylab= 'No. CpG', breaks= 20)
dev.off()

# GLM 
# ---

glmbs<- glmBS.2(bsobj= cpg.2, bs= cpg.2@design$bs, contrast= args$contrast, family= 'binomial', nthreads= args$nprocs)
row.names(glmbs)<- row.names(cpg.2@cnt_met)
write.table(glmbs, file= file.path(args$indir, 'glmbs.txt'), row.names= FALSE, col.names= TRUE, sep= '\t', quote= FALSE)

## MEMO: Estimate is "oxBS - BS" (reverse alphanumeric order)
pneg<- glmbs$pvalue[which(glmbs$estimate > 0)] ## Pvalues of negative estimates.
ppos<- glmbs$pvalue[which(glmbs$estimate < 0)]

pdf(file.path(args$indir, 'hist_pvals_glm.pdf'), width= 18/2.54, height= 10/2.54, pointsize= 8)
par(las= 0, mfrow= c(1, 2))
    h1<- hist(pneg, border= 'transparent', add= FALSE, col= 'grey70', freq= FALSE, breaks= 20, xlab= 'P-value', main= "P-values from GLM binom")
    h2<- hist(ppos, border= 'red', freq= FALSE, add= TRUE, breaks= 20)
    points(pch= 17, col= 'blue', x= 0.05, y= 0)
    
    h1<- hist(pneg, border= 'transparent', add= FALSE, col= 'grey70', freq= TRUE, breaks= 20, xlab= 'P-value', main= sprintf("%s P-values from GLM binom\nn= %s", args$indir, nrow(glmbs)))
    h2<- hist(ppos, border= 'red', freq= TRUE, add= TRUE, breaks= 20)
    points(pch= 17, col= 'blue', x= 0.05, y= 0)
    legend('topleft', legend=c('Pval +ve %', 'Pval -ve %'), col= c('red', 'grey70'), pch= 15, bg= 'white')
dev.off()

quit(save= 'no')
