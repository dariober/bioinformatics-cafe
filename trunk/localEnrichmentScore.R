#! /usr/bin/env Rscript

done<- suppressWarnings(suppressMessages(require(argparse)))
if(done == FALSE){
    cat('\nPlease install the "argparse" package. Open an R session and execute:\n\n')
    cat('> install.packages("argparse")\n\n')
    cat('Once you are at it, install also the data.table package, if not already installed:\n\n')
    cat('> install.packages("data.table")\n\n')
    quit(save= 'no', status= 1)
}
done<- suppressWarnings(suppressMessages(require(data.table)))
if(done == FALSE){
    cat('Please install the "data.table" package. Open an R session and execute:\n')
    cat('install.packages("data.table")\n\n')
    quit(save= 'no', status= 1)
}

#done<- suppressWarnings(suppressMessages(require(aod)))
#if(done == FALSE){
#    cat('Please install the "data.table" package. Open an R session and execute:\n')
#    cat('install.packages("data.table")\n\n')
#    quit(save= 'no', status= 1)
#}

VERSION<- '0.1.0a'

docstring<- sprintf("DESCRIPTION \\n\\
Compare local enrichment between pairs of files, e.g. condition 1 vs 2, to \\n\\
compute log fold change and associated significance. Alternatively, if only one \\n\\
condition is given, compute the overal significance of enrichment in a set of \\n\\
input files.\\n\\n\\
\\
Enrichment is calculated for each pair first by comparing the number of reads \\n\\
in target and in background by chi^2 test like this: \\n\\n\\
\\
        input | ctrl \\n\\
-------+------+----- \\n\\
target |  1,1 | 1,2  \\n\\
-------+------+----- \\n\\
 flank |  2,1 | 2,2  \\n\\
\\
\\n\\
\\
A final pvalue is obtained by combining p-values from each pair via Stouffer method. \\n\\
Input files are typically produced by localEnrichmentBed.py. \\n\\n\\
NB: The direction of change tested is always [input - control] if this \\n\\
difference is < 0 than the chi^2 pvalue is returned as negative. \\n\\n\\
\\
OUTPUT \\n\\
1) A *.rank.bed file with tested regions ranked by rank-product, combined p-value \\n\\
(as -log10) and average log2fc. Note that log2fc is the mean of the individual \\n\\
logFCs.\\n\\
2) A *.score.bed file with enrichment for each pair of files and rank. This file \\n\\
produced only if control file(s) are given. \\n\\n\\
EXAMPLE \\n\\
\\n\\
## Compare two pairs of pull-downs and controls:\\n\\
localEnrichmentScore.R -i BG4-1.leb.bed  BG4-2.leb.bed -c ctrl-1.leb.bed  ctrl-2.leb.bed -o bg4-ctrl \\n\\
\\n\\
## Combine evidence of enrichment in two replicates:\\n\\
localEnrichmentScore.R -i BG4-1.leb.bed  BG4-2.leb.bed -o bg4-ctrl \\n\\
\\n\\
Version %s", VERSION)


parser<- ArgumentParser(description= docstring, formatter_class= 'argparse.RawTextHelpFormatter')

parser$add_argument("-i", "--input", help= "Input files generate by localEnrichmentBed.py", nargs= '+', required= TRUE)
parser$add_argument("-c", "--ctrl", help= "Control file (e.g. the input). It must be one file (for all inputs) or one for each input", nargs= '+')
parser$add_argument("-o", "--out", help= "Basename for output files", required= TRUE)
parser$add_argument("-ps", "--pseudocnt", help= "Pseudo-counts to add to each cell to avoid zeros", type= 'integer', default= 0)
parser$add_argument("-r", "--rankby", help= "Rank by this variable. Default log10_pval", default= 'log10_pval', choices= c("log10_pval", "log2fc"))
parser$add_argument("-m", "--method", help= "_Not implemented yet_ Method to detect differences (only chisq supported.)", default= 'chisq')

xargs<- parser$parse_args()

# -----------------------------------------------------------------------------

gmean<- function(x){
    exp(mean(log(x)))
}

Stouffer.comb <- function(p, w) {
    # Modified from Wikipedia
    # p: is a vector of p-values
    # w: Vector of weights
    if (missing(w)) {
    w <- rep(1, length(p))/length(p)
    } else {
    if (length(w) != length(p))
       stop("Length of p and w must equal!")
    }
    Zi <- qnorm(1-p) 
    Z  <- sum(w*Zi)/sqrt(sum(w^2))
    if(!is.na(Z) & Z > 7){
        p.val <- -pnorm(Z, log.p= TRUE)
    } else {
        p.val <- 1-pnorm(Z)
    }
    return(p.val)
}

log10pvalToPvalForLogFC<- function(log2fc, log10pval, pmin= 1e-16, pmax= 0.999999){
    ## Convert pvalues from -log10() format to original scale 0 to 1.
    ## If log2fc is -ve make the pvalue 1 (or close to), meaning that
    ## there is no evidence at all of enrichment.
    ## log2fc: Vector of logFC
    ## log10pval: Vector of -log10(pvalues)
    ## pmin: Cap min pvalue: If a pvalue is lowewr than this, reset it to pmin.
    ## pmax: Cap max pvalue: If a pvalue is more than this, reset it.
    ## Example:
    ## log10pvalToPvalForLogFC(c(1, 2, -1), c(2, 3, 1))
    ## [1] 0.010000 0.001000 0.999999
    ## 
    stopifnot(length(log2fc) == length(log10pval))
    stopifnot(pmin >= 0 & pmin <= 1)
    stopifnot(pmax >= 0 & pmax <= 1)
    stopifnot(pmax > pmin)
    stopifnot(log10pval >= 0 | is.na(log10pval)) ## NAs are not explicitly handled!!
    pval<- 10^(-log10pval)
    newp<- ifelse(log2fc < 0, pmax, pval)  ## Set -ve pvalue to max pval
    newp<- ifelse(newp < pmin, pmin, newp) ## Reset pval if too close to 0 or 1.
    newp<- ifelse(newp > pmax, pmax, newp)
    return(newp)
}

z.score<- function(x){
    ## Return z-scores for vector x
    z<- (x - mean(x))/sd(x)
    return(z)
}

localZ<- function(x, y, nbins= 10){
    # Return z-score of each y for each window around x
    # x:
    #   Vector of intensities. Typically logCPM (x-axis in MA plot). 
    # y:
    #   Vector to be z-score'd. Typically logFC (y-axis in MA plot).
    # nbins:
    #   Divide the x vector into this many bins, each of them containing the
    #   same number of datapoints. z-scores are calculated within each bin.
    #
    # See google.code/bioinformatics-misc for this function.
    # -------------------------------------------------------------------------
    ## Number of datapoints surrounding the target point
    n<- round(length(x) / nbins, 0) 

    ## Record the order of input x (and y)
    xorder<- order(x)
    
    ## Sort y by x
    yorder<- y[xorder]

    ## Initialize vector of z-scores zvec
    zvec<- vector(length= length(x))
    ## For each point p in x,
    xleft<- floor(n/2)
    xright<- ceiling(n/2)
    for (i in 1:length(y)){
        if (i <= xleft){
            slice<- yorder[1:n]
        } else if ((i + xright) >= length(x)) {
            slice<- yorder[((length(x) - n)+1) : length(x)]
        } else {
            slice<- yorder[((i - xleft)+1) : (i + xright)]
        }
        zvec[i]<- (yorder[i] - mean(slice)) / sd(slice)
    }
    ## Return zvec in the original order of x (and y)
    zvec_ordered<- zvec[order(xorder)]
    return(zvec_ordered)
}

# -----------------------------------------------------------------------------
if(xargs$method == 'chisq'){
    if(length(xargs$ctrl) != 0 & length(xargs$ctrl) != 1 & length(xargs$ctrl) != length(xargs$input)){
        cat('\nPlease provide either one control file or as many control files (possibly repeated) as inputs.\n\n')
        quit(save= 'no', status= 1)
    }
    if(length(xargs$ctrl) == 1){
        controls<- rep(xargs$ctrl, length(xargs$input))
    } else {
        controls<- xargs$ctrl
    }
}
# Read data
# =============================================================================
lebfiles<- c(xargs$input, unique(xargs$ctrl))
leb<- fread(lebfiles[1])
loci<- leb[, list(chrom, start, end, targetID)]
leb$library_id<- lebfiles[1]
for(f in lebfiles[2:length(lebfiles)]){
    tmp<- fread(f)
    if(all(loci$targetID == tmp$targetID) == FALSE){
        cat('\nInput files do not have the same loci!\n\n')
        quit(save= 'no', status= 1)
    }
    tmp$library_id<- f
    leb<- rbindlist(list(leb, tmp))
}
rm(tmp)

## Add pseudocounts
leb[, target_cnt := target_cnt + xargs$pseudocnt]
leb[, flank_cnt := flank_cnt + xargs$pseudocnt]

leb[, type := ifelse(library_id %in% xargs$ctrl, 'ctrl', 'pd')]
# =============================================================================
# Only input files
if(is.null(xargs$ctrl)){
    ## Do not compare, just give a combined pvalue and logFC
    ## Do not create the score.bed file as this would be just the input files
    ## concatenated.
    ## Make sure the output format is consistent with the chisq output below
    leb[, ranked := list(rank(-log2fc)), by= library_id]
    leb$pval<- log10pvalToPvalForLogFC(leb$log2fc, leb$log10_pval)
    rankProdDT<- leb[, list(avgLog2fc= mean(log2fc), log10_combp= -log10(Stouffer.comb(pval)),
        RP= gmean(ranked)), by= list(chrom, start, end, targetID)]
    write.table(rankProdDT, file= paste(xargs$out, 'rank.bed', sep= '.'),
        row.names= FALSE, col.names= TRUE, sep= '\t', quote= FALSE)
    quit(save= 'no')
}

if(xargs$method %in% c('binomial', 'betabin')){
    
    cat('\n\n Sorry not implemented yet!\n')
    quit(save= 'n')

    # GLM
    # ==========================================================================
    log10_pval<- rep(NA, nrow(loci))
    for(i in 1: nrow(loci)){
        tid<- loci$targetID[i]
        dat<- leb[targetID == tid, ]
        if(xargs$method == 'binomial'){
            gp<- summary(glm(cbind(target_cnt, flank_cnt) ~ type, data= dat, family= binomial()))
            log10_pval[i]<- -log10(gp$coefficients[2,4])
        } else if (xargs$method == 'betabin') {
            gp<- summary(betabin(cbind(target_cnt, flank_cnt) ~ type, ~1, data= dat))
            log10_pval[i]<- -log10(gp@Coef[2,4])
        }
    }
    loci[, log10_pval := log10_pval]
    write.table(loci, file= stdout(), row.names= FALSE, col.names= TRUE, sep= '\t', quote= FALSE)
    quit()
} else if(xargs$method == 'chisq'){
    # Reshape to matrices
    # ==========================================================================
    target<- dcast.data.table(leb, targetID ~ library_id, value.var= 'target_cnt')
    flank<- dcast.data.table(leb, targetID ~ library_id, value.var= 'flank_cnt')
    stopifnot(target$targetID == flank$targetID)
    stopifnot(target$targetID == loci$targetID)
    target<- as.matrix(target[, targetID := NULL])
    rownames(target)<- loci$targetID
    flank<- as.matrix(flank[, targetID := NULL])
    rownames(flank)<- loci$targetID
    stopifnot(colnames(target) == colnames(flank))
    cnames<- colnames(target)
    
    # Evaluate enrichment
    # ==========================================================================
    first<- TRUE
    for(i in 1:length(controls)){
        ## For each pair of files compare target_cnt and falnk_cnt in input vs ctrl
        ctrl<- controls[i]
        ctrli<- which(cnames == ctrl)
        pd<- xargs$input[i]
        pdi<- which(cnames == pd)    
        pval<- rep(NA, nrow(target))
        log2fc<- rep(NA, nrow(target))
        for(j in 1:nrow(target)){
            m<- matrix(c(target[j, pdi], target[j, ctrli],
                         flank[j, pdi], flank[j, ctrli]), nrow= 2, byrow= TRUE)
            colnames(m)<- c(pd, ctrl)
            rownames(m)<- c('target', 'flank')
            #
            # Matrix m is a contingency table like this:
            #         input | ctrl
            # -------+------+-----
            # target |  1,1 | 1,2 
            # -------+------+-----
            #  flank |  2,1 | 2,2
            #
            if(sum(m) > 0){
                log2fc[j]<- log2( (m[1,1] / m[2,1]) /
                                  (m[1,2] / m[2,2])
                                )
                pval[j]<- chisq.test(m)$p.value
            } else {
                pval[j]<- NA
                log2fc[j]<- NA
            }
        }
        tmp<- data.table(loci, input= pd, control= ctrl, pval= pval, log2fc)
        if(first){
            enrich<- tmp
            first<- FALSE
        } else {
            enrich<- rbindlist(list(enrich, tmp))
        }
    }
    enrich[, log10_pval := -log10(pval)]

    ## For combining pvals:
    newp<- log10pvalToPvalForLogFC(enrich$log2fc, enrich$log10_pval)
    
    enrich[, log10_pval := ifelse(log2fc < 0, -log10_pval, log10_pval)] ## Change sign of pvals
    if(xargs$rankby == 'log10_pval'){
        enrich[, ranked := list(rank(-log10_pval)), by= input]
    } else if (xargs$rankby == 'log2fcl'){
        enrich[, ranked := list(rank(-log2fc)), by= input]    
    } else {
        stop(sprintf("Unsopported option to rank: '%s'", xargs$rankby))
    }
#    newp<- ifelse(enrich$log2fc < 0, 1, enrich$pval)  ## Set -ve pval to 1
#    newp<- ifelse(newp < 1e-16, 1e-16, newp)          ## Reset pval if too close to 0 or 1.
#    newp<- ifelse(newp > 0.999999, 0.999999, newp)
    enrich$pval<- newp
    rm(newp)
    
    rankProdDT<- enrich[, list(avgLog2fc= mean(log2fc), log10_combp= -log10(Stouffer.comb(pval)), RP= gmean(ranked)), by= targetID]
    rankProdDT<- merge(loci, rankProdDT, by= 'targetID')
    rankProdDT<- rankProdDT[, list(chrom, start, end, targetID, avgLog2fc, log10_combp, RP) ]
    enrich[, pval := NULL]

    options(scipen= 9) # This is to make ranks as integers    
    write.table(enrich, file= paste(xargs$out, 'score.bed', sep= '.'), row.names= FALSE, col.names= TRUE, sep= '\t', quote= FALSE)
    write.table(rankProdDT, file= paste(xargs$out, 'rank.bed', sep= '.'), row.names= FALSE, col.names= TRUE, sep= '\t', quote= FALSE)
} else {
    cat(sprintf('\nInvalid method! Got "%s"\n\n', xargs$method))
}
warnings()
quit(save= 'no')