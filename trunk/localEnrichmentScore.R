#! /usr/bin/env Rscript

done<- suppressWarnings(suppressMessages(require(argparse)))
if(done == FALSE){
    cat('Please install the "argparse" package. Open an R session and execute:\n')
    cat('install.packages("argparse")\n\n')
    quit(save= 'no', status= 1)
}
done<- suppressWarnings(suppressMessages(require(data.table)))
if(done == FALSE){
    cat('Please install the "data.table" package. Open an R session and execute:\n')
    cat('install.packages("data.table")\n\n')
    quit(save= 'no', status= 1)
}

done<- suppressWarnings(suppressMessages(require(aod)))
if(done == FALSE){
    cat('Please install the "data.table" package. Open an R session and execute:\n')
    cat('install.packages("data.table")\n\n')
    quit(save= 'no', status= 1)
}

VERSION<- '0.1.0a'

docstring<- sprintf("DESCRIPTION \\n\\
Compare local enrichment between pairs of files, e.g. condition 1 vs 2, to \\n\\
compute log fold change and associated significance. Alternatively, compute the \\n\\
overal significance of enrichment in a set of input files\\n\\n\\
\\
Enrichment is calculated for each pair first by comparing the number of reads \\n\\
in target and in background by chi^2 test. \\n\\n\\
\\
A final pvalue is obtained by combining p-values from each pair via Stouffer method. \\n\\
Input files are typically produced by localEnrichmentBed.py. \\n\\n\\
NB: The direction of change tested is always [input - control] if this \\n\\
difference is < 0 than the chi^2 pvalue is returned as negative. \\n\\n\\
\\
OUTPUT \\n\\
1) rank.bed file with tested regiosn ranked by rank-product, combined p-value \\n\\
(as -log10) and average log2fc. Note that log2fc is the mean of the individual \\n\\
logFCs.\\n\\
2) score.bed file with enrichment for each pair of files and rank. This file \\n\\
produced only if control file(s) are given. \\n\\n\\
Version %s", VERSION)


parser<- ArgumentParser(description= docstring, formatter_class= 'argparse.RawTextHelpFormatter')

parser$add_argument("-i", "--input", help= "Input files generate by localEnrichmentBed.py", nargs= '+', required= TRUE)
parser$add_argument("-c", "--ctrl", help= "Control file (e.g. the input). It must be one file (for all inputs) or one for each input", nargs= '+')
parser$add_argument("-o", "--out", help= "Basename for output files", required= TRUE)
parser$add_argument("-m", "--method", help= "_Not implemented yet_ Method to detect differences", default= 'chisq')

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
        # Recompute in log space
        p.val<- Stouffer.log(p, w)[['p.value']]
    } else {
        p.val <- 1-pnorm(Z)
    }
    return(p.val)
}

Stouffer.log<- function(p, w) {
    # Stouffer method modified for log space. Don't use directly,
    # call Stouffer.comb() instead
    if (missing(w)) {
        w <- rep(1, length(p))/length(p)
    } else {
        if (length(w) != length(p))
        stop("Length of p and w must equal!")
    }
    logp<- log(p)     ## log transform pvalues
    Zi <- -qnorm(logp, log.p= TRUE) ## Work in log space
    Z  <- sum(w*Zi)/sqrt(sum(w^2))
    p.val <- -pnorm(Z, log.p= TRUE)
    out.p<- c(Z= Z, p.value= p.val)
    return( out.p )
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
            if(sum(m) > 0){
                log2fc[j]<- log2( (m[1,1] / m[2,1]) / (m[1,2] / m[2,2]) )
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
    enrich[, ranked := list(rank(-log2fc)), by= input]

#    newp<- ifelse(enrich$log2fc < 0, 1, enrich$pval)  ## Set -ve pval to 1
#    newp<- ifelse(newp < 1e-16, 1e-16, newp)          ## Reset pval if too close to 0 or 1.
#    newp<- ifelse(newp > 0.999999, 0.999999, newp)
    enrich$pval<- newp
    rm(newp)
    
    rankProdDT<- enrich[, list(avgLog2fc= mean(log2fc), log10_combp= -log10(Stouffer.comb(pval)), RP= gmean(ranked)), by= targetID]
    rankProdDT<- merge(loci, rankProdDT, by= 'targetID')
    rankProdDT<- rankProdDT[, list(chrom, start, end, targetID, avgLog2fc, log10_combp, RP) ]
    enrich[, pval := NULL]
    
    write.table(enrich, file= paste(xargs$out, 'score.bed', sep= '.'), row.names= FALSE, col.names= TRUE, sep= '\t', quote= FALSE)
    write.table(rankProdDT, file= paste(xargs$out, 'rank.bed', sep= '.'), row.names= FALSE, col.names= TRUE, sep= '\t', quote= FALSE)
} else {
    cat(sprintf('\nInvalid method! Got "%s"\n\n', xargs$method))
}
quit(save= 'no')