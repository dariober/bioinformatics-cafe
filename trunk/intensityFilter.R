# -----------------------------------------------------------------------------
# Functions to apply intensity dependent filter.
# Typically applied to differential expression output (edgeR, DEseq) where
# large fold changes are seen at low intensity genes.
#
# See SEQanswers, post #34 in http://seqanswers.com/forums/showthread.php?t=9998
#
## Example usage:
## Assuming a table of differentially expressed genes (logFC) and intensities (logCPM)
## has been generated with edgeR using topTags():
#
# detable<- topTags(lrt, n= nrow(d))$table
# 
"
source('~/svn_checkout/bioinformatics-misc/intensityFilter.R')
set.seed(1234)
logCPM<- rnorm(n= 10000, mean= 5, sd= 1)
set.seed(12345)
logFC<- rnorm(n= 10000, mean= 0, sd= 0.1)^3
smoothScatter(logCPM, logFC, nrpoints= 1000)
z3<- localZ.3(logCPM, logFC, nbins= 10)
z<- localZ.2(logCPM, logFC, nbins= 10)

points(logCPM, logFC, col= ifelse(abs(z3) > 2, 'red', NA), pch= 19)
"
## Calculate z-score. See below differece between localZ and localZ.2:
# detable$zScore<- localZ(detable$logCPM, detable$logFC, nbins= 20)
# -----------------------------------------------------------------------------

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


## -----------------------------------------------------------------------------
## DEPRECATED
## -----------------------------------------------------------------------------

#localZ.2<- function(x, y, nbins= 10){
#    # Same purpose and args as localZ but using a moving window centered around each
#    # data point. I.e. z-score for gene a is computed using the genes surrounding a.
#    # It is more accurate than localZ but much slower (half minute for 15000
#    # datapoints).
#    # -------------------------------------------------------------------------
#    ## Number of datapoints surrounding the target point
##    n<- round(length(x) / nbins, 0) 
#    xrank<- rank(x)
#    zx<- vector(length= length(y))
#    for (i in 1:length(xrank)){
#        if(i %% 1000 == 0){
#            print(i)
#        }
#        slice<- y[getLocalSlice(x= xrank, target= xrank[i], n= n)]
#        zscore<- (y[i] - mean(slice))/sd(slice)
#        zx[i]<- zscore
#    }
#    return(zx)
#}
#
#localZ.quantile<- function(x, y, nbins= 10){
#    # Return z-score of y for each bin along x
#    # x:
#    #   Vector of intensities. Typically logCPM (x-axis in MA plot). 
#    # y:
#    #   Vector to be normalized. Typically logFC (y-axis in MA plot).
#    # nbins:
#    #   Divide the x vector into this many bins, each of them containing the
#    #   same number of datapoints. z-scores are calculated within each bin.
#    # --------------------------------------------------------------------------
#    if (!is.numeric(x)){
#        stop("Vector x must be numeric")
#    }
#    if (!is.numeric(y)){
#        stop("Vector y must be numeric")
#    }
#    if (!is.numeric(nbins) | length(nbins) != 1 | nbins < 2){
#        stop("nbins must be one single integer equal or greater than 2")
#    }
#    if(length(x) != length(y)){
#        stop("Length of x and y differ")
#    }
#    qq<- quantile(x, p= seq(0, 1, length.out= nbins))
#    zx<- vector(length= length(y))
#    for( i in 1:(length(qq)-1) ){
#        if ( i == (length(qq)-1) ){
#            qidx<- which(x >= qq[i])
#        }
#        else {
#            qidx<- which(x >= qq[i] & x < qq[i+1])
#        }
#        zx[qidx]<- z.score(y[qidx])
#    }
#    if (length(zx) != length(y)){
#        stop("Unexpected length of z-score vector")
#    }
#    return(zx)
#}
#
#
#getLocalSlice<- function(x, target, n= NULL){
#    # Get datapoints in x surrounding the datapoint target
#    # x:
#    #   Vector of datapoint rank to slice. Typically the intensity or rank(logCPM).
#    # target:
#    #   The datapoint rank in x around which to get the slice
#    # n:
#    #   The number of x datapoints around target to return. Default to
#    #   length(x)/10. 
#    #
#    # Return: Logical vector of the same length as x with TRUE for the elements
#    # of x included in the slice.
#    # getLocalSlice(x= 1:100, target= 20, n= 10)
#    #
#    # TODO:
#    # Pass max(x) to arg. No need to compute max for each gene (length(x) very quick)
#    # Need to sort x and access by index slicing x[i:j]
#    # --------------------------------------------------------------------------
#    if (n > max(x)){ ## Quick for 15k slow for 150k. Possible avoid
#        stop("Number of local datapoints (n) greater than max rank")
#    }
#    if (target %in% x == FALSE){  ## c.ca. 11sec for 15k iters, AVOID!
#        stop("Datapoint target not in vector x")
#    }
#    if(is.null(n)){  ## Neglegible
#        n= length(x)/10
#    }
#    n<- round(n, 0)
#    xleft<- target - floor(n/2)
#    xright<- target + ceiling(n/2)
#    if (xleft < 1){
#        xright <- n + 1
#        xleft<- 1
#    } else if(xright > max(x)) {
#        xleft<- max(x) - n
#        xright<- max(x)
#    }
#    slice<- x %in% xleft:xright  ## 12s for 15k
#    return(slice)
#}
