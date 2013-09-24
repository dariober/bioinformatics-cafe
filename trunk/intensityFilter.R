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
## Calculate z-score. See below differece between localZ and localZ.2:
# detable$zScore<- localZ(detable$logCPM, detable$logFC, nbins= 20)
# detable$zScore<- localZ.2(detable$logCPM, detable$logFC, nbins= 20)
# -----------------------------------------------------------------------------

z.score<- function(x){
    ## Return z-scores for vector x
    z<- (x - mean(x))/sd(x)
    return(z)
}

getLocalSlice<- function(x, target, n= NULL){
    # Get datapoints in x surrounding the datapoint target
    # x:
    #   Vector of datapoint rank to slice. Typically the intensity or rank(logCPM).
    # target:
    #   The datapoint rank in x around which to get the slice
    # n:
    #   The number of x datapoints around target to return. Default to
    #   length(x)/10. 
    #
    # Return: Logical vector of the same length as x with TRUE for the elements
    # of x included in the slice.
    # getLocalSlice(x= 1:100, target= 20, n= 10)
    # --------------------------------------------------------------------------
    if (n > max(x)){
        stop("Number of local datapoints (n) greater than max rank")
    }
    if (target %in% x == FALSE){
        stop("Datapoint target not in vector x")
    }
    if(is.null(n)){
        n= length(x)/10
    }
    n<- round(n, 0)
    xleft<- target - floor(n/2)
    xright<- target + ceiling(n/2)
    if (xleft < 1){
        xright <- n + 1
        xleft<- 1
    } else if(xright > max(x)) {
        xleft<- max(x) - n
        xright<- max(x)
    }
    slice<- x %in% xleft:xright
    return(slice)
}

localZ.2<- function(x, y, nbins= 10){
    # Same purpose and args as localZ but using a moving window centered around each
    # data point. I.e. z-score for gene a is computed using the genes surrounding a.
    # It is more accurate than localZ but much slower (half minute for 15000
    # datapoints).
    # -------------------------------------------------------------------------
    ## Number of datapoints surrounding the target point
    n<- round(length(x) / nbins, 0) 
    xrank<- rank(x)
    zx<- vector(length= length(y))
    for (i in 1:length(xrank)){
        if(i %% 1000 == 0){
            print(i)
        }
        slice<- y[getLocalSlice(x= xrank, target= xrank[i], n= n)]
        zscore<- (y[i] - mean(slice))/sd(slice)
        zx[i]<- zscore
    }
    return(zx)
}

localZ<- function(x, y, nbins= 10){
    # Return z-score of y for each bin along x
    # x:
    #   Vector of intensities. Typically logCPM (x-axis in MA plot). 
    # y:
    #   Vector to be normalized. Typically logFC (y-axis in MA plot).
    # nbins:
    #   Divide the x vector into this many bins, each of them containing the
    #   same number of datapoints. z-scores are calculated within each bin.
    # --------------------------------------------------------------------------
    if (!is.numeric(x)){
        stop("Vector x must be numeric")
    }
    if (!is.numeric(y)){
        stop("Vector y must be numeric")
    }
    if (!is.numeric(nbins) | length(nbins) != 1 | nbins < 2){
        stop("nbins must be one single integer equal or greater than 2")
    }
    if(length(x) != length(y)){
        stop("Length of x and y differ")
    }
    qq<- quantile(x, p= seq(0, 1, length.out= nbins))
    zx<- vector(length= length(y))
    for( i in 1:(length(qq)-1) ){
        if ( i == (length(qq)-1) ){
            qidx<- which(x >= qq[i])
        }
        else {
            qidx<- which(x >= qq[i] & x < qq[i+1])
        }
        zx[qidx]<- z.score(y[qidx])
    }
    if (length(zx) != length(y)){
        stop("Unexpected length of z-score vector")
    }
    return(zx)
}
