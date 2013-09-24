# -----------------------------------------------------------------------------
# Functions to apply intensity dependent filter.
# Typically applied to differential expression output (edgeR, DEseq) where
# large fold changes are seen at low intensity genes.
#
# See SEQanswers, post #34 in http://seqanswers.com/forums/showthread.php?t=9998
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
    #
    #
    #
    xrank<- rank(x)
    zx<- vector(length= length(y))
    for (i in xrank){
        slice<- y[getLocalSlice(x= xrank, target= xrank[i], n= nbins)]
        zscore<- (y[i] - mean(slice))/sd(slice)
        zx[i]<- zscore
    }
}

localZ<- function(x, y, nbins= 10){
    # Return z-score of y for each bin along x
    # x:
    #   Vector of intensities. Typically logCPM (x-axis in MA plot). 
    # y:
    #   Vector to be normalized. Typically logFC (y-axis in MA plot).
    # nbins:
    #   Divide x vector into this many bins, each of them containing the
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
