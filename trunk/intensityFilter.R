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
}

intensityFilter(x, y, nbins= 10){
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
    qq<- quantile(x, p= seq(0, 1, length.out= nbins))
    zx<- vector(length= length(x))
    for( i in 1:(length(qq)-1) ){
        if ( i == (length(qq)-1) ){
            qidx<- which(x >= qq[i])
        }
        else {
            qidx<- which(x >= qq[i] & x < qq[i+1])
        }
        zx[qidx]<- z.score(x[qidx])
    }
    if (length(zx) != length(x)){
        stop("Unexpected length of z-score vector")
    }
    return(zx)
}
