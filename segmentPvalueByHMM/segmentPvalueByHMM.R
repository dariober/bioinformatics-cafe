#! /usr/bin/env Rscript

## TODO:
## * If a chrom is fully demethylated or methylated the state might be wrongly assigned!!
## * Print results as they are produced by each loop of HMMrunner()

VERSION<- '0.2.0'
done<- suppressWarnings(suppressMessages(require(RHmm)))
if(done == FALSE){
    cat('Please install the "RHmm" package.\n')
    quit(save= 'no', status= 1)
}
done<- suppressWarnings(suppressMessages(require(data.table)))
if(done == FALSE){
    cat('\nPlease install the "data.table" package. Open an R session and execute:\n\n')
    cat('> install.packages("data.table")\n\n')
    quit(save= 'no', status= 1)
}

## FUNCTIONS
## ---------

readBedFile<- function(filename, pIndex, chromIndex= 1, header= FALSE){
    # Read bed file where p-values to segment are in column index `pIndex` and chromosomes in 
    # column `chromIndex`
    write('Reading in data file... ', stderr())
    if(filename == '-'){
        bed<- data.table(read.table(file('stdin'), header= header, sep= '\t', stringsAsFactors= FALSE)[, c(chromIndex, pIndex)])
    } else if(grepl('\\.gz$', filename)){
        bed<- fread(sprintf('gunzip -c %s', filename), select= c(chromIndex, pIndex), header= header, sep= '\t', showProgress= FALSE)
    } else {
        bed<- fread(filename, select= c(chromIndex, pIndex), header= header, sep= '\t', showProgress= FALSE)
    }
    setnames(bed, names(bed)[1], 'chrom')
    setnames(bed, names(bed)[2], 'pvals')
    write(sprintf('Done: %s rows', nrow(bed)), stderr())
    return(bed)
}

recodeOnePvalue<- function(p, recodePval= c(0.1, 0.05, 0.001)){
    ## Recode pvalues to be discrete
    stopifnot(recodePval >= 0 & recodePval <= 1)
    stopifnot(sort(recodePval, decreasing= TRUE) == recodePval)
    stopifnot(length(recodePval) == 3)
    if (is.na(p)){
        o<- NA
    } else if (p >= recodePval[1] ){
        o<- 0
    } else if (p < recodePval[1] & p >= recodePval[2]) {
        o<- 1
    } else if (p < recodePval[2] & p >= recodePval[3]) {
        o<- 2
    } else if (p< recodePval[3]){
        o<- 3
    } else {
        stop('Unexpected p-value')
    }
    return(o)
}
recodePvals<- function(p, recodePval= c(0.1, 0.05, 0.001)){
    stopifnot(is.numeric(p))
    stopifnot(min(p, na.rm= TRUE) >= 0)
    stopifnot(max(p, na.rm= TRUE) <= 1)
    recoded<- sapply(p, recodeOnePvalue, recodePval=recodePval)
    return(recoded)    
}


prepareChromChunks<- function(chrSize, maxChunkSize){
    # Generate chunks of approximately equal size spanning chrSize.
    # Each chunks <= maxChunkSize
    # chrSize:
    #   The size (of the chrom) to be split
    # maxChunkSize:
    #   Max size of each chunks
    # Return:
    #   list with components:
    #       * starts: where chunks start
    #       * ends: where chunks end
    # Example:
    # prepareChromChunks(50, 20)
    # prepareChromChunks(20, 20)
    # prepareChromChunks(100, 33)
    #    $starts
    #    1 26 51 76
    #    $ends
    #    25  50  75 100
    #
    nchunks<- ceiling(chrSize / maxChunkSize)
    chunksize<- floor(chrSize  / nchunks)
    if (chrSize <= maxChunkSize){
        ## Chrom fits in one chunk
        starts<- 1
        ends<- chrSize
    } else {
        steps<- seq(0, chunksize*nchunks, by= chunksize) ## seq(1, nrow(Chrom), by= chunksize)
        starts<- steps[1 : (length(steps)-1)] + 1
        ends<- c(steps[2 : (length(steps)-1)], chrSize)
        ## TODO: Last chunk might be larger larger than maxChunkSize. E.g. prepareChromChunks(100, 17)
        stopifnot(all(starts[2:nchunks] - ends[1:(nchunks-1)] == 1)) ## Chunks are contiguous without gaps and overlaps
    }
    ## Some checks chunks are correct:
    stopifnot(nchunks == length(starts))
    stopifnot(
#        all( (ends - starts) <= maxChunkSize ), ## Chunks are small enough
        starts[1] == 1, ## Start from first position
        ends[nchunks] == chrSize, ## Finish end of chorm
        length(starts) == length(ends)
    )
    return(list(starts= starts, ends= ends))
}

recodeOneState<- function(states, obs, levels){
    # Recode the vector of unique states accoring to pvalues and recoding levels
    # If HMM gives a stretch of only one state (e.g. only demethylated) you
    # want to make this state low pval or high pval
    stopifnot(length(unique(states)) == 1)
    med<- median(obs, na.rm= TRUE)
    if(med > 0.1) {
        # Educated guess for median. If above this threshold call the state "high pval"
        newx<- rep(levels[length(levels)], length(states))
    } else {
        newx<- rep(levels[1], length(states))
    }
    return(newx)
}

recodeStates<- function(states, obs, levels){
    # Rename the HMM decoded states to match the given ones. Renaming is done by sorting
    # the decoded states by ascending obs
    # states:
    #   Vector of states to be renamed. Typically from output of viterbi(hmmfit, obs)$states
    # obs:
    #   Vector of observation. It will be used to get the mean of each state and sort
    #   states accordingly
    # levels:
    #   Vector of levels of length equal to the number of unique states to be used for
    #   renaming. Must be sorted by increasing pvalue. E.g. c('M', 'u') for M (low p),
    #   u (high p).
    # Return:
    #    Vector of renamed states
    # Example:
    #   states<- c(2,2,2,1,1,1,2,2)
    #   obs<- c(0.1, 0.15, 0.2, 0.9, 0.8, 0.85, 0.1, 0.2)
    #   levels<- c('L', 'H')
    #   recodeStates(states, obs, levels) #<<< c(L, L, L, H, H, H, L, L)
    # -------------------------------------------------------------------------
    stopifnot(length(states) == length(obs))
    if(length(unique(states)) == 1){
        newx<- recodeOneState(states, obs, levels)
    } else {
        avg<- data.table(states, obs)
        avg<- avg[, list(stateavg= mean(obs, na.rm= TRUE)), by= states]
        avg<- avg[order(stateavg), ]
        avg$newStates<- levels
        newx<- levels[match(states, avg$states)]
    }
    return(newx)
}

HMMChrom<- function(chromBed, MAXPOS= 40000){
    # Fit HMM thorugh chromosome.
    # Input chromBed is data.table with columns:
    stopifnot(c('pvals', 'recode_P') %in% names(chromBed))
    hmmChrom<- data.table(chrom= NA, pvals= NA, state= NA, postM= NA, state_id= NA)[0,]
    startsEnds<- prepareChromChunks(nrow(chromBed), MAXPOS)
    starts<- startsEnds$starts
    ends<- startsEnds$ends

    hmmNeeded<- TRUE
    if(length(chromBed$recode_P) < 500000){
        ## Fit the model to entire chrom, if not too big:
        hmmfit<- HMMFit(chromBed$recode_P, nStates= 2, dis= 'DISCRETE')
        hmmNeeded<- FALSE
    }
    for (j in 1:length(starts)){
        write(sprintf('from: %s to: %s length: %s', starts[j], ends[j], ends[j] - starts[j] + 1), stderr())
        dat<- chromBed[starts[j]:ends[j],]
        recode_P<- as.numeric(dat$recode_P)
        if(hmmNeeded){
            hmmfit<- HMMFit(recode_P, nStates= 2, dis= 'DISCRETE')
        }           
        ## DATA DECODING
        vit<- viterbi(hmmfit, recode_P)
     
        ## Recode to have low pval Methylated
        xstates<- recodeStates(vit$states, dat$pvals, levels= c('M', 'u'))  

        ## POSTERIOR DECODING: Get the probability of each cytosine to belong to the M or u state
        fbLog <- forwardBackward(hmmfit, dat$recode_P)
        prob<- fbLog$Gamma # Matrix of probs. Rows: Observations (loci), Cols: States
        # We need Need to pick the column corresponding to state 'M'.
        stopifnot(ncol(prob) == 2)
        if(xstates[1] == 'M'){
            if(vit$states[1] == 1){
                mx<- 1
            } else {
                mx<- 2
            }
        } else if(xstates[1] == 'u'){
            if(vit$states[1] == 1){
                mx<- 2
            } else {
                mx<- 1
            }
        } else {
            stop('Unexpected condition')
        }
        hmmOutSub<- data.table(chrom= dat$chrom, pvals= dat$pvals, state= xstates, postM= prob[,mx])
        xrle<- rle(hmmOutSub$state)
        state_id<- paste0(hmmOutSub$state, rep(1:length(xrle$values), times= xrle$length))
        hmmOutSub[, state_id := state_id]
        hmmChrom<- rbindlist(list(hmmChrom, hmmOutSub))
    }
    return(hmmChrom)
}

HMMrunner<- function(bed, recodePval= c(0.1, 0.05, 0.001)){
    # Run HMM. bed is a bed file with columns 'chrom' and 'pvals'
    bed[, recode_P := recodePvals(pvals, recodePval= recodePval)]
    hmmOut<- data.table(state= NA, postM= NA)[0,]
    outHeader<- TRUE
    fromID<- 1
    for (xchr in unique(bed$chrom)) {
        write(sprintf('Chrom: %s', xchr), stderr())
        xbed<- bed[chrom == xchr, ]
        hmmChrom<- HMMChrom(xbed)
        ## At the end of each chromn write to results
        write.table(hmmChrom, stdout(), sep= '\t', col.names= outHeader, row.names= FALSE, quote= FALSE)
        outHeader<- FALSE    
        hmmOut<- rbindlist(list(hmmOut, hmmChrom[, list(state, postM)]))
    }
    return(hmmOut)
}

# END_OF_FUNCTIONS <- Don't change this string it is used to source in run_test.R
# ==============================================================================
# If you just want to use these functions in via source(...) exit after having sourced.
if(interactive()){
    options(show.error.messages=FALSE)
    on.exit(options(show.error.messages=TRUE))
    stop()
}

# ==============================================================================

done<- suppressWarnings(suppressMessages(require(argparse)))
if(done == FALSE){
    cat('\nPlease install the "argparse" package. Open an R session and execute:\n\n')
    cat('> install.packages("argparse")\n\n')
    quit(save= 'no', status= 1)
}


docstring<- sprintf("DESCRIPTION \\n\\
Segment pvalues in input file in two states M: stretches of low pvalues; u: otherwise.\\n\\
Input file must be tab separated with a column of chromosomes and a column of p-values. \\n\\
\\n\\
OUTPUT:\\n\\
Tab separated file to stdout with columns: \\n\\
chrom, input pvalue, recoded pvalue, state, posterior prob of state M, state_id.\\n\\
\\n\\
Version %s", VERSION)
parser<- ArgumentParser(description= docstring, formatter_class= 'argparse.RawTextHelpFormatter')
parser$add_argument("-i", "--input", help= "Input file, can be gzip'd, use - to read from stdin. Must have chrom names in column 1.", required= TRUE)
parser$add_argument("-p", "--pIndex", help= "Column index with raw pvalues. 1-based.", required= TRUE, type= 'integer')
parser$add_argument("-c", "--chromIndex", help= "Column index with chromosome. Default 1", default= 1, type= 'integer')
parser$add_argument("-H", "--header", help= "Input file has header", action= 'store_true')
parser$add_argument("-P", "--recodePval", help= "List of three floats to recode pvalues in discrete categories. Default: 0.1 0.05 0.001.\\n\\
I.e. pvals recoded as 0 if >= 0.1; 1 if 0.1 < p <= 0.5; 2 if 0.05 < p <= 0.001; 3 if p < 0.001.\\n\\
To recode to less then 3 categories terminate the list with 0 e.g. `0.05 0 0`\\n\\
to classify pvalues in >=0.05 or < 0.05", nargs= 3, type= 'double', default= c(0.1, 0.05, 0.001))

xargs<- parser$parse_args()

# ==============================================================================

bed<- readBedFile(xargs$input, xargs$pIndex, chromIndex= xargs$chromIndex, header= xargs$header)
states<- HMMrunner(bed, recodePval= xargs$recodePval)
quit(save= 'no')
