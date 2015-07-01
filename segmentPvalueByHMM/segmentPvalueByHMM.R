#! /usr/bin/env Rscript

## TODO:
## * If a chrom is fully demethylated or methylated the state might be wrongly assigned!!
## * Print results as they are produced by each loop of HMMrunner()

VERSION<- '0.1.0'
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

## Recode pvalues to be discrete
recodeOnePvalue<- function(p){
    if (is.na(p)){
        o<- NA
    } else if (p >= 0.1 ){
        o<- 0
    } else if (p < 0.1 & p >= 0.05) {
        o<- 1
    } else if (p < 0.05 & p >= 0.001) {
        o<- 2
    } else if (p< 0.001){
        o<- 3
    } else {
        stop('Unexpected p-value')
    }
    return(o)
}
recodePvals<- function(p){
    stopifnot(is.numeric(p))
    stopifnot(min(p, na.rm= TRUE) >= 0)
    stopifnot(max(p, na.rm= TRUE) <= 1)
    recoded<- sapply(p, recodeOnePvalue)
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
    #   renaming
    # Return:
    #    Vector of renamed states
    # Example:
    #   states<- c(2,2,2,1,1,1,2,2)
    #   obs<- c(0.1, 0.15, 0.2, 0.9, 0.8, 0.85, 0.1, 0.2)
    #   levels<- c('L', 'H')
    #   recodeStates(states, obs, levels) #<<< c(L, L, L, H, H, H, L, L)
    # -------------------------------------------------------------------------
    stopifnot(length(states) == length(obs))
    avg<- data.table(states, obs)
    avg<- avg[, list(stateavg= mean(obs, na.rm= TRUE)), by= states]
    avg<- avg[order(stateavg), ]
    avg$newStates<- levels
    newx<- levels[match(states, avg$states)]
    return(newx)
}

HMMrunner<- function(bed, MAXPOS= 40000){
    # Run HMM. bed is a bed file with columns 'chrom' and 'pvals'
    bed[, recode_P := recodePvals(pvals)]
    states<- vector()
    posteriors<- vector()
    hmmOut<- data.table(state= NA, postM= NA)[0,]
    for (xchr in unique(bed$chrom)) {
        ## Loop through each chromosome to fit HMM.
        ## Need to split chroms in chunks due to memory limit
        xbed<- bed[chrom == xchr, ]
        startsEnds<- prepareChromChunks(nrow(xbed), MAXPOS)
        starts<- startsEnds$starts
        ends<- startsEnds$ends
    
        hmmNeeded<- TRUE
        if(length(xbed$recode_P) < 500000){
            ## Fit the model to entire chrom, if not too big:
            hmmfit<- HMMFit(xbed$recode_P, nStates= 2, dis= 'DISCRETE')
            hmmNeeded<- FALSE
        }
        # print(hmmfit$HMM$transMat)
        for (j in 1:length(starts)){
            write(sprintf('%s %s %s %s', xchr, starts[j], ends[j], ends[j] - starts[j]), stderr())
            dat<- xbed[starts[j]:ends[j],]
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
        hmmOutSub<- data.table(state= xstates, posteriors= prob[,mx])
        
        ## Write to stout as results come through
        write.table(cbind(xbed, hmmOutSub), stdout(), sep= '\t', col.names= TRUE, row.names= FALSE, quote= FALSE)

        hmmOut<- rbindlist(list(hmmOut, hmmOutSub))
        }
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
chrom, input pvalue, recoded pvalue, state, posterior prob of state M.\\n\\
\\n\\
Version %s", VERSION)
parser<- ArgumentParser(description= docstring, formatter_class= 'argparse.RawTextHelpFormatter')
parser$add_argument("-i", "--input", help= "Input file, can be gzip'd, use - to read from stdin. Must have chrom names in column 1.", required= TRUE)
parser$add_argument("-p", "--pIndex", help= "Column index with raw pvalues. 1-based.", required= TRUE, type= 'integer')
parser$add_argument("-c", "--chromIndex", help= "Column index with chromosome. Default 1", default= 1, type= 'integer')
parser$add_argument("-H", "--header", help= "Input file has header", action= 'store_true')

xargs<- parser$parse_args()

# ==============================================================================

bed<- readBedFile(xargs$input, xargs$pIndex, chromIndex= xargs$chromIndex, header= xargs$header)
states<- HMMrunner(bed)
quit(save= 'no')
