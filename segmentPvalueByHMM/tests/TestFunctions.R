test.canReadBedFile<- function(){
    pIndex<- 4
    bed<- readBedFile('tests/input-1.bed', pIndex)
    checkTrue(all(bed[, pvals] == c(0.1, 0.2, 0.3)))

    bed<- readBedFile('tests/input-1.bed', pIndex, header= TRUE)
    checkTrue(all(bed[, pvals] == c(0.2, 0.3)))
}

test.canRecodePvaluesToDiscrete<- function(){
    pvals<- c(0, 0.9, 0.07, 0.01, NA, 0.7, 0.3)
    expt<-  c(3, 0, 1, 2, NA, 0, 0)
    rpval<- recodePvals(pvals, recodePval= c(0.1, 0.05, 0.001))
    checkTrue(all(expt == rpval, na.rm= TRUE))

    rpval<- recodePvals(pvals, recodePval= c(0.3, 0, 0)) # Two categories: >= 0.3 or < 0.3
    expt<-  c(1, 0, 1, 1, NA, 0, 0)
    checkTrue(all(expt == rpval, na.rm= TRUE))
    
    checkException(recodePvals(-1))
    checkException(recodePvals(1.1))
    checkException(recodePvals(c(0, 'Bla')))
    checkException(recodePvals(pvals, recodePval= c(0.1, 0.05, -1))) # -1 not allowed
}

test.canDivideChromInChunks<- function(){
    chunks<- prepareChromChunks(100, 100)
    checkEquals(chunks$starts, 1)
    checkEquals(chunks$end, 100)
    
    chunks<- prepareChromChunks(100, 10)
    expt<- seq(0, 90, by= 10)+1
    checkEquals(expt, chunks$starts)
    expt<- seq(10, 100, by= 10)
    checkEquals(expt, chunks$end)
}

test.canRecodeStates<- function(){
    expt<- c('L', 'L', 'L', 'H', 'H', 'H', 'L', 'L')
    states<- c(2,2,2,1,1,1,2,2)
    obs<- c(0.1, 0.15, 0.2, 0.9, 0.8, 0.85, 0.1, 0.2)
    levels<- c('L', 'H')
    recoded<- recodeStates(states, obs, levels)
    checkEquals( expt, recoded )
    
    # Swith states, recoding should stay the same
    states<- c(1,1,1,2,2,2,1,1)
    recoded<- recodeStates(states, obs, levels)
    checkEquals( expt, recoded )
}

test.canRecodeOneState<- function(){
    # A stretch of only high pvals.
    expt<- rep('u', 10)
    states<- rep(1, 10)
    obs<- runif(n= 10, min= 0.8, max= 1)

    recoded<- recodeStates(states, obs, c('M', 'u'))
    checkEquals( expt, recoded )
    
    ## Only small pvals
    expt<- rep('M', 10)
    states<- rep(1, 10)
    obs<- runif(n= 10, min= 0, max= 0.1)

    recoded<- recodeStates(states, obs, c('M', 'u'))
    checkEquals(expt, recoded)
}

test.canRunHmmAndGetStates<- function(){
    bed<- data.table(chrom= rep('chr1', 100), pvals= c(
            runif(n= 30, min= 0.2, max= 1),
            runif(n= 10, min= 0, max= 0.01),
            runif(n= 60, min= 0.2, max= 1)))
    states<- HMMrunner(bed)
    expt<- c(rep('u', 30), rep('M', 10), rep('u', 60))
    checkEquals( expt, states$state )
    
    bed2<- bed[, chrom := 'chr2']
    bed3<- bed[, chrom := 'chr3']
    bbed<- rbind(bed, bed2, bed3)
    states<- HMMrunner(bbed)
    expt2<- c(expt, expt, expt)
    checkEquals( expt2, states$state )

    # Put low pvals first
    bed<- data.table(chrom= rep('chr1', 100), pvals= c(
            runif(n= 10, min= 0, max= 0.01),
            runif(n= 80, min= 0.2, max= 1),
            runif(n= 10, min= 0, max= 0.05)))
    states<- HMMrunner(bed)
    expt<- c(rep('M', 10), rep('u', 80), rep('M', 10))
    checkEquals( expt, states$state )
}

test.canRunHmmAndGetPosteriors<- function(){
    bed<- data.table(chrom= rep('chr1', 100), pvals= c(
            runif(n= 30, min= 0.2, max= 1),
            runif(n= 10, min= 0, max= 0.01),
            runif(n= 60, min= 0.2, max= 1)))
    states<- HMMrunner(bed)

    post<- states[, list(avg_p= mean(postM)), by= state]
    checkTrue(post[state == 'M', avg_p] > post[state == 'u', avg_p])
    checkTrue(all(states[state == 'M', postM] > 0.8))
    checkTrue(all(states[state == 'u', postM] < 0.5))

    # Put low pvals first
    bed<- data.table(chrom= rep('chr1', 100), pvals= c(
            runif(n= 10, min= 0, max= 0.01),
            runif(n= 80, min= 0.2, max= 1),
            runif(n= 10, min= 0, max= 0.05)))
    states<- HMMrunner(bed)
    post<- states[, list(avg_p= mean(postM)), by= state]
    checkTrue(post[state == 'M', avg_p] > post[state == 'u', avg_p])
    checkTrue(all(states[state == 'M', postM] > 0.8))
    checkTrue(all(states[state == 'u', postM] < 0.5))
}

test.canRunScript<- function(){
    cmd<- './segmentPvalueByHMM.R -i tests/input-2.txt.gz -p 2'
    exitCode<- system(cmd)
    checkEquals(0, exitCode)
    
    ## Can read from stdin:
    cmd<- 'gunzip -c tests/input-2.txt.gz  | ./segmentPvalueByHMM.R -i - -p 2'
    exitCode<- system(cmd)
    checkEquals(0, exitCode)

    ## Multiple chroms:
    cmd<- './segmentPvalueByHMM.R -i tests/input-3.txt -p 2'
    exitCode<- system(cmd)
    checkEquals(0, exitCode)

    ## Custom pvalue thresholds:
    cmd<- './segmentPvalueByHMM.R -i tests/input-2.txt.gz -p 2 -P 0.01 0 0'
    exitCode<- system(cmd)
    checkEquals(0, exitCode)
}

