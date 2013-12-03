readMarkDuplicates<- function(x){
    # Read output of picard MarkDuplicates
    # x:
    #   Metrics file from markDuplicates
    # Return:
    #   List with compenents metrics and hist. 
    # --------------------------------------------------------------------------
    mdup<- readLines(x)
    
    ## Read summary metrics:
    headerIdx<- grep('## METRICS CLASS\tnet.sf.picard.sam.DuplicationMetrics', mdup) + 1
    met<- mdup[headerIdx:(headerIdx + 1)]
    header<- strsplit(met[1], '\t')[[1]]
    metrics<- as.data.frame(t(strsplit(met[2], '\t')[[1]]), stringsAsFactors= FALSE)
    colnames(metrics)<- header
    for(i in 2:ncol(metrics)){
        metrics[, i]<- as.numeric(metrics[, i])
    }
    
    ## Histogram
    header<- grep('BIN\tVALUE', mdup)
    if(length(header) > 0){
        histDup<- mdup[header:length(mdup)]
        histDup<- t(sapply(strsplit(histDup, '\t'), function(x) x[1:2]))
        colnames(histDup)<- histDup[1,]
        histDup<- as.data.frame(histDup[2:nrow(histDup),], stringsAsFactors= FALSE, colClasses= c('numeric', 'numeric'))
        histDup<- histDup[complete.cases(histDup),]
    } else {
        histDup<- NULL
    }
    ##
    mdupList<- list(metrics= metrics, hist= histDup)
    return(mdupList)
}

