readMarkDuplicates<- function(x){
    # Read output of picard MarkDuplicates
    # x:
    #   Metrics file from markDuplicates
    # Return:
    #   List with compenents metrics and hist. 
    # --------------------------------------------------------------------------
    for (i in 1:length(x)){
        mdup<- readLines(x[i])
        
        ## Read summary metrics:
        headerIdx<- grep('## METRICS CLASS\tnet.sf.picard.sam.DuplicationMetrics', mdup) + 1
        met<- mdup[headerIdx:(headerIdx + 1)]
        header<- strsplit(met[1], '\t')[[1]]
        metrics<- as.data.frame(t(strsplit(met[2], '\t')[[1]]), stringsAsFactors= FALSE)
        colnames(metrics)<- header
        for(j in 2:ncol(metrics)){
            metrics[, j]<- as.numeric(metrics[, j])
        }
        metrics$filename<- x[i]
        ## Histogram
        header<- grep('BIN\tVALUE', mdup)
        if(length(header) > 0){
            histDup<- mdup[header:length(mdup)]
            histDup<- t(sapply(strsplit(histDup, '\t'), function(y) y[1:2]))
            colnames(histDup)<- histDup[1,]
            histDup<- as.data.frame(histDup[2:nrow(histDup),], stringsAsFactors= FALSE)
            histDup<- histDup[complete.cases(histDup),]
            histDup$filename<- x[i]
            histDup$BIN<- as.numeric(histDup$BIN)
            histDup$VALUE<- as.numeric(histDup$VALUE)
        } else {
            histDup<- NULL
        }
        if(i == 1){
            metricsCat<- metrics
            histDupCat<- histDup
        } else {
            metricsCat<- rbind(metricsCat, metrics)
            histDupCat<- rbind(histDupCat, histDup)
        }
    }
    ##
    mdupList<- list(metrics= metricsCat, hist= histDupCat)
    return(mdupList)
}

