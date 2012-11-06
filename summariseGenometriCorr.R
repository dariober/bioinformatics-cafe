summariseGenometriCorr<- function(x, xslot= "awhole"){
    # Summarize the statistics in the output object of GenometriCorrelation()
    #
    #                       Arguments
    # x:
    #   Output of GenometriCorrelation
    # awhole:
    #   Slot in x to summarize (e.g. 'chr1' or 'awhole')
    #
    #                       Values
    # Statistics as dataframe
    
   exlude_stats<- c("relative.distances.data",
                    "projection.test",
                    "absolute.min.distance.data",
                    "absolute.inter.reference.distance.data",
                    "relative.distances.ecdf.deviation.area.null.list",
                    "scaled.absolute.min.distance.sum.null.list",
                    "jaccard.measure.null.list")
    
    dataset= x[[xslot]]
    for(ptest in dataset["projection.test"]){
        dataset<- c(dataset, ptest)
    }
    stats<- names(dataset)[!names(dataset) %in% exlude_stats]
    stats_list<- as.data.frame(t(as.matrix(unlist(dataset[stats]))))
    return(stats_list)
}
