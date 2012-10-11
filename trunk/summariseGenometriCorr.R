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
    
    stats<-   c("query.population",
                "reference.population",
            #    "relative.distances.data",
                "relative.distances.ks.p.value",
                "relative.distances.ecdf.deviation.area",
                "relative.distances.ecdf.area.correlation",
            #    "projection.test",
            #    "query.reference.intersection",
            #    "query.reference.union",
            #    "jaccard.measure",
                "projection.test.p.value",
                "projection.test.lower.tail",
            #    "absolute.min.distance.data",
            #    "absolute.inter.reference.distance.data",
            #    "scaled.absolute.min.distance.sum",
            #    "relative.distances.ecdf.deviation.area.null.list",
            #    "scaled.absolute.min.distance.sum.null.list",
            #    "jaccard.measure.null.list",
            #    "relative.distances.ecdf.deviation.area.p.value",
                "scaled.absolute.min.distance.sum.p.value",
                "scaled.absolute.min.distance.sum.lower.tail",
                "jaccard.measure.p.value",
                "jaccard.measure.lower.tail")
    
    dataset= x[[xslot]]
    stats_list<- as.data.frame(t(as.matrix(unlist(dataset[stats]))))
    return(stats_list)
}