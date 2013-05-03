## source('~/svn_checkout/bioinformatics-misc/BSdata.R')
read.bsdata<- function(x, ...){
    "Read the concatenated output of mpileup2methylation.py to a list similar to
    limma data structure.
    x:
        File to read. Header line absent (use header=TRUE otherwise).
        Columns are expected to be:
        chrom	start	end	pct.met	cnt.met	tot_reads	strand	library_id
    ...:
        Further arguments passed to read.table e.g. header
    "
    BSlist<- list()
    bsdata<- read.table(x, sep= '\t', colClasses= c('character', 'integer', 'integer', 'numeric', 'integer', 'integer', 'character', 'character'), ...)
    names(bsdata)[1:8]<- c('chrom', 'start', 'end', 'pct.met', 'cnt.met', 'tot_reads', 'strand', 'library_id')
    bsdata$locus<- paste(bsdata$chrom, bsdata$start, bsdata$end, sep= '_')
    
    ## Get union of all positions. Make it BED format
    allPos<- unique(bsdata[, c('chrom', 'start', 'end', 'locus', 'strand')])
    allPos$score<- '.'
    BSlist$loci<- allPos[, c('chrom', 'start', 'end', 'locus', 'score', 'strand')]
    rm(allPos)
    
    ## Matrix of total counts
    tot_reads<- reshape(bsdata[, c('locus', 'library_id', 'tot_reads')], timevar= 'library_id', v.names= 'tot_reads', idvar= 'locus', direction= 'wide')
    rownames(tot_reads)<- tot_reads$locus
    tot_reads<- as.matrix(tot_reads[, 2:ncol(tot_reads)])
    colnames(tot_reads)<- sub('^tot_reads\\.', '', colnames(tot_reads), perl= TRUE)
    BSlist$tot_reads<- tot_reads
    BSlist$tot_reads[is.na(BSlist$tot_reads)]<- 0
    rm(tot_reads)
    
    ## Matrix of methylated counts
    cnt.met<- reshape(bsdata[, c('locus', 'library_id', 'cnt.met')], timevar= 'library_id', v.names= 'cnt.met', idvar= 'locus', direction= 'wide')
    rownames(cnt.met)<- cnt.met$locus
    cnt.met<- as.matrix(cnt.met[, 2:ncol(cnt.met)])
    colnames(cnt.met)<- sub('^cnt.met\\.', '', colnames(cnt.met), perl= TRUE)
    BSlist$cnt_met<- cnt.met
    BSlist$cnt_met[is.na(BSlist$cnt_met)]<- 0
    rm(cnt.met)
    if (!all(rownames(BSlist$cnt_met) == rownames(BSlist$tot_reads))){
        stop('Unexpected row names')
    }
    if (!all(BSlist$locus == rownames(BSlist$tot_reads))){
        stop('Unexpected locus names or orders')
    }
    ## Matrix of percentage methylated
    BSlist$pct_met<- 100*(BSlist$cnt_met / BSlist$tot_reads)
    return(BSlist)
}

