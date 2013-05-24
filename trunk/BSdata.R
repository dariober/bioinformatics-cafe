## source('~/svn_checkout/bioinformatics-misc/BSdata.R')

### Objject to store BS data in format similar to limma

BSdata<- setClass("BSdata", representation(
                                        loci= 'data.frame',
                                        tot_reads= 'matrix',
                                        cnt_met= 'matrix',
                                        pct_met= 'matrix',
                                        design= 'data.frame'))
                                        
read.bsdata<- function(x, ...){
    "Read the concatenated output of mpileup2methylation.py to a BSdata object
    x:
        File to read. Header line absent (use header=TRUE otherwise).
        Columns are expected to be:
        chrom	start	end	pct.met	cnt.met	tot_reads	strand	library_id
    ...:
        Further arguments passed to read.table e.g. header
    "
    BSobj<- BSdata()
    bsdata<- read.table(x, sep= '\t', colClasses= c('character', 'integer', 'integer', 'numeric', 'integer', 'integer', 'character', 'character'), ...)
    names(bsdata)[1:8]<- c('chrom', 'start', 'end', 'pct.met', 'cnt.met', 'tot_reads', 'strand', 'library_id')
    bsdata$locus<- paste(bsdata$chrom, bsdata$start, bsdata$end, sep= '_')
    
    ## Get union of all positions. Make it BED format
    allPos<- unique(bsdata[, c('chrom', 'start', 'end', 'locus', 'strand')])
    allPos$score<- '.'
    BSobj@loci<- allPos[, c('chrom', 'start', 'end', 'locus', 'score', 'strand')] ## order(allPos$chrom, allPos$start, allPos$end)
    rm(allPos)
    
    ## Matrix of total counts
    tot_reads<- reshape(bsdata[, c('locus', 'library_id', 'tot_reads')], timevar= 'library_id', v.names= 'tot_reads', idvar= 'locus', direction= 'wide')
    rownames(tot_reads)<- tot_reads$locus
    tot_reads<- as.matrix(tot_reads[, 2:ncol(tot_reads)])
    colnames(tot_reads)<- sub('^tot_reads\\.', '', colnames(tot_reads), perl= TRUE)
    BSobj@tot_reads<- tot_reads
    BSobj@tot_reads[is.na(BSobj@tot_reads)]<- 0
    rm(tot_reads)
    
    ## Matrix of methylated counts
    cnt.met<- reshape(bsdata[, c('locus', 'library_id', 'cnt.met')], timevar= 'library_id', v.names= 'cnt.met', idvar= 'locus', direction= 'wide')
    rownames(cnt.met)<- cnt.met$locus
    cnt.met<- as.matrix(cnt.met[, 2:ncol(cnt.met)])
    colnames(cnt.met)<- sub('^cnt.met\\.', '', colnames(cnt.met), perl= TRUE)
    BSobj@cnt_met<- cnt.met
    BSobj@cnt_met[is.na(BSobj@cnt_met)]<- 0
    rm(cnt.met)
    if (!all(rownames(BSobj@cnt_met) == rownames(BSobj@tot_reads))){
        stop('Unexpected row names')
    }
    if (!all(BSobj@loci$locus == rownames(BSobj@tot_reads))){
        stop('Unexpected locus names or orders')
    }
    ## Matrix of percentage methylated
    BSobj@pct_met<- 100*(BSobj@cnt_met / BSobj@tot_reads)
    
    ## Library names extracted from column headers go to design
    BSobj@design<- data.frame(index= 1:ncol(BSobj@tot_reads), library_id= colnames(BSobj@tot_reads), bs= rep(NA, ncol(BSobj@tot_reads)))
    
    return(BSobj)
}

makeBSdataLong<- function(bsobj){
    "Make a long format of bsdata obj"
    if(class(bsobj) != "BSdata"){
        stop('Object os not of class BSdata')
    }
    longf<- data.frame(locus= rep(rownames(bsobj@cnt_met), ncol(bsobj@cnt_met)),
        cnt_M= c(bsobj@cnt_met),
        cnt_m= c(bsobj@tot_reads - bsobj@cnt_met),
        tot_reads= c(bsobj@tot_reads),
        library_id= rep(colnames(bsobj@cnt_met), each= nrow(bsobj@cnt_met))
    )
    longf<- merge(longf, bsobj@design[, c('library_id', 'bs')], by.x= c('library_id'), by.y= c('library_id'), sort= FALSE)
    longf<- longf[order(longf$locus, longf$library_id),]
    return(longf)
}

checkMatrices<- function(bsobj){
    "Return TRUE if a number of checks on BSdata matrices are passed"
    censor<- TRUE
    for(s in slotNames(bsobj)){
        m<- slot(bsobj, s)
        if(class(m) != 'matrix'){
            next
        }
        if(min(m, na.rm= TRUE) < 0){
            stop(sprintf('Value < 0 found in slot "%s"', s))
        }
        if(censor){
            rnames<- rownames(m)
            cnames<- colnames(m)
            censor<- FALSE
        } else {
            if(all(rnames == rownames(m)) != TRUE){
                stop("Non-matching row names")
            }
            if(all(cnames == colnames(m)) != TRUE){
                stop("Non-matching column names")
            }
        }
    }
    return(TRUE)
}

BSdataApply<- function(bsobj, slots= c('tot_reads', 'cnt_met', 'pct_met'), FUN, ...){
    "Apply the function FUN to the *matrices* in BSdata object.
    Returns a BSdata object with the matrices apply'd
    FUN:
        A function returning a matrix
    ...:
        Further args passed to FUN
    ___________________________________________________________________________
    Example:
        BSdataApply(bsdata, FUN= function(x) x[1:10,])
    "
    if(all(slots %in% slotNames(bsobj)) == FALSE){
       stop('Some slot names not found in slotNames(bsobj)') 
    }
    bsapp<- bsobj
    for(s in slots){
        m<- slot(bsobj, s)
        if(class(m) != 'matrix'){
            stop(sprintf('Slot "%s" is not a matrix', s))
        }
        slot(bsapp, s)<- FUN(m, ...)
    }
    checkMatrices(bsapp)
    return(bsapp)
}

BSmatrix2begraph<- function(bsobj, m, outdir= '.', header= FALSE){
    "Produce a bedGraph file for each column (library) of a matrix in a BSdata object
    bsobj:
        A BSdata object    
    m:
        Slot name of the matrix to convert to bedGraph
    outdir:
        Output directory for the files. Default to cwdir.
        Files named as <column-name>.<matrix-name>.bedGraph
    Example:
        BSmatrix2begraph(bsobj, 'pct_met')
    Assumptions:
        The row name of the matrix has match in bsobj@loci$locus from which chrom
        coordinates can be extracted
    "
    bdg<- slot(bsobj, m)
    if(all(rownames(bdg) %in% bsobj@loci$locus) == FALSE){
        stop('Not all rownames in the bsobj could be extracted from bsobj@loci$locus')
    }
    libraries<- colnames(bdg)
    locus<- rownames(bdg)
    bdg<- as.data.frame(bdg)
    bdg$locus<- locus
    bdg<- merge(bdg, bsobj@loci[, c('locus', 'chrom', 'start', 'end')], by.x= c('locus'), by.y= c('locus'))
    bdg<- bdg[order(bdg$chrom, bdg$start, bdg$end),]
    dir.create(outdir, showWarnings = FALSE, recursive = FALSE)
    for(x in libraries){
        outfile<- paste(x, m, 'bedGraph', sep= '.')
        write.table(bdg[, c('chrom', 'start', 'end', x)], file.path(outdir, outfile), row.names= FALSE, col.names= header, sep= '\t', quote= FALSE)
    }
}



