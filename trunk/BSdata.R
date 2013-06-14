## source('~/svn_checkout/bioinformatics-misc/BSdata.R')

library(ff)
library(ffbase)

### Object to store BS data in format similar to limma

setOldClass("ff_matrix")

BSdata<- setClass("BSdata", representation(
                                        loci= 'data.frame',
                                        tot_reads= 'ff_matrix',
                                        cnt_met= 'ff_matrix',
                                        pct_met= 'ff_matrix',
                                        design= 'data.frame'))
                                        
read.bsdata<- function(infiles, library_ids= FALSE, ...){
    # ---------------------------------------------------
    # Read the concatenated output of mpileup2methylation.py to a BSdata object
    # infiles:
    #     Vector of files to read. Header line absent (use header=TRUE otherwise).
    #     Columns are expected to be:
    #     chrom	start	end	pct.met	cnt.met	tot_reads	strand	[library_id]
    # library_ids:
    #     Vector of (unique) library IDs assigned to each file in x
    # save.name:
    #     Name of file where the image will be saved by calling ffsave.image().
    # ...:
    #     Further arguments passed to read.table e.g. header
    #
    # MEMO: use Sys.glob('*.bedGraph.gz') to create a vectior of files with given glob pattern.
    # ---------------------------------------------------
  
    BSobj<- BSdata()
    
    if(length(library_ids)!=length(unique(library_ids))){
        stop(sprintf('Duplicate library_ids found in %s', library_ids))
    }    
    if(length(infiles)!=length(unique(infiles))){
        stop(sprintf('Duplicate input files found in %s', infiles))
    } 
    if(length(infiles)!=length(library_ids)){
        stop('Number of input files does not equal number of library IDs')
    }
    colClasses<- c('factor', 'integer', 'integer', 'numeric', 'integer', 'integer', 'factor')
    colNames<-  c('chrom', 'start', 'end', 'pct.met', 'cnt.met', 'tot_reads', 'strand')
    ftempl= FALSE
    for(i in seq(1, length(infiles))){
        xfile<- infiles[i]
        libid<- library_ids[i]
	cat(sprintf('Reading file: %s; ID: %s\n', xfile, libid))
	if(!ftempl){
 	    ## Read in the first file as "template"
	    bsdata<- read.table.ffdf(x= NULL, file= xfile, sep= '\t', col.names= colNames, colClasses= colClasses, ...)
            bsdata$library_id<- as.ff(as.factor(rep(libid, nrow(bsdata))), vmode= 'short')
            ftempl<- TRUE
	} else {
	    newdata<- read.table.ffdf(x= NULL, file= xfile, sep= '\t', col.names= colNames, colClasses= colClasses, ...)
            newdata$library_id<- as.ff(as.factor(rep(libid, nrow(newdata))), vmode= 'short')
            #bsdata<- read.table.ffdf(x= bsdata, file= xfile, sep= '\t', colClasses= colClasses, col.names= colNames, ...)
            bsdata<- ffdfappend(bsdata, newdata, recode= FALSE, adjustvmode= FALSE)
	}
	# bsdata<- read.table(x, sep= '\t', colClasses= c('character', 'integer', 'integer', 'numeric', 'integer', 'integer', 'character', 'character'), ...)	
    }
    #rm(newdata)
    # bsdata<- read.table(x, sep= '\t', colClasses= c('character', 'integer', 'integer', 'numeric', 'integer', 'integer', 'character', 'character'), ...)
    # names(bsdata)[1:8]<- c('chrom', 'start', 'end', 'pct.met', 'cnt.met', 'tot_reads', 'strand', 'library_id')
    cat('Generating union of all positions...\n')
    bsdata$locus<- as.ff(as.factor( paste(bsdata$chrom[1:nrow(bsdata)], bsdata$start[1:nrow(bsdata)], bsdata$end[1:nrow(bsdata)], sep= '_')), vmode= 'integer') ## paste(bsdata$chrom, bsdata$start, bsdata$end, sep= '_')
    ## Get union of all positions. Make it BED fo#rmat
    allPos<- unique(bsdata[, c('chrom', 'start', 'end', 'locus', 'strand')])
    allPos$score<- '.'
    BSobj@loci<- allPos[, c('chrom', 'start', 'end', 'locus', 'score', 'strand')] ## order(allPos$chrom, allPos$start, allPos$end)
    #rm(allPos)
    ## Matrix of total counts
    cat('Creating matrix of total counts...\n')
    tot_reads<- reshape(bsdata[, c('locus', 'library_id', 'tot_reads')], timevar= 'library_id', v.names= 'tot_reads', idvar= 'locus', direction= 'wide')
    rownames(tot_reads)<- tot_reads$locus
    tot_reads<- as.matrix(tot_reads[, 2:ncol(tot_reads)], storage.mode= 'integer')
    colnames(tot_reads)<- sub('^tot_reads\\.', '', colnames(tot_reads), perl= TRUE)
    tot_reads[is.na(tot_reads)]<- as.integer(0)
    tot_reads<- as.ff(tot_reads)
	BSobj@tot_reads<- tot_reads
    nr<- nrow(BSobj@tot_reads)
    nc<- ncol(BSobj@tot_reads)
    #rm(tot_reads)
       
    
    ## Matrix of methylated counts
    cat('Creating matrix of methylated counts...\n')
    cnt.met<- reshape(bsdata[, c('locus', 'library_id', 'cnt.met')], timevar= 'library_id', v.names= 'cnt.met', idvar= 'locus', direction= 'wide')
    rownames(cnt.met)<- cnt.met$locus
    cnt.met<- as.matrix(cnt.met[, 2:ncol(cnt.met)], storage.mode= 'integer')
    colnames(cnt.met)<- sub('^cnt.met\\.', '', colnames(cnt.met), perl= TRUE)
    cnt.met[is.na(cnt.met)]<- as.integer(0)
    cnt.met<- as.ff(cnt.met)
	BSobj@cnt_met<- cnt.met
    #rm(cnt.met)
    if (!all(rownames(BSobj@cnt_met) == rownames(BSobj@tot_reads))){
        stop('Unexpected row names')
    }
    if (!all(BSobj@loci$locus == rownames(BSobj@tot_reads))){
        stop('Unexpected locus names or orders')
    }
    
    ## Matrix of percentage methylated
    pct_met<- as.ff(100*(BSobj@cnt_met[1:nr, 1:nc] / BSobj@tot_reads[1:nr, 1:nc]))
	BSobj@pct_met<- pct_met
    
    ## Library names extracted from column headers go to design
    BSobj@design<- data.frame(index= 1:ncol(BSobj@tot_reads), library_id= colnames(BSobj@tot_reads), bs= rep(NA, ncol(BSobj@tot_reads)))

    return(BSobj)
}

save.BSdata<- function(BSobj, outname){
	"Save BSdata object
	BSobj:
		Object to save
	outname:
		Output name passed to ffsave. Do not add '.Rdata' (ffsave will do it)
	
	To reload: ffload(BSobj)
	"
	bsobj<- bsobj<- deparse(substitute(BSobj))   ## Return string with name of BSobj. eg if (BSobj= cpg, ...) return string "cpg"
	assign(paste(bsobj, 'tot_reads', sep= '.'), BSobj@tot_reads)  ## Create a variable fpr each ff obj e.g cpg.tot_reads
	assign(paste(bsobj, 'cnt_met', sep= '.'), BSobj@cnt_met)
	assign(paste(bsobj, 'pct_met', sep= '.'), BSobj@pct_met)
	ffsave(list= c(bsobj, paste(bsobj, 'tot_reads', sep= '.'), paste(bsobj, 'cnt_met', sep= '.'), paste(bsobj, 'pct_met', sep= '.')), file= outname)
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
        if(class(m)[1] != 'matrix'){
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
        if(class(m)[1] != 'ff_matrix'){
            stop(sprintf('Slot "%s" is not a matrix', s))
        }
        slot(bsapp, s)<- as.ff(FUN(m, ...))
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



