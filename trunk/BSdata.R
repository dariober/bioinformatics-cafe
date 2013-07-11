## source('~/svn_checkout/bioinformatics-misc/BSdata.R')

library(ff)
library(ffbase)
# library(bigmemory)
# library(foreach)
# library(doMC)

setOldClass("ff_matrix")
setOldClass("ffdf")


name.reader<- function(m){
	"Get colnames from input file (matrix), skipping first column (this is row names)"
	dfsample<- read.table(m, sep= '\t', nrows= 1, colClasses= 'character')
	colNames<- as.character(dfsample[1,2:ncol(dfsample)])
	return(colNames)
}


BSdata<- setClass("BSdata", representation(
## Object to store BS data in format similar to limma
                                        loci= 'ffdf',
                                        tot_reads= 'ffdf',
                                        cnt_met= 'ffdf',
                                        pct_met= 'ffdf',
                                        design= 'data.frame'))

read.bsdata<- function(prefix, gzip= TRUE, mat= c('loci', 'pct_met', 'cnt_met', 'tot_reads'), ...){
	"Read the output files from BSmatrix.py and put them in a BSdata object

	ARGS:
	prefix:
		Prefix of the input files produced by BSmatrix.py. It is expected to find
			<prefix>.loci.gz
			<prefix>.cnt_met.mat.gz
			<prefix>.pct_met.mat.gz
			<prefix>.tot_reads.mat.gz	
		Input files can be either gzipped or not.
	gzip:
		Logical. Are the input files gzipped? (I.e. ending in .gz). Default TRUE
	mat:
		Which matrices should be read? Default is all of them:
		c('loci', 'pct_met', 'cnt_met', 'tot_reads').
		For speed of reading consider leaving out pct_met
	...:
		Further arguments passed to read.table, e.g. nrows.
	VALUE:
		BSdata object
	"
    BSobj<- BSdata()

	if( gzip ){
		gz= '.gz'
	} else {
		gz= ''
	}
	
	## File suffixes
	sfx<- list(
		loci= paste('.loci', gz, sep= ''),
		pct_met= paste('.pct_met.mat', gz, sep= ''),
		cnt_met= paste('.cnt_met.mat', gz, sep= ''),
		tot_reads= paste('.tot_reads.mat', gz, sep= '')
	)
	inf<- sfx
	for(i in 1:length(inf)){
		"Input files"
		inf[[i]]<- paste(prefix, inf[[i]], sep= '')
	}
	for(f in inf){
		if( !file.exists(f) ){
			stop(sprintf('File "%s" not found', f))
		}
	}
	dfsample<- read.table(inf$cnt_met, sep= '\t', nrows= 1)

	if ('loci' %in% mat){
		cat(sprintf('Reading file "%s"...\n', inf$loci))
		BSobj@loci<- read.table.ffdf(file= inf$loci, header= TRUE, ...)
	}
	if ('tot_reads' %in% mat){
		cat(sprintf('Reading total file "%s"...\n', inf$tot_reads))
		BSobj@tot_reads<- read.table.ffdf(file= inf$tot_reads, header= TRUE, colClasses= c('factor', rep('integer', ncol(dfsample)-1)), ...) 
		row.names(BSobj@tot_reads)<- BSobj@tot_reads$locus
		BSobj@tot_reads<- as.ffdf(BSobj@tot_reads[, 2:ncol(BSobj@tot_reads)])
	}

	if ('cnt_met' %in% mat){
		cat(sprintf('Reading file matrix "%s"...\n', inf$cnt_met))
		BSobj@cnt_met<- read.table.ffdf(file= inf$cnt_met, header= TRUE, colClasses= c('factor', rep('integer', ncol(dfsample)-1)), ...) 
		row.names(BSobj@cnt_met)<- BSobj@cnt_met$locus
		BSobj@cnt_met<- as.ffdf(BSobj@cnt_met[, 2:ncol(BSobj@cnt_met)])
	}
	
	if ('pct_met' %in% mat){
		cat(sprintf('Reading file matrix "%s"...\n', inf$pct_met))
		BSobj@pct_met<- read.table.ffdf(file= inf$pct_met, header= TRUE, colClasses= c('factor', rep('double', ncol(dfsample)-1)), ...) 
		row.names(BSobj@pct_met)<- BSobj@pct_met$locus
		BSobj@pct_met<- as.ffdf(BSobj@pct_met[, 2:ncol(BSobj@pct_met)])
	}
	## Some checks:
	#stopifnot(ncol(BSobj@tot_reads) == ncol(BSobj@cnt_met))
	#stopifnot(ncol(BSobj@cnt_met) == ncol(BSobj@pct_met))
	
	#stopifnot(nrow(BSobj@tot_reads) == nrow(BSobj@cnt_met))
	#stopifnot(nrow(BSobj@cnt_met) == nrow(BSobj@pct_met))
	#stopifnot(nrow(BSobj@loci) == nrow(BSobj@cnt_met))
	
	#stopifnot(colnames(BSobj@tot_reads) == colnames(BSobj@cnt_met))
	#stopifnot(colnames(BSobj@cnt_met) == colnames(BSobj@pct_met))

	## design df
	BSobj@design<- data.frame(
		library_id= colnames(BSobj@tot_reads),
		bs= rep(NA, ncol(BSobj@tot_reads))
	)
	
	cat(sprintf('\n%s samples\n', ncol(BSobj@tot_reads)))
	cat(sprintf('%s loci\n', nrow(BSobj@loci)))
	
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
        if(class(m)[1] != 'ffdf'){
            stop(sprintf('Slot "%s" is not a matrix', s))
        }
        slot(bsapp, s)<- as.ffdf(FUN(m, ...)) ## MEMO: use as.ff() if obj is ff_matrix
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

#read.bsdata<- function(prefix, gzip= TRUE, ...){
#	"Read the output files from BSmatrix.py and put them in a BSdata object
#
#	ARGS:
#	prefix:
#		Prefix of the input files produced by BSmatrix.py. It is expected to find
#			<prefix>.loci.gz
#			<prefix>.cnt_met.mat.gz
#			<prefix>.pct_met.mat.gz
#			<prefix>.tot_reads.mat.gz	
#		Input files can be either gzipped or not.
#	gzip:
#		Logical. Are the input files gzipped? (I.e. ending in .gz). Default TRUE
#	...:
#		Further arguments passed to read.table, e.g. nrows. (Not implemented)
#	VALUE:
#		BSdata object
#	"
#    BSobj<- BSdata()
#
#	if( gzip ){
#		gz= '.gz'
#	} else {
#		gz= ''
#	}
#	
#	## File suffixes
#	sfx<- list(
#		loci= paste('.loci', gz, sep= ''),
#		pct_met= paste('.pct_met.mat', gz, sep= ''),
#		cnt_met= paste('.cnt_met.mat', gz, sep= ''),
#		tot_reads= paste('.tot_reads.mat', gz, sep= '')
#	)
#	inf<- sfx
#	for(i in 1:length(inf)){
#		"Input files"
#		inf[[i]]<- paste(prefix, inf[[i]], sep= '')
#	}
#	for(f in inf){
#		if( !file.exists(f) ){
#			stop(sprintf('File "%s" not found', f))
#		}
#	}
#	cat(sprintf('Reading file "%s"...\n', inf$loci))
#	BSobj@loci<- read.table(file= inf$loci, header= TRUE, sep= '\t', colClasses= c('character', 'integer', 'integer', 'character', 'character', 'character'), nrows= -1)
#
#	cat(sprintf('Reading total file "%s"...\n', inf$tot_reads))
#	colNames<- name.reader(inf$tot_reads)
#	BSobj@tot_reads<- read.big.matrix(filename= inf$tot_reads, header= FALSE, skip= 1, sep= '\t', type= 'short', has.row.names= TRUE, col.names= colNames)
#
#	cat(sprintf('Reading file matrix "%s"...\n', inf$cnt_met))
#	colNames<- name.reader(inf$cnt_met)
#	BSobj@cnt_met<- read.big.matrix(filename= inf$cnt_met, header= FALSE, skip= 1, sep= '\t', type= 'short', has.row.names= TRUE, col.names= colNames)
#
#	cat(sprintf('Reading file matrix "%s"...\n', inf$pct_met))
#	colNames<- name.reader(inf$pct_met)
#	BSobj@pct_met<- read.big.matrix(filename= inf$pct_met, header= FALSE, skip= 1, sep= '\t', type= 'double', has.row.names= TRUE, col.names= colNames)
#
#		## Some checks:
##	stopifnot(ncol(BSobj@tot_reads) == ncol(BSobj@cnt_met))
##	stopifnot(ncol(BSobj@cnt_met) == ncol(BSobj@pct_met))
#	
##	stopifnot(nrow(BSobj@tot_reads) == nrow(BSobj@cnt_met))
##	stopifnot(nrow(BSobj@cnt_met) == nrow(BSobj@pct_met))
##	stopifnot(nrow(BSobj@loci) == nrow(BSobj@cnt_met))
#	
##	stopifnot(colnames(BSobj@tot_reads) == colnames(BSobj@cnt_met))
##	stopifnot(colnames(BSobj@cnt_met) == colnames(BSobj@pct_met))
#
#	## design df
#	BSobj@design<- data.frame(
#		library_id= colNames,
#		bs= rep(NA, length(colNames))
#	)
#	
##	cat(sprintf('\n%s samples\n', ncol(BSobj@tot_reads)))
##	cat(sprintf('%s loci\n', nrow(BSobj@loci)))
#	
#	return(BSobj)
#}


#BSdata<- setClass("BSdata", representation(
## Object to store BS data in format similar to limma
#                                        loci= 'data.frame',
#                                        tot_reads= 'big.matrix',
#                                        cnt_met= 'big.matrix',
#                                        pct_met= 'big.matrix',
#                                        design= 'data.frame'))