#!/usr/bin/env Rscript

library(tools)
# -----------------------------------------------------------------------------
# This script template read by pycoverage.RPlot
#
# TODO:
# - Read args from config file.
# -----------------------------------------------------------------------------

recycle<- function(x, y){
    ## Recycles or trim vector y to be of the same length as x
    ## Test:
    ## x<- c('a', 'b', 'c', 'd')
    ## y<- c('1', '2')
    ## ylong<- recycle(x, y) #>>>  c('1', '2', '1', '2')
    if(is.null(y)){
        yext<- rep(NA, length(x))
        return(yext)
    }
    if(length(y) >= length(x)){
        return(y[1:length(x)])
    } else{
        yext<- rep(y, times= floor(length(x) / length(y)))
        if(length(yext) < length(x)){
            yext<- c(yext, y[1:(length(x) %%%% length(y))])
        }
        return(yext)
    }
}

makeTransparent<- function(someColor, alpha=100){
    "Given a colour name (e.g. 'red'), make it transparent.
    someColor:
    Vector of colour names to make transparent e.g. c('red', 'blue')
    alpha:
    Alpha transparency. 100 fully opaque, 0 fully transparent.
    Credit: http://stackoverflow.com/questions/8047668/transparent-equivalent-of-given-color
    "
    newColor<-col2rgb(someColor)
    apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
    blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

change_colour<- function(x, ref_idx, up_colour){
    ## Function to pass to apply to reset colour 
    to_update<- as.numeric(x[ref_idx])
    x[to_update]<- up_colour
    return(x)
}

reshape_mcov<- function(mcov){
    ## Reshape mcov (file *.grp.bed.txt) to databse long form
    # Column indexes of the counts
    count_pos<- list(
        Z= grep('\\.Z$', names(mcov), perl= TRUE), ## Indexes of columns with sum of A+C+G+T+N
        A= grep('\\.A$', names(mcov), perl= TRUE),
        C= grep('\\.C$', names(mcov), perl= TRUE),
        G= grep('\\.G$', names(mcov), perl= TRUE),
        T= grep('\\.T$', names(mcov), perl= TRUE)
    )
    ## A check: all the counts above have the same length:
    if(length(unique((sapply(count_pos, length)))) != 1){
        stop('An error occured while processing the data.')
    }
    ## Reshape to database long format
    bam_names<- grep('\\.Z$', names(mcov), perl= TRUE, value= TRUE)
    bam_names<- sub('\\.Z$', '', bam_names)
    mcov_long<- mcov[rep(1:nrow(mcov), length(bam_names)), c('chrom', 'start', 'end')]
    mcov_long$file_name<- rep(bam_names, each= nrow(mcov))
    for(b in c('A', 'C', 'G', 'T', 'Z')){
        mcov_long[[b]]<- as.vector(as.matrix((mcov[, count_pos[[b]]])))
    }
    mcov_long$feature<- 'coverage'
    mcov_long$name<- NA
    mcov_long$strand<- '.'
    return(mcov_long)
}

make_colour_df<- function(mcov_long, col_nuc, col_track){
    ## Produce a colour dataframe on the basis of mcov_long and colour vector in col_nuc
    col_df<- data.frame(
        mcov_long[, c('chrom', 'start', 'end')],
        colA= rep(ifelse(is.null(col_nuc), makeTransparent('green', 95), col_nuc[1]), nrow(mcov_long)),
        colC= rep(ifelse(is.null(col_nuc), makeTransparent('blue', 95), col_nuc[2]), nrow(mcov_long)),
        colG= rep(ifelse(is.null(col_nuc), makeTransparent('orange', 95), col_nuc[3]), nrow(mcov_long)),
        colT= rep(ifelse(is.null(col_nuc), makeTransparent('red', 95), col_nuc[4]), nrow(mcov_long)),
        colZ= rep(col_track, nrow(mcov_long)),
        stringsAsFactors= FALSE
        )
    ## If the range is too wide, use only one colour:
    region_size<- max(mcov_long$end) - min(mcov_long$start)
    if(region_size > %(maxres)s) {
        col_df[,c('colA', 'colC', 'colG', 'colT', 'colZ')]<- col_track
    }
    return(col_df)
} 

update_colour_df<- function(colour_df, refbases, col_track){
    ## Update colour dataframe according to sequence. Where read nuc == reference set
    ## colour to 'grey' (or something)
    ## If refbases has no rows there will be no difference.
    ## ------------------------------------------------------------------------------
    if('%(col_all)s' == 'False'){
        ## Where you have single bases, set base to the reference
        colour_df<- merge(colour_df, refbases, by.x= c('chrom', 'start', 'end'), by.y= c('chrom', 'start', 'end'), all.x= TRUE)
        ## Now, update colour where reference == base
        cols_idx<- match(c('colA', 'colC', 'colG', 'colT'), colnames(colour_df)) ## Positions of colour columns
        ## Add a column which says which column index should be updated
        colour_df$up_idx<- ifelse(colour_df$base == 'A', cols_idx[1],
                    ifelse(colour_df$base == 'C', cols_idx[2],
                        ifelse(colour_df$base == 'G', cols_idx[3],
                           ifelse(colour_df$base == 'T', cols_idx[4], NA
                            )
                        )
                    )       
                )
        ## Position of the newly added column
        ref_idx<- which(names(colour_df) == 'up_idx')
        ## Reset colour to N where sequenced base matches reference base
        newcolours<- t(apply(colour_df, 1, change_colour, ref_idx, col_track))
        colour_df<- data.frame(newcolours, stringsAsFactors= FALSE)
        colour_df$start<- as.integer(colour_df$start)
        colour_df$end<- as.integer(colour_df$end)        
    } else {
        colour_df$base<- NA
    }
    return(colour_df)
}

setPlotHeights<- function(heights){
    "Set sensible (?) values for top and bottom panel given the vector of heights
    heights:
        Numeric vector, expected to be of the same length as inputlist.
    Return:
        Numeric vector of length inputlist + 2 to be passed to layout()"
    top<- max(heights) / 8
    bottom<- max(heights) / 2
    return(c(top, heights, bottom))
}

makePlotName<- function(data_df, xlim){
    "Make a plot name from chromosome name extracted from data_df and region limits
    from xlim.
    "
    chrom<- unique(data_df$chrom)
    xstart<- formatC(xlim[1] + 1, big.mark= ',', format= 'd')
    xend<- formatC(xlim[2], big.mark= ',', format= 'd')
    plotname<- plotname<- sprintf('%%s:%%s-%%s', chrom, xstart, xend)
    return(plotname)
}

regname2plotname<- function(regname){
    "DEPRECATED
    -----------
    Convert regname 'chrom_start_end[_name]'
    to a plotname with format 'chr:start-end [name]'
    Examples:
        regname2plotname('chr7_5566757_5566829_ACTB_2') #>>> 'chr7:5,566,757-5,566,829 ACTB_2'
        regname2plotname('chr7_5566757_5566829_ACTB')   #>>> 'chr7:5,566,757-5,566,829 ACTB'
        regname2plotname('chr7_5566757_5566829')        #>>> 'chr7:5,566,757-5,566,829'
    "
    oripen<- options('scipen')
    options(scipen= 99)
    pname<- strsplit(regname, '_')[[1]]
    chrom<- pname[1]
    start<- formatC(as.numeric(pname[2]), big.mark= ',', format= 'd')
    end<- formatC(as.numeric(pname[3]), big.mark= ',', format= 'd')
    options(scipen= oripen)
    if(length(pname) > 3){
        rname<- paste(pname[4:length(pname)], collapse= '_')
        plotname<- sprintf('%%s:%%s-%%s %%s', chrom, start, end, rname)
    } else {
        plotname<- sprintf('%%s:%%s-%%s', chrom, start, end)
    }
    return(plotname)
}

cex.for.height<- function(text, height){
    ## text: string of text
    ## width: desired width
    w<- strheight(text, cex= 1)
    cex.out<- height / w
    return(cex.out)
}

getFeatureExtremes<- function(pdata){
    "Extract from df pdata the extreme coordinates of a feature.
    pdata:
        Data frame with at least columns start, end, name, strand
    Returns:
        Data frame with columns start, and, name, strand
    "
    pMin<- aggregate(pdata[, 'start'], by= list(name= pdata$name, strand= pdata$strand), min)
    names(pMin)[ncol(pMin)]<- 'start'
    pMax<- aggregate(pdata[, 'end'], by= list(name= pdata$name, strand= pdata$strand), max)
    names(pMax)[ncol(pMax)]<- 'end'
    pextr<- merge(pMin, pMax, sort= FALSE, by.x= c('name', 'strand'), by.y= c('name', 'strand'))
    return(pextr)
}

filename2tracktype<- function(filename){
    "Determine what track type (coverage or annotation) should be assigned to
    this file given the extension
    "
    filename= sub('\\.gz', '', filename, perl= TRUE)
    ext<- file_ext(filename)
    if(ext %%in%% c('bedGraph', 'bedgraph', 'bam')){
        return('coverage')
    }
    else {
        return('annotation')
    }
}

transparent.border<- function(pdata, xlim, r){
    "Decide whether the border of rect() should be transparent (returns TRUE) or
    not (return FALSE). Decision is based on the amount of datapoints in pdata
    spanning xlim.
    pdata:
        data frame with columns start and end which will be plotted by rect()
    xlim:
        Vector of two ints with the range spanned by the x-axis
    r:
        Threshold ratio to decide whether to make border transparent.
        If pdata/xlim > r -> TRUE (border is transparent because data is dense enough)
    
    pdata<- data.frame(start= c(0, 20, 40), end= c(10, 30, 60)) ## span<- 40
    xlim<- c(0, 1000)
    r<- 1/20
        
    "
    ## pdata span: Intervals are expected to be non-overlapping (true for bam and
    ## bedgraph usually).
    ## 0    10
    ## 15   20
    span<- sum(pdata$end - pdata$start)
    xspan<- xlim[2] -xlim[1]
    xr<- span/xspan
    if(xr < r){
        "Data is sparse, border not transparent"
        return(FALSE)
    } else {
        return(TRUE)
    }
}

# ------------------------------------------------------------------------------
# Intial settings
# ------------------------------------------------------------------------------
inputlist<- c(%(inputlist)s)
nonbam<- '%(nonbam)s'
mpileup_grp_bed_txt<- '%(mcov)s'
inputlist<- c(%(inputlist)s)
cex<- %(cex)s
cex_axis<- %(cex_axis)s * cex
cex_names<- %(cex_names)s * cex 
cex_seq<- %(cex_seq)s * cex
cex_range<- %(cex_range)s * cex
cex_ann<- %(cex_ann)s * cex
col_nuc<- c(%(col_nuc)s)
col_track<- recycle(inputlist, c(%(col_track)s))
snames<- recycle(inputlist, c(%(names)s))        ## c() evaluates to NULL.
col_names<- recycle(inputlist, c(%(col_names)s))
bg<- recycle( inputlist, c(%(bg)s) )
ymax<- recycle( inputlist, c(%(ymax)s) )
ymin<- recycle( inputlist, c(%(ymin)s) )
xlim<- c(%(xlim1)s - 1, %(xlim2)s)
ylab<- recycle( inputlist, c(%(ylab)s) )
vheights<- as.numeric(recycle( inputlist, c(%(vheights)s) ))
mar<- c(%(mar)s)
regname<- '%(regname)s'
pwidth<- %(pwidth)s
pheight<- %(pheight)s

# ------------------------------------------------------------------------------
# DATA INPUT
# Coverage and annotation will all be in "data_df"
# ------------------------------------------------------------------------------

colClasses<- c('character', 'integer', 'integer', 'character', 'integer', 'integer', 'integer', 'integer', 'numeric', 'character', 'character', 'character')
colNames<- c('chrom', 'start', 'end', 'file_name', 'A', 'C', 'G', 'T', 'Z', 'feature', 'name', 'strand')
# NON BAM FILES
if( nonbam != '' ){
    data_df<- read.table(nonbam, header= FALSE, col.names= colNames, sep= '\t', stringsAsFactors= FALSE, comment.char= '', colClasses= colClasses)   
}
# BAM FILES

if( mpileup_grp_bed_txt != ''){
    header<- read.table(mpileup_grp_bed_txt, header= FALSE, sep= '\t', stringsAsFactors= FALSE, nrows= 1, comment.char= '')   
    cc<- c('character', 'integer', 'integer', rep('numeric', length(header)-3))
    mcov<- read.table(mpileup_grp_bed_txt, header= FALSE, sep= '\t', stringsAsFactors= FALSE, skip= 1, comment.char= '', colClasses= cc)
    names(mcov)<- header
    mcov2<- reshape_mcov(mcov) ## Long format
    if( nonbam != ''){
        data_df<- rbind(data_df, mcov2)
    } else {
        data_df<- mcov2
    }
    rm(mcov2)
    rm(mcov)
}
## Max of all scores, 0 if there are no rows to plot
if(nrow(data_df[which(is.na(data_df$Z) == FALSE), ]) > 0){
    max_Z<- max(data_df$Z, na.rm= TRUE)
    min_Z<- min(data_df$Z, na.rm= TRUE)
} else {
    max_Z<- 0
    min_Z<- 0
}

## Reference bases
refbases<- read.table('%(refbases)s', header= TRUE, sep= '\t', stringsAsFactors= FALSE, comment.char= '',
    colClasses= c('character', 'integer', 'integer', 'character'))

# ------------------------------------------------------------------------------
# Plotting
# ------------------------------------------------------------------------------
samples<- unique(data_df$file_name)
nplots<- length(inputlist)

# LAYOUT ----------------------------------------------------------------------
# For layout() you need to know what files are `coverage` and what are not. Include
# Duplicate files! 
# Use file extensions to decide whether a file is to be plotted as coverage or
# annotation
plot_type<- data.frame(file_name= inputlist, feature= sapply(inputlist, filename2tracktype))
plot_type$order<- 1:length(inputlist)

if(all(is.na(vheights))){
    cov_height= 1
    ann_height= 1/6
    top_height= 1/8
    bottom_height= 1/2
    plot_type$height<- ifelse(plot_type$feature == 'coverage', cov_height, ann_height)
    plot_heights<- c(top_height, plot_type$height, bottom_height)
} else {
    plot_type$height<- vheights
    plot_heights<- setPlotHeights(vheights)
}
lay.mat<- as.matrix(c(1, plot_type$order + 1, max(plot_type$order) + 2), ncol= 1)

# ----------------------------- PLOT SIZE -------------------------------------

ncoverage<- nrow(plot_type[which(plot_type$feature == 'coverage'),]) ## No. plots which are coverage
nannotation<- nrow(plot_type[which(plot_type$feature != 'coverage'),]) ## No. plots whichare not coverage

if(pheight <= 0){
    #Get sensible default values for height
    pheight= pwidth/5 + (pwidth/4 * ncoverage) + (pwidth/8 * nannotation)
}

pdf('%(pdffile)s', width= pwidth/2.54, height= pheight/2.54, pointsize= %(psize)s)
layout(lay.mat, heights= plot_heights)

## TOP PANEL
## ---------
par(xaxt= 'n', yaxt= 'n', bty= 'n', mar= mar)
plot(0, ylim= c(0,100), xlim= xlim, ylab= '', xlab= '', type= 'n')
plotname<- makePlotName(data_df, xlim)
text(x= mean(xlim), y= 10, labels= plotname, cex= cex.for.height(plotname, 80), adj= c(0.5,0))
## MAIN PANELS
## ----------

for(i in 1:nrow(plot_type)){
    file_name<- plot_type$file_name[i]
    type<- plot_type$feature[i]
    pymax<- ymax[i] ## This can be the keyword 'max', 'indiv' or a float
    pymin<- ymin[i] ## This can be 'min' or a float.
    pdata<- data_df[which(data_df$file_name == file_name), ]
    if(is.na(snames[i])){
        libname<- basename(file_name)
    } else {
        libname<- snames[i]
    }
    col4track<- ifelse(col_track[i] == '', ifelse(type == 'coverage', 'grey', 'firebrick4'), col_track[i])
    if(type == 'coverage'){
        if(nrow(pdata) > 0){
            col_df<- make_colour_df(pdata, col_nuc, col4track)
            ## Update colour dataframe according to reference
            col_df<- update_colour_df(col_df, refbases, col4track)
            pdata<- unique(merge(pdata, col_df, by.x<- c('chrom', 'start', 'end'), by.y<- c('chrom', 'start', 'end'), sort= FALSE))
            ## This is to have base colours stacked on top of each others
            Z<- pdata$Z
            A<- pdata$A
            C<- pdata$C + A
            G<- pdata$G + C
            T<- pdata$T + G
            min_z<- min(pdata$Z) ## Min and max for plot
            max_z<- max(pdata$Z)
        } else {
            ## If there is no data just get min & max for plotting.
            min_z<- 0
            max_z<- 0            
        }
        ## Set maximum for y-axt
        ## ---------------------
        if(pymax == 'indiv'){
            pymax<- max_z
        } else if (pymax == 'max') {
            pymax<- max_Z
        } else {
            pymax<- as.numeric(pymax)
        }
        ## Set minimum for y-axt
        ## ---------------------
        if(pymin == 'indiv'){
            pymin<- min_z
        } else if(pymin == 'min'){
            pymin<- ifelse(min_Z < 0, min_Z, 0)
        } else {
            pymin<- as.numeric(pymin)
        }
        ## Set up plot
        ## ----------
        par(las= 1, mar= mar, bty= 'l', xaxt= 'n', yaxt= 's', mgp= c(3, 0.7, 0))
        plot(x= 0, type= 'n', xlab= '', ylab= '', ylim= c(pymin, pymax), xlim= xlim)
        mtext(side= 2, line= 3, text= ylab[i], cex= par('cex') * cex_axis, las= 0)
        rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col= bg[i], border= 'transparent')
        if('%(nogrid)s' == 'False'){
            grid(col= 'darkgrey')
        }
        if(nrow(pdata) > 0){
            border<- ifelse(transparent.border(pdata, xlim, 1/10), 'transparent', col4track)
            rect(xleft= pdata$start, ybottom= rep(0, nrow(pdata)), xright= pdata$end, ytop= pdata$A, col= pdata$colA, border= border)
            rect(xleft= pdata$start, ybottom= A,             xright= pdata$end, ytop= C, col= pdata$colC, border= border)
            rect(xleft= pdata$start, ybottom= C,             xright= pdata$end, ytop= G, col= pdata$colG, border= border)
            rect(xleft= pdata$start, ybottom= G,             xright= pdata$end, ytop= T, col= pdata$colT, border= border)
            rect(xleft= pdata$start, ybottom= T,             xright= pdata$end, ytop= Z, col= pdata$colZ, border= border)
        }
        text(x= par('usr')[1] + ((par('usr')[2] - par('usr')[1])*0.01), y= par('usr')[4] * 1, adj= c(0,1), labels= libname, col= col_names[i], cex= cex_names)
    } else {
        ## If type is non-coverage (annotation)
        ## ------------------------------------
        par(xaxt= 'n', yaxt= 'n', bty= 'n', mar= mar)
        plot(0, type= 'n', ylim= c(0, 100), xlim= xlim, xlab= '', ylab= '')
        rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col= makeTransparent('blue', 20), border= 'transparent')
        offs<- 40
        thick_bottom<- offs - 35
        thick_top<-    offs + 35
        thin_bottom<-  offs - 20
        thin_top<-     offs + 20
        if(nrow(pdata) > 0){
            fextr<- getFeatureExtremes(pdata)
            rect(xleft= pdata$start,
                        xright= pdata$end,
                        ybottom= ifelse(pdata$feature == 'CDS', thick_bottom, thin_bottom),
                        ytop= ifelse(pdata$feature == 'CDS', thick_top, thin_top),
                        col= col4track, border= col4track)
            segments(y0= offs, y1= offs, x0= fextr$start, x1= fextr$end, col= col4track)
            text(x= rowMeans(fextr[, c('start', 'end')]), y= thick_top + 10, labels= paste(fextr$name, ifelse(fextr$strand == '.', '', fextr$strand)), adj= c(0.5,0), col= '%(col_text_ann)s', cex= cex_ann)
        }
        text(x= par('usr')[1] + ((par('usr')[2] - par('usr')[1])*0.01), y= thick_top + 10, adj= c(0,0), labels= libname, col= col_names[i], cex= cex_names) #par('usr')[4] * 1
    }
}
## BOTTOM PANEL
## -----------
## x-axis labels, tickmarks and range: Note very low level
baseline<- 90 ## Annotate bottom panel from this y-coord.
par(xaxt= 'n', yaxt= 'n', bty= 'n', mar= c(mar[1], mar[2], 0, mar[4]), tcl= 0.5)
plot(0, type= 'n', ylim= c(0, 100), xlim= xlim, xlab= '', ylab= '')
par(xaxt= 's')
x<- axis(side= 3, labels= FALSE, tick= FALSE)
segments(x0= x, x1= x, y0= 110, y1= 95) ## Give 110 to make sure it goes all the way to the top
text(x= x, y= baseline, labels= formatC(x, format= 'd', big.mark= ','), cex= cex_axis, adj= c(1, 1), xpd= NA)
w<- strheight(formatC(x[1], format= 'd', big.mark= ','), cex= cex_axis)
## 
## Range
#options(scipen= 99) ## Do not use exp notation for large ranges
xrange<- x[length(x)] - x[1]
baseline2<- baseline - (w * 2)
text(y= baseline2, x= mean(c(x[1], x[length(x)])), labels= paste(formatC(xrange, format= 'd', big.mark= ','), 'bp'), cex= cex_range, adj= c(0.5, 1))
text(y= baseline2, x= c(x[1], x[length(x)]), labels= '|', cex= cex_range, adj= c(0.5, 1))
w2<- strheight(formatC(paste(xrange, 'bp')), cex= cex_range)
## Sequence annotation
if(nrow(refbases) > 0){
    baseline3<- baseline2 - (w2 * 2)
    text(x= refbases$end, y= baseline3, labels= refbases$base, adj= c(1, 1), col= '%(col_seq)s', family= 'mono', font= 1, cex= cex_seq)
}
dev.off()
