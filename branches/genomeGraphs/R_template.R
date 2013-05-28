#!/usr/bin/env Rscript

library(tools)
# -----------------------------------------------------------------------------
# This script template read by pycoverage.RPlot
#
# MEMO for column naming:
# A, C, G, T: Count of bases matching on forward. N: N matched on forward
# a, c, g, t: Coun of bases matching on reverse. n: N matched on reverse
# Z: sum of ACGTN
# z: sum of acgtn
# z + Z= depth of coverage. 
#
# TODO:
# -----------------------------------------------------------------------------

recycle<- function(n, y){
    ## Recycle or trim vector y to be of length n
    ## Test:
    ## x<- c('a', 'b', 'c', 'd')
    ## y<- c('1', '2')
    ## ylong<- recycle(length(x), y) #>>>  c('1', '2', '1', '2')
    if(is.null(y)){
        yext<- rep(NA, n)
        return(yext)
    }
    if(length(y) >= n){
        return(y[1:n])
    } else{
        yext<- rep(y, times= floor(n / length(y)))
        if(length(yext) < n){
            yext<- c(yext, y[1:(n %%%% length(y))])
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

reshape_mcov<- function(mcov, bases= count_header){
    ## Reshape mcov (file *.grp.bed.txt) to databse long form
    # Column indexes of the counts
    count_pos<- list()
    for(x in count_header){
        count_pos[[x]]<- grep(paste('\\.', x, '$', sep= ''), names(mcov), perl= TRUE)
    }
    ## A check: all the counts above have the same length:
    if(length(unique((sapply(count_pos, length)))) != 1){
        stop('An error occured while processing the data.')
    }
    ## Reshape to database long format
    bam_names<- grep('\\.Z$', names(mcov), perl= TRUE, value= TRUE)
    bam_names<- sub('\\.Z$', '', bam_names)
    mcov_long<- mcov[rep(1:nrow(mcov), length(bam_names)), c('chrom', 'start', 'end')]
    mcov_long$file_name<- rep(bam_names, each= nrow(mcov))
    for(b in bases){
        mcov_long[[b]]<- as.vector(as.matrix((mcov[, count_pos[[b]]])))
    }
    mcov_long$feature<- 'coverage'
    mcov_long$name<- NA
    mcov_long$strand<- '.'
    return(mcov_long)
}
update_row_col_df<- function(x, colour_schema, col_df){
    "Update the row x from col_df to the colour in col_match.
    base_idx:
        Index where column base is to be found
    col_df should look
    like:
               base     A     a    C    c      G      g   T   t    N    n
col_mismatch      A green green blue blue orange orange red red grey grey
col_mismatch.1    G green green blue blue orange orange red red grey grey
    "
    ## Positions to update (where 'base' matches colour_schema)
    update_row<- x
    base_idx<- which(names(col_df) == 'base')
    refbase<- toupper(x[base_idx])
    up_idx<- which(toupper(names(col_df)) == refbase)
    update_row[up_idx]<- colour_schema$col_match[which(colour_schema$base == refbase)]
    return(list(update_row))
}

make_colour_df<- function(pdata, colour_schema, refbases){
    "Create a color dataframe: Each row has the colour for each row in pdata,
    given the colour coding in colour_schema.
    refbases:
        dataframe with the reference base at each bed position
    "
    col_df<- data.frame(rbind(colour_schema$col_mismatch), stringsAsFactors= FALSE)
    names(col_df)<- colour_schema$base
    col_df<- col_df[rep(1, nrow(pdata)), ]
    col_df$chrom<- pdata$chrom  ## Coords to update on the bases of refbases
    col_df$start<- pdata$start
    col_df$end<- pdata$end
    col_df<- merge(col_df, refbases, by.x= c('chrom', 'start', 'end'), by.y= c('chrom', 'start', 'end'), all.x= TRUE, sort= FALSE)
    return(col_df)
}
update_col_df<- function(col_df, refbases){
    col_df<- apply(col_df, 1, update_row_col_df, colour_schema, col_df) ## This return a list
    col_df<- do.call( rbind, lapply(col_df, function(u) do.call(rbind, u)) ) ## List to df see http://stackoverflow.com/questions/8799990/converting-given-list-into-dataframe
    return(data.frame(col_df, stringsAsFactors= FALSE))
}

setPlotHeights<- function(heights, topCoef= NULL, bottomCoef= NULL){
    "Set sensible (?) values for top and bottom panel given the vector of heights
    heights:
        Numeric vector, expected to be of the same length as inputlist.
    Return:
        Numeric vector of length inputlist + 2 to be passed to layout()"
    top<-    sum(heights) * topCoef ## max(heights) / 8
    bottom<- sum(heights) * bottomCoef ## max(heights) / 2
    return(c(top, heights, bottom))
}

makePlotName<- function(chrom, xlim){
    "Make a plot name from chromosome name extracted from data_df and region limits
    from xlim.
    "
    chrom<- unique(chrom)
    xstart<- formatC(xlim[1] + 1, big.mark= ',', format= 'd')
    xend<- formatC(xlim[2], big.mark= ',', format= 'd')
    plotname<- plotname<- sprintf('%%s:%%s-%%s', chrom, xstart, xend)
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

axisRich<- function(side= 1, ...){
    "Increase the number of tickmarks on plot.
    side, ...:
        Arguments passed to axis()        
    "
    if(side %%in%% c(1, 3)){
        lim<- par('usr')[2] 
    } else if(side %%in%% c(2, 4)){
        lim<- par('usr')[4]
    } else {
        stop('Invalid side')
    }
    x<- axis(side= 1, labels= FALSE, tick= FALSE)
    mid<- (x[2] - x[1])/2
    halfmid<- mid/2
    x2<- sort(c(x, x + halfmid, x + mid, x + mid + halfmid))
    x2<- x2[which(x2 < lim)]
    axis(side= side, at= x2, ...)
    return(x2)
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

lines.bdg<- function(bdg.start, bdg.end, bdg.score, ybottom= min(bdg.score), ...){
    "Draw lines for bedgraph.
    bdg.start,
    bdg.end,
    bdg.score:
        Vectors of start, end and score of bedgrah features (i.e. 2nd, 3rd and 4th
        column of a bedgraph file).
    ybottom:
        When bdg features have gaps, drop vertical lines down to this baseline.
        (Typically 0 or `min(score)`)
    ...:
        Arguments passed to _segments_.
    Return:
        Indexes of where gaps start 
        Side effect of drawing segments for bedgraph profile.
    Example:
        start<- c(seq(0, 30, by= 5), seq(80, 120, by= 10), seq(150, 250, by= 25))
        end<- c(seq(5, 35, by= 5), seq(90, 130, by= 10), seq(175, 275, by= 25))
        score<- (1:length(start))^2
        plot(x= start, y= score, type= 'n', xlim= c(min(start), max(end)))
        lines.bdg(start, end, score, col= 'red', lwd= 2)
    "
    # Difference between one bdg feature and the next.
    nrows<- length(bdg.start)
    disc<- bdg.start[2:nrows] - bdg.end[1:(nrows-1)] 
    ## Where gaps are:
    gaps<- which(diff(disc) > 0) + 1
    nogaps<- which(!1:nrows %%in%% gaps)
    ## Horiz. intervals
    segments(x0= bdg.start, x1= bdg.end, y0= bdg.score, y1= bdg.score, ...)
    ## Vertical lines for adjacent features
    segments(x0= bdg.end[nogaps], x1= bdg.end[nogaps],
             y0= bdg.score[nogaps], y1= c(bdg.score[2:nrows], NA)[nogaps], ...)
    ## Vertical lines for discontinuities
    segments(x0= bdg.end[gaps], x1= bdg.end[gaps], y0= rep(ybottom, length(gaps)), y1= bdg.score[gaps], ...)
    segments(x0= bdg.start[gaps+1], x1= bdg.start[gaps+1], y0= rep(ybottom, length(gaps)), y1= bdg.score[gaps+1], ...)
    return(gaps)
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
col_nuc<- recycle(4, c(%(col_nuc)s))
no_col_bases<- ifelse("%(no_col_bases)s" == 'False', FALSE, TRUE)
col_line<- recycle(length(inputlist), c(%(col_line)s))
lwd<- as.numeric(recycle(length(inputlist), c(%(lwd)s)))
col_text_ann<- recycle(length(inputlist), c(%(col_text_ann)s))
col_track<- recycle(length(inputlist), c(%(col_track)s))
col_track_rev<- recycle(length(inputlist), c(%(col_track_rev)s))
snames<- recycle(length(inputlist), c(%(names)s))        ## c() evaluates to NULL.
col_names<- recycle(length(inputlist), c(%(col_names)s))
bg<- recycle(length(inputlist), c(%(bg)s) )
col_grid<- recycle(length(inputlist), c(%(col_grid)s) )
ymax<- recycle(length(inputlist), c(%(ymax)s) )
ymin<- recycle(length(inputlist), c(%(ymin)s) )
regLim<- c(%(bstart)s - 1, %(bend)s) ## These are the interval extremes as found on the input --bed, before slop 
xlim<- c(%(xlim1)s - 1, %(xlim2)s)   ## Coords after slop
chrom<- '%(chrom)s'
ylab<- recycle(length(inputlist), c(%(ylab)s) )
vheights<- as.numeric(recycle(length(inputlist), c(%(vheights)s) ))
title<- '%(title)s'
mar<- c(%(mar)s)
pwidth<- %(pwidth)s
pheight<- %(pheight)s
maxseq<- %(maxseq)s
rcode<- recycle(length(inputlist), c(%(rcode)s) )

# ------------------------------------------------------------------------------
# DATA INPUT
# Coverage and annotation will all be in "data_df"
# ------------------------------------------------------------------------------

# NON BAM FILES
# -------------
count_header<- c(%(count_header)s)
colClasses<- c('character', 'integer', 'integer', 'character', rep('integer', length(count_header)-1), 'numeric', 'character', 'character', 'character')
colNames<- c('chrom', 'start', 'end', 'file_name', count_header, 'feature', 'name', 'strand')
if( nonbam != '' ){
    data_df<- read.table(nonbam, header= FALSE, col.names= colNames, sep= '\t', stringsAsFactors= FALSE, comment.char= '', colClasses= colClasses)   
}

# BAM FILES
# ---------
if( mpileup_grp_bed_txt != ''){
    header<- read.table(mpileup_grp_bed_txt, header= FALSE, sep= '\t', stringsAsFactors= FALSE, nrows= 1, comment.char= '')   
    cc<- c('character', 'integer', 'integer', rep('numeric', length(header) - 3))
    mcov<- read.table(mpileup_grp_bed_txt, header= FALSE, sep= '\t', stringsAsFactors= FALSE, skip= 1, comment.char= '', colClasses= cc)
    names(mcov)<- header
    mcov2<- reshape_mcov(mcov, count_header) ## Long format
    if( nonbam != ''){
        data_df<- rbind(data_df, mcov2)
    } else {
        data_df<- mcov2
    }
    rm(mcov2)
    rm(mcov)
}
data_df$totZ<- data_df$Z + data_df$z ## Total depth= sum of ACTGN on plus and minus

## Max of all scores, 0 if there are no rows to plot
if(nrow(data_df[which(is.na(data_df$totZ) == FALSE), ]) > 0){
    max_Z<- max(data_df$totZ, na.rm= TRUE)
    min_Z<- min(data_df$totZ, na.rm= TRUE)
} else {
    max_Z<- 0
    min_Z<- 0
}

## Reference bases
refbases<- read.table('%(refbases)s', header= TRUE, sep= '\t', stringsAsFactors= FALSE, comment.char= '',
    colClasses= c('character', 'integer', 'integer', 'character'))
refbases$base<- toupper(refbases$base)

## Plotname: region in named after the original bed feature, not the slopped one. Change-me?
if(title == ''){
    plotname<- makePlotName(chrom, regLim)
} else if (grepl('^:region:.*', title, perl= TRUE)){
    mytitle<- sub('^:region:', '', title, perl= TRUE)
    plotname<- makePlotName(chrom, regLim)
    plotname<- paste(plotname, mytitle, sep= '')
} else {
    plotname<- title
}

## Colours for indvidual nucleotides
defaulColours<- c(makeTransparent('green', 80), makeTransparent('blue', 80), makeTransparent('orange', 80), makeTransparent('red', 80))
for(i in 1:length(col_nuc)){
    if(col_nuc[i] == ''){
        col_nuc[i]<- defaulColours[i]
    }
}
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
#    top_height= 1/8
#    bottom_height= 1/2
    plot_type$height<- ifelse(plot_type$feature == 'coverage', cov_height, ann_height)
#    plot_heights<- c(top_height, plot_type$height, bottom_height)
} else {
    plot_type$height<- vheights
#    plot_heights<- setPlotHeights(vheights)
}
plot_heights<- setPlotHeights(plot_type$height, topCoef= 0.1 / (log2(nplots)+1), bottomCoef= 0.25 / (log2(nplots)+1))
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
par(xaxt= 'n', yaxt= 'n', bty= 'n', mar= mar, xaxs= 'i')
plot(0, ylim= c(0,100), xlim= xlim, ylab= '', xlab= '', type= 'n')
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
    if(col_track_rev[i] == '' && type == 'coverage'){
        col4track_rev<- 'pink'
    } else if (col_track_rev[i] == '' && type != 'coverage'){
        col4track_rev<- 'firebrick4'
    } else if (col_track_rev[i] == 'NA'){
        col4track_rev<- col4track
    } else {
        col4track_rev<- col_track_rev[i]
    }
    if(type == 'coverage'){
        if(nrow(pdata) > 0){
            ## Need to decide which colour schema to use:
            if((xlim[2] - xlim[1]) < maxseq & no_col_bases == FALSE){
                colour_schema<- data.frame(
                    base= c('A', 'a', 'C', 'c', 'G', 'g', 'T', 't', 'N', 'n'),
                    col_match= rep(col4track, 10),
                    col_mismatch= c(rep(col_nuc, each= 2), col4track, col4track),
                    stringsAsFactors= FALSE
                )
            } else {
                colour_schema<- data.frame(
                    base= c('A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n'),
                    col_match= rep(c(col4track, col4track_rev), each= 5),
                    col_mismatch= rep(c(col4track, col4track_rev), each= 5),
                    stringsAsFactors= FALSE
                )
            }
            col_df<- make_colour_df(pdata, colour_schema, refbases)
            col_df<- update_col_df(col_df, refbases)
            min_z<- min(pdata$totZ) ## Min and max for plot
            max_z<- max(pdata$totZ)
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
        ## -----------
        par(las= 1, mar= mar, bty= 'l', xaxt= 'n', yaxt= 's', mgp= c(3, 0.7, 0), xaxs= 'i')
        plot(x= 0, type= 'n', xlab= '', ylab= '', ylim= c(pymin, pymax), xlim= xlim, cex.axis= cex_axis)
        mtext(side= 2, line= 3, text= ylab[i], cex= par('cex') * cex_axis, las= 0)
        rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col= bg[i], border= 'transparent')
        grid(col= col_grid[i], lwd= 0.5, lty= 'dotted')
        
        ## Draw plot 
        ## ---------
        if(nrow(pdata) > 0){
            border<- ifelse(transparent.border(pdata, xlim, 1/10), 'transparent', col4track)
            ytop<-    rep(0, nrow(pdata))
            ybottom<- rep(0, nrow(pdata))
            if(file_ext(file_name) == 'bam'){
                ## For BAM files
                for(x in colour_schema$base){
                    ## Stack bases one on top of the other
                    ytop<- ytop + pdata[[x]]
                    rect(xleft= pdata$start, ybottom= ybottom, xright= pdata$end, ytop= ytop, col= col_df[[x]], border= border)
                    ybottom<- ytop
                }
            } else {
                ## For bedgraph:
                rect(xleft= pdata$start, ybottom= ybottom, xright= pdata$end, ytop= pdata$totZ, col= col4track, border= border)
            }
            lines.bdg(pdata$start, pdata$end, bdg.score= pdata$totZ, ybottom= ybottom, col= col_line[i], lwd= lwd[i])
        }
        text(x= par('usr')[1] + ((par('usr')[2] - par('usr')[1])*0.01), y= par('usr')[4] * 1, adj= c(0,1), labels= libname, col= col_names[i], cex= cex_names)
    } else {
        ## If type is non-coverage (annotation)
        ## ------------------------------------
        par(xaxt= 'n', yaxt= 'n', bty= 'n', mar= mar, xaxs= 'i')
        plot(0, type= 'n', ylim= c(0, 100), xlim= xlim, xlab= '', ylab= '', lwd= 0.2)
        rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col= makeTransparent('blue', 20), border= 'transparent')
        offs<- 35
        thick_bottom<- offs - 30
        thick_top<-    offs + 30
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
            lab<- paste(fextr$name, ifelse(fextr$strand == '.', '', fextr$strand))
            text(x= rowMeans(fextr[, c('start', 'end')]), y= thick_top + 2, labels= lab, adj= c(0.5,0), col= col_text_ann[i], cex= cex.for.height(lab, 98 - (thick_top + 2)))
        }
        text(x= par('usr')[1] + ((par('usr')[2] - par('usr')[1])*0.01), y= thick_top + 10, adj= c(0,0), labels= libname, col= col_names[i], cex= cex_names) #par('usr')[4] * 1
    }
    points(x= regLim, y= rep(par('usr')[3], 2), cex= 2, pch= 17, col= 'red')
    eval(parse(text= rcode[i]))
}
## BOTTOM PANEL
## -----------
## x-axis labels, tickmarks and range: Note very low level
baseline<- 85 ## Annotate bottom panel from this y-coord.
par(xaxt= 'n', yaxt= 'n', bty= 'n', mar= c(mar[1], mar[2], 0, mar[4]), tcl= 0.5, xaxs= 'i')
plot(0, type= 'n', ylim= c(0, 100), xlim= xlim, xlab= '', ylab= '')
par(xaxt= 's', xaxs= 'i')
x<- axis(side= 3, labels= FALSE, tick= FALSE)
xrich<- axisRich(side= 3, tick= FALSE, labels= FALSE)
segments(x0= xrich, x1= xrich, y0= 110, y1= 95, col= 'grey20', lwd= 0.2) ## Give 110 to make sure it goes all the way to the top
segments(x0= x, x1= x, y0= 110, y1= 90) ## Give 110 to make sure it goes all the way to the top
text(x= x, y= baseline, labels= formatC(x, format= 'd', big.mark= ','), cex= cex_axis, adj= c(1, 1), xpd= NA)
w<- strheight(formatC(x[1], format= 'd', big.mark= ','), cex= cex_axis)
## 
## Range
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

# -----------------------------------------------------------------------------
#                                  TRITUME
# ----------------------------------------------------------------------------- 
#
#addGapsXY<- function(x, y, gap_size= 1, fill.y= NA){
#    "Identify gaps in x and repeat x,min(y) where gaps occur.
#    x:
#        Vector of x-coordinates where gaps are to be identifies
#    y:
#        Vector of y-coordinates
#    gap_size:
#        Distance between adjacent x-points to consider them as gap. 
#    fill.y:
#        Value to use for y in coorespondence of x gaps.
#    Return:
#        Dataframe with x, y (coords) with NA where gaps occur in y-coords in gaps
#        are filled with .
#    Examples:
#        ## xy data with gaps
#        xv<- c(1:10, 21:30, 35:40, 50:60)
#        yv<- round(runif(length(xv), min= 0, max= 100))
#        ## Fill and mark gaps
#        (xdat<- addGapsXY(xv, yv))
#    "
#    gaps<- which(diff(x) > gap_size) ## Where gaps are.
#    gap_start<- gaps + 0.1
#    gap_idx<- order(c(1:length(x), gap_start))
#    gap_val<- sort(c(1:length(x), gap_start))
#    x_fill<- c(x, x[gaps], x[gaps+1])[gap_idx]  ## x plus x-values to repeat
#    y_fill<- c(y, rep(fill.y, length(gaps)*2))[gap_idx]
#    gap_pos<- c(rep(TRUE, length(y)), rep(c(FALSE, FALSE), length(gaps)))[gap_idx]
#    xdat<- data.frame(x= x_fill, y= y_fill, gaps= gap_pos, stringsAsFactors= FALSE)
#    return(xdat)
#}
#
#gap.points<- function(xdat, ybottom= par('usr')[3], type= 'l', ...){
#    "Plot datapoints in dataframe xdat so that gaps go
#    down to base y-axis. Plotting done by lines first and segments then.
#    xdat:
#        Dataframe produced by addGapsXY
#    ybottom:
#        Draw vertical drops down to this y-coord (default down to y-axis bottom
#        limit)
#    ...:
#        Further arguments passed to points.
#    Examples:
#        ## xy data with gaps
#        xv<- c(1:10, 21:30, 35:40, 50:60)
#        yv<- round(runif(length(xv), min= 0, max= 100))
#        ## Fill and mark gaps
#        xdat<- addGapsXY(xv, yv)
#        plot(xv, yv, type= 'l', ylim= c(0, max(yv)), col= 'grey95')
#        gap.points(xdat, col= 'red', lwd= 2)
#    "
#    lines(xdat$x, xdat$y, ...)
#    if(type %%in%% c('l', 'o', 'b')){
#        gaps<- which(xdat$gaps == FALSE)
#        xna<- xdat[c(which(xdat$gaps == FALSE) - 1, which(xdat$gaps == FALSE) + 2), ]
#        xna<- xna[xna$gaps,]
#        segments(x0= c(xdat$x[gaps-1], xdat$x[gaps+1]),
#                 x1= c(xdat$x[gaps-1], xdat$x[gaps+1]),
#                 y0= c(xdat$y[gaps-1], xdat$y[gaps+1]),
#                 y1= rep(ybottom, length(gaps)*2),
#                 ...
#        )
#    }
#}
