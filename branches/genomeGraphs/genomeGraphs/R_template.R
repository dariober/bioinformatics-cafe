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

recode<- function(codes){
    "Recode numeric indexes of --overplot to numbers from 1 to n with no gaps
    Return:
        Vector of recoded indexes in the same order as codes
    Example:
        recode(c(10, 3, 3, 8, 10))
        [1] 3 1 1 2 3
    This is necessary because layout() doesn't like gaps.
    "
    recode_df<- data.frame(op= codes, pos= 1:length(codes), rank= rank(codes, ties.method= 'min'))
    rank_idx<- data.frame(rank= sort(unique(recode_df$rank)), idx= 1:length(unique(recode_df$rank)))
    recode_df<- merge(recode_df, rank_idx, sort= FALSE)
    recode_df<- recode_df[order(recode_df$pos),]
    return(recode_df$idx)
}

setOverplotting<- function(overplot, n){
    "Extend or create overplot vector to be of length n.
    If any element in overplot == 'NA' (as string not as N/A), then all the elements
    will be 1:n.
    If not enough elements are in overplot, they will be extended with increments 
    "
    if(length(overplot) > n){
        overplot<- overplot[1:n]
    }
    if(any(overplot == 'NA')){
        overplot<- 1:n
    } else {
        overplot<- recode(as.numeric(overplot))
        diffop<- n - length(overplot)
        if(diffop > 0){
            xp<- seq(1, diffop) + max(overplot)
            overplot<- c(overplot, xp)
        }
        overplot<- overplot[1:n]
    }
    return(overplot)
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


setHeight<- function(params, pwidth){
    "Get sensible height (in cm) given the plotting paramaters and figure width.
    params:
        Data frame from which the number and type of plot panels can be extracted.
    pwidth:
        Width of the plot.
    Return:
        Float of the height in cm to be passed to pdf()
    "
    ncoverage<- length(unique(plot_params$overplot[which(plot_params$feature == 'coverage')])) ## No. plots which are coverage
    nannotation<- length(unique(plot_params$overplot[which(plot_params$feature != 'coverage')])) ## No. plots which are not coverage
    
    if(ncoverage == 0){
        # With tall plots squash it a bit more
        pheight= pwidth/10 * nannotation
    } else {
        pheight= pwidth/5 + (pwidth/6 * ncoverage) + (pwidth/10 * nannotation)
    }
    return(pheight)
}

setPlotHeights<- function(heights, plot_params, topCoef= NULL, bottomCoef= NULL){
    "Set sensible (?) values for top and bottom panel given the vector of heights
    heights:
        Numeric vector, expected to be of the same length as inputlist.
    Return:
        Numeric vector of length inputlist + 2 to be passed to layout()"
    ncoverage<- length(unique(plot_params$overplot[which(plot_params$feature == 'coverage')])) ## No. plots which are coverage
    nannotation<- length(unique(plot_params$overplot[which(plot_params$feature != 'coverage')])) ## No. plots which are not coverage
    pp<- data.frame(unique(plot_params[, c('overplot', 'feature')]), heights)
    if(nannotation > 0){
        ref<- mean(pp$heights[which(pp$feature == 'annotation')]) 
        top<- ref * 0.75
        bottom<- ref * 2
    } else {
        ref<- mean(pp$heights[which(pp$feature == 'coverage')]) / 4
        top<- ref * 0.5
        bottom<- ref * 2            
    }
#    if(is.null(topCoef)){
#        topCoef<- 0.1 / (log2(length(unique(plot_params$overplot)))+1)
#    }
#    if(is.null(bottomCoef)){
#        bottomCoef<- 0.25 / (log2(length(unique(plot_params$overplot)))+1)
#    }
#    top<-    sum(heights) * topCoef ## max(heights) / 8
#    bottom<- sum(heights) * bottomCoef ## max(heights) / 2
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
    ## height: desired height
    w<- strheight(text, cex= 1)
    cex.out<- height / w
    return(cex.out)
}

cex.for.width<- function(text, width){
    ## text: string of text
    ## width: desired width
    w<- strwidth(text, cex= 1)
    cex.out<- width / w
    return(cex.out)
}

cex.fit<- function(text, xspan, yspan){
    "Returns the (smallest) cex value that fits text to the given x and y sizes.
    The fit is not perfect (see example) but should suite most cases.
    Example:
        plot(1:10, type= 'n'); abline(v= c(1, 6), h= c(5, 6))
        cx<- cex.fit('Masticabrodo', 5, 1)
        text(x= 1, y= 5, labels= 'Masticabrodo', cex= cx, adj= c(0,0))

        cx<- cex.fit('Masticabrodo', 4, 1)
        text(x= 6, y= 5, labels= 'Masticabrodo', cex= cx, adj= c(0,0))
    "
    cex.y<- cex.for.height(text, yspan)
    cex.x<- cex.for.width(text, xspan)
    return(min(c(cex.y, cex.x)))
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
    ## Don't use this because pdf files become a pain to open 
#    segments(x0= bdg.end[gaps], x1= bdg.end[gaps], y0= rep(ybottom, length(gaps)), y1= bdg.score[gaps], ...)
#    segments(x0= bdg.start[gaps+1], x1= bdg.start[gaps+1], y0= rep(ybottom, length(gaps)), y1= bdg.score[gaps+1], ...)
    return(gaps)
}


# ------------------------------------------------------------------------------
# Intial parameters
# ------------------------------------------------------------------------------
#
# plot_params has all the parameters that apply to each plot individually
# -----------------------------------------------------------------------
plot_params<- data.frame(file_name= c(%(inputlist)s), stringsAsFactors= FALSE)
plot_params$col_line<- recycle(nrow(plot_params), c(%(col_line)s))
plot_params$lwd<- as.numeric(recycle(nrow(plot_params), c(%(lwd)s)))
plot_params$col_text_ann<- recycle(nrow(plot_params), c(%(col_text_ann)s))
plot_params$col_track<- recycle(nrow(plot_params), c(%(col_track)s))
plot_params$col_track_rev<- recycle(nrow(plot_params), c(%(col_track_rev)s))
plot_params$snames<- recycle(nrow(plot_params), c(%(names)s))        ## c() evaluates to NULL.
plot_params$col_names<- recycle(nrow(plot_params), c(%(col_names)s))
plot_params$col_mark<- recycle(nrow(plot_params), c(%(col_mark)s))
plot_params$bg<- recycle(nrow(plot_params), c(%(bg)s) )
plot_params$col_grid<- recycle(nrow(plot_params), c(%(col_grid)s) )
plot_params$ymax<- recycle(nrow(plot_params), c(%(ymax)s) )
plot_params$ymin<- recycle(nrow(plot_params), c(%(ymin)s) )
plot_params$ylab<- recycle(nrow(plot_params), c(%(ylab)s) )
plot_params$cex_lab<- as.numeric(recycle(nrow(plot_params), c(%(cex_lab)s) ))
plot_params$vheights<- as.numeric(recycle(nrow(plot_params), c(%(vheights)s) ))
plot_params$rcode<- recycle(nrow(plot_params), c(%(rcode)s) )
plot_params$overplot<- setOverplotting(c(%(overplot)s), nrow(plot_params) )

plot_params$feature<- sapply(plot_params$file_name, filename2tracktype)
plot_params$input_order<- 1:nrow(plot_params)
plot_params<- plot_params[order(plot_params$overplot),]

# Check params
# ------------
pp<- unique(plot_params[,c('overplot', 'feature')])
if(nrow(pp) != length(unique(plot_params$overplot))){
    ## Check you are not overplotting coverage & annotation
    print(plot_params[,c('file_name', 'overplot', 'feature')])
    stop('I connot overplot annotation and coverage')
}

# Global settings
# ---------------
nonbam<- '%(nonbam)s'
mpileup_grp_bed_txt<- '%(mcov)s'
cex<- %(cex)s
cex_axis<- %(cex_axis)s * cex
cex_names<- %(cex_names)s * cex 
cex_seq<- %(cex_seq)s * cex
## cex_range<- %%(cex_range)s * cex
col_nuc<- recycle(4, c(%(col_nuc)s))
no_col_bases<- ifelse("%(no_col_bases)s" == 'False', FALSE, TRUE)
regLim<- c(%(bstart)s - 1, %(bend)s) ## These are the interval extremes as found on the input --bed, before slop 
xlim<- c(%(xlim1)s - 1, %(xlim2)s)   ## Coords after slop
chrom<- '%(chrom)s'
mar_heights<- as.numeric(c(%(mar_heights)s))
title<- '%(title)s'
cex_title<- %(cex_title)s
mar<- c(%(mar)s)
pwidth<- %(pwidth)s
pheight<- %(pheight)s
maxseq<- %(maxseq)s

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
## -------------------------------------------------
if(nrow(data_df[which(is.na(data_df$totZ) == FALSE), ]) > 0){
    max_Z<- max(data_df$totZ, na.rm= TRUE)
    min_Z<- min(data_df$totZ, na.rm= TRUE)
} else {
    max_Z<- 0
    min_Z<- 0
}

## Reference bases
## ---------------
refbases<- read.table('%(refbases)s', header= TRUE, sep= '\t', stringsAsFactors= FALSE, comment.char= '',
    colClasses= c('character', 'integer', 'integer', 'character'))
refbases$base<- toupper(refbases$base)

## Plotname: region in named after the original bed feature, not the slopped one. Change-me?
## -----------------------------------------------------------------------------------------
if(title == 'None'){
    plotname<- makePlotName(chrom, regLim)
} else {
    plotname<- gsub(':region:', makePlotName(chrom, regLim), title, fixed= TRUE)
}
#} else if (grepl('^:region:.*', title, perl= TRUE)){
#    mytitle<- sub('^:region:', '', title, perl= TRUE)
#    plotname<- makePlotName(chrom, regLim)
#    plotname<- paste(plotname, mytitle, sep= '')
#} else {
#    plotname<- title
#}

## Colours for indvidual nucleotides
## ---------------------------------
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
nplots<- nrow(plot_params)

# Min and Max of each group of plots:
data_df<- unique(merge(data_df, plot_params[, c('file_name', 'overplot')], by.x= 'file_name', by.y= 'file_name', sort= FALSE))

# LAYOUT ----------------------------------------------------------------------
# For layout() you need to know what files are `coverage` and what are not. Include
# Duplicate files! 
# Use file extensions to decide whether a file is to be plotted as coverage or
# annotation

# print(plot_params)

if(all(is.na(plot_params$vheights))){
    cov_height= 1
    ann_height= 1/6
    plot_H<- ifelse(plot_params$feature == 'coverage', cov_height, ann_height)
} else {
    plot_H<- plot_params$vheights
}
x<- aggregate(1:length(plot_H), by= list(plot_params$overplot), min)$x
plot_H<- plot_H[x]

if(all(mar_heights < 0)){
    plot_heights<- setPlotHeights(plot_H, plot_params)
} else {
    plot_heights<- c(mar_heights[1], plot_H, mar_heights[2])
}
lay.mat<- as.matrix(c(1, unique(plot_params$overplot) + 1, max(plot_params$overplot) + 2), ncol= 1)

# ----------------------------- PLOT SIZE -------------------------------------

if(pheight <= 0){
    pheight<- setHeight(plot_params, pwidth)
}
pdf('%(pdffile)s', width= pwidth/2.54, height= pheight/2.54, pointsize= %(psize)s)
# png('myplot.png', width= pwidth*100, height= pheight*100, res= 480, pointsize= 7)
layout(lay.mat, heights= plot_heights)

## TOP PANEL
## ---------
par(xaxt= 'n', yaxt= 'n', bty= 'n', mar= mar, xaxs= 'i')
plot(0, ylim= c(0,100), xlim= xlim, ylab= '', xlab= '', type= 'n')

cex.title.1<- cex.for.height(plotname, 80)
cex.title.2<- cex.for.width(plotname, (xlim[2] - xlim[1]) * 0.75)
cex.title.3<- min(c(cex.title.1, cex.title.2)) * cex_title
title_H<- strheight(plotname, cex= cex.title.3)
title_W<- strwidth(plotname, cex= cex.title.3)
if(title_H > 100 | title_W > (xlim[2] - xlim[1])){
    ## If the widht or height of the title exceeds the limit, resize it.
    cex.h<- cex.for.height(plotname, 100)
    cex.w<- cex.for.width(plotname, (xlim[2] - xlim[1]))
    cex.title.4<- min(c(cex.h, cex.w)) 
} else {
    cex.title.4<- ifelse(strheight(plotname, cex= cex.title.3) > 100, cex.for.height(plotname, 100), cex.title.3)
}
text(x= mean(xlim), y= 10, labels= plotname, cex= cex.title.4, adj= c(0.5,0))

## MAIN PANELS
## ----------
drawNewPlot<- TRUE
this_plot<- NA ## plot_params$overplot[1]
current_plot<- NA
for(i in 1:nrow(plot_params)){
    this_plot<- plot_params$overplot[i]
    file_name<- plot_params$file_name[i]
    pdata<- data_df[which(data_df$file_name == file_name & data_df$overplot == this_plot), ]
    type<- plot_params$feature[i]
    if(is.na(this_plot) | is.na(current_plot) | this_plot != current_plot){
        drawNewPlot<- TRUE
        current_plot<- this_plot
    } else {
        drawNewPlot<- FALSE
    }
    pymax<- plot_params$ymax[i] ## This can be the keyword 'max', 'indiv' or a float
    pymin<- plot_params$ymin[i] ## This can be 'min' or a float.
    if(is.na(plot_params$snames[i])){
        libname<- basename(file_name)
    } else {
        libname<- plot_params$snames[i]
    }
    col4track<- ifelse(plot_params$col_track[i] == '', ifelse(type == 'coverage', makeTransparent('grey', 90), 'firebrick4'), plot_params$col_track[i])
    if(plot_params$col_track_rev[i] == '' && type == 'coverage'){
        col4track_rev<- 'pink'
    } else if (plot_params$col_track_rev[i] == '' && type != 'coverage'){
        col4track_rev<- 'firebrick4'
    } else if (plot_params$col_track_rev[i] == 'NA'){
        col4track_rev<- col4track
    } else {
        col4track_rev<- plot_params$col_track_rev[i]
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
            min_z<- min(data_df$totZ[which(data_df$overplot == this_plot)], na.rm= TRUE) ## Min and max for plot
            max_z<- max(data_df$totZ[which(data_df$overplot == this_plot)], na.rm= TRUE)
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
        if(drawNewPlot){
            par(las= 1, mar= mar, bty= 'l', xaxt= 'n', yaxt= 's', mgp= c(3, 0.7, 0), xaxs= 'i')
            plot(x= 0, type= 'n', xlab= '', ylab= '', ylim= c(pymin, pymax), xlim= xlim, yaxt= 'n') # 
            ## Reset cex_axis if output is going to be too big
            cex_yaxis<- cex_axis
            yax<- axis(side= 2, labels= FALSE, tick= TRUE, lwd= 0.5)
            maxw<- max(sapply(yax, strwidth, cex= cex_axis))
            marWidth<- par('mar')[2] * strwidth('M') ## Width of the margin in user coords
            if(maxw > marWidth){
                maxLab<- yax[which(sapply(yax, strwidth, cex= cex_axis) == maxw)][1]
                cex_yaxis<- cex.for.width(text= maxLab, marWidth)
            }
            yax<- axis(side= 2, labels= TRUE, tick= FALSE, cex.axis= cex_yaxis)
            cex_lab<- plot_params$cex_lab[i]
            labBase<- par('usr')[1] - maxw - strwidth('W', cex= cex_yaxis)*1.5 ## Where lab is going to sit in user coords
            if(cex_lab < 0){
                cex_lab<- cex_yaxis
            } else {
                cex_lab<-  cex * cex_lab
            }
            ylabYpos<- par('usr')[4] - ( (par('usr')[4] - par('usr')[3])/2 )
            text(x= par('usr')[1] - marWidth, y= ylabYpos,
                 labels= plot_params$ylab[i], cex= cex_lab, srt= 90, xpd= NA, adj= c(0.5, 0))
            rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col= plot_params$bg[i], border= 'transparent')
            grid(col= plot_params$col_grid[i], lwd= 0.5, lty= 'dotted')
            text(x= par('usr')[1] + ((par('usr')[2] - par('usr')[1])*0.01), y= par('usr')[4] * 1, adj= c(0,1), labels= libname, col= plot_params$col_names[i], cex= cex_names)
            ## Marks
            points(x= regLim, y= rep(par('usr')[3], 2), cex= 2, pch= 17, col= plot_params$col_mark[i])
            ## Additional features
            eval(parse(text= plot_params$rcode[i]))
        }
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
            lines.bdg(pdata$start, pdata$end, bdg.score= pdata$totZ, ybottom= ybottom, col= plot_params$col_line[i], lwd= plot_params$lwd[i])
        }
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
            cexHeight<- cex.for.height(lab, 98 - (thick_top + 2))
            cex.ann<- ifelse(cexHeight < cex_names, cexHeight, cex_names) ## Choose the smallest cex to avoid names too big or names overlapping features
            text(x= rowMeans(fextr[, c('start', 'end')]), y= thick_top + 2, labels= lab, adj= c(0.5,0), col= plot_params$col_text_ann[i], cex= cex.ann)
        }
        cexHeight<- cex.for.height(libname, 98 - (thick_top + 2))
        cx<- ifelse(cexHeight < cex_names, cexHeight, cex_names) ## Choose the smallest cex to avoid names too big or names overlapping features
        text(x= par('usr')[1] + ((par('usr')[2] - par('usr')[1])*0.01), y= 100, adj= c(0,1), labels= libname, col= plot_params$col_names[i], cex= cx)
        points(x= regLim, y= rep(par('usr')[3], 2), cex= 1, pch= 17, col= 'red')
    }
}
## BOTTOM PANEL
## -----------
## x-axis labels, tickmarks and range: Note very low level
par(xaxt= 'n', yaxt= 'n', bty= 'n', mar= c(mar[1], mar[2], 0, mar[4]), tcl= 0.5, xaxs= 'i')
plot(0, type= 'n', ylim= c(0, 100), xlim= xlim, xlab= '', ylab= '')
par(xaxt= 's', xaxs= 'i')
strht<- strheight('1', cex= cex_axis)
if((strht * 6) > 95){
    cex_axis<- cex.for.height(text= '1', height= 95 / 6)
    strht<- strheight('1', cex= cex_axis)
    if(cex_seq > cex_axis){
        cex_seq<- cex_axis
    }
}
x<- axis(side= 3, labels= FALSE, tick= FALSE)
xrich<- axisRich(side= 3, tick= FALSE, labels= FALSE)
segments(x0= xrich, x1= xrich, y0= 110, y1= 100 - (strht*0.5), col= 'grey20', lwd= 0.2) ## Give 110 to make sure it goes all the way to the top
segments(x0= x, x1= x, y0= 110, y1= 100 - (strht *0.7)) ## Give 110 to make sure it goes all the way to the top
text(x= x, y= 100 - (strht*1.2), labels= formatC(x, format= 'd', big.mark= ','), cex= cex_axis, adj= c(1, 1), xpd= NA)
## 
## Range
xrange<- x[length(x)] - x[1]
text(y= 100 - (strht * 3), x= mean(c(x[1], x[length(x)])), labels= paste(formatC(xrange, format= 'd', big.mark= ','), 'bp'), cex= cex_axis, adj= c(0.5, 1))
text(y= 100 - (strht * 3), x= c(x[1], x[length(x)]), labels= '|', cex= cex_axis, adj= c(0.5, 1))
## Sequence annotation
if(nrow(refbases) > 0){
    charw<- strwidth('A', cex= cex_seq,  family= 'mono', font= 1)
    if(charw > 1){
        ## Reset string size if too wide
        cex_seq<- cex.for.width('A', 1)        
    }
    text(x= refbases$end, y= 100 - (strht * 5), labels= refbases$base, adj= c(1, 1), col= '%(col_seq)s', family= 'mono', font= 1, cex= cex_seq)
}
# abline(h= 100 - (strht * 5) - strht) ## Bottom line of text
dev.off()