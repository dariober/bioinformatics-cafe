#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# This script template read by pycoverage.RPlot
#
# TODO:
# - Window compression to potentially large bed coverage files.
# - Y-axis labelling
# - Read args from config file.
# -----------------------------------------------------------------------------

recycle<- function(x, y){
    ## Recycles or trim vector y to be of the same length of x
    if(is.null(y)){
        yext<- rep(NA, length(x))
        return(yext)
    }
    if(length(y) >= length(x)){
        return(y[1:length(x)])
    } else{
        yext<- c(rep(y, times= floor(length(x) / length(y))), y[1:(length(x) %%%% length(y))] )
        return(yext)
    }
}

makeTransparent<-function(someColor, alpha=100){
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
    mcov_long$strand<- NA
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
vheights<- as.numeric(recycle( inputlist, c(%(vheights)s) ))
mar<- c(%(mar)s)

# ------------------------------------------------------------------------------
# DATA INPUT
# ------------------------------------------------------------------------------

# NON BAM FILES

if( nonbam != '' ){
    data_df<- read.table(nonbam, header= TRUE, sep= '\t', stringsAsFactors= FALSE, comment.char= '')   
}
# BAM FILES

if( mpileup_grp_bed_txt != ''){
    header<- read.table(mpileup_grp_bed_txt, header= FALSE, sep= '\t', stringsAsFactors= FALSE, nrows= 1, comment.char= '')   
    mcov<- read.table(mpileup_grp_bed_txt, header= FALSE, sep= '\t', stringsAsFactors= FALSE, skip= 1, comment.char= '')
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

## Reference bases

refbases<- read.table('%(refbases)s', header= TRUE, sep= '\t', stringsAsFactors= FALSE, comment.char= '')

# ------------------------------------------------------------------------------
# Plotting
# ------------------------------------------------------------------------------
samples<- unique(data_df$file_name)
nplots<- length(inputlist)

pwidth<- %(pwidth)s
pheight<- %(pheight)s
if(pheight <= 0){
    #Get sensible default values for height
    pheight<- (pwidth * 0.35) + ((pwidth/10) / nplots)
}

# LAYOUT ----------------------------------------------------------------------
# For layout() you need to know what files are `coverage` and what are not. Include
# Duplicate files!
plot_type<- unique(data_df[, c('file_name', 'feature')])
plot_type$feature<- ifelse(plot_type$feature == 'coverage', 'coverage', 'annotation')
plot_type<- unique(plot_type)
plot_type<- merge(plot_type, data.frame(file_name= inputlist, order=1:length(inputlist), stringsAsFactors= FALSE))
plot_type<- plot_type[order(plot_type$order),]

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

pdf('%(pdffile)s', width= pwidth/2.54, height= pheight/2.54, pointsize= %(psize)s)
layout(lay.mat, heights= plot_heights)
## Top panel
par(xaxt= 'n', yaxt= 'n', bty= 'n', mar= mar)
plot(0, ylim= c(0,100), xlim= c(%(xlim1)s, %(xlim2)s), ylab= '', xlab= '')
for(i in 1:nrow(plot_type)){
    file_name<- plot_type$file_name[i]
    type<- plot_type$feature[i]
    pdata<- data_df[which(data_df$file_name == file_name), ]
    if(is.na(snames[i])){
        libname<- basename(file_name)
    } else {
        libname<- snames[i]
    }
    col4track<- ifelse(col_track[i] == '', ifelse(type == 'coverage', 'grey', 'firebrick4'), col_track[i])
    if(type == 'coverage'){
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
        ## Set maximum for y-axt
        ## ---------------------
        pymax<- ymax[i]
        if(pymax == 'max'){
            pymax<- max(data_df$Z, na.rm= TRUE)
        } else if(pymax == 'indiv'){
            pymax<- max(pdata$Z)
        } else {
            pymax<- as.numeric(pymax)
        }
        pymin<- ymin[i]
        if(pymin == 'min'){
            pymin<- ifelse(min(data_df$Z, na.rm= TRUE) < 0, min(data_df$Z, na.rm= TRUE), 0)
        } else {
            pymin<- as.numeric(pymin)
        }
        par(las= 1, mar= mar, bty= 'l', xaxt= 'n', yaxt= 's', mgp= c(3, 0.7, 0))
        plot(x= 0, type= 'n', xlab= '', ylab= '', ylim= c(pymin, pymax), xlim= c(%(xlim1)s, %(xlim2)s), cex.axis= cex_axis)
        rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col= bg[i], border= 'transparent')
        if('%(nogrid)s' == 'False'){
            grid(col= 'darkgrey')
        }
        rect(xleft= pdata$start, ybottom= rep(0, nrow(pdata)), xright= pdata$end, ytop= pdata$A, col= pdata$colA, border= 'transparent')
        rect(xleft= pdata$start, ybottom= A,             xright= pdata$end, ytop= C, col= pdata$colC, border= 'transparent')
        rect(xleft= pdata$start, ybottom= C,             xright= pdata$end, ytop= G, col= pdata$colG, border= 'transparent')
        rect(xleft= pdata$start, ybottom= G,             xright= pdata$end, ytop= T, col= pdata$colT, border= 'transparent')
        rect(xleft= pdata$start, ybottom= T,             xright= pdata$end, ytop= Z, col= pdata$colZ, border= 'transparent')
        text(x= par('usr')[1] + ((par('usr')[2] - par('usr')[1])*0.01), y= par('usr')[4] * 1, adj= c(0,1), labels= libname, col= col_names[i], cex= cex_names)
    } else {
        par(xaxt= 'n', yaxt= 'n', bty= 'n', mar= mar)
        plot(0, type= 'n', ylim= c(0, 100), xlim= c(%(xlim1)s, %(xlim2)s), xlab= '', ylab= '')
        rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col= makeTransparent('blue', 20), border= 'transparent')
        thick_bottom<- 30 - 15
        thick_top<-    30 + 15
        thin_bottom<-  30 - 5
        thin_top<-     30 + 5
        rect(xleft= pdata$start,
                    xright= pdata$end,
                    ybottom= ifelse(pdata$feature == 'CDS', thick_bottom, thin_bottom),
                    ytop= ifelse(pdata$feature == 'CDS', thick_top, thin_top),
                    col= col4track, border= col4track)
        text(x= rowMeans(pdata[, c('start', 'end')]), y= thick_top + 10, labels= paste(pdata$name, pdata$strand), adj= c(0.5,0), col= '%(col_text_ann)s', cex= cex_ann)
        text(x= par('usr')[1] + ((par('usr')[2] - par('usr')[1])*0.01), y= par('usr')[4] * 1, adj= c(0,1), labels= libname, col= col_names[i], cex= cex_names)
    }
}

## x-axis labels, tickmarks and range: Note very low level
baseline<- 90 ## Annotate bottom panel from this y-coord.
par(xaxt= 'n', yaxt= 'n', bty= 'n', mar= c(mar[1], mar[2], 0, mar[4]), tcl= 0.5)
plot(0, type= 'n', ylim= c(0, 100), xlim= c(%(xlim1)s, %(xlim2)s), xlab= '', ylab= '')
par(xaxt= 's')
x<- axis(side= 3, labels= FALSE, tick= FALSE)
segments(x0= x, x1= x, y0= 110, y1= 95) ## Give 110 to make sure it goes all the way to the top
text(x= x, y= baseline, labels= formatC(x, format= 'd', big.mark= ','), cex= cex_axis, adj= c(1, 1), xpd= NA)
w<- strheight(formatC(x[1], format= 'd', big.mark= ','), cex= cex_axis)
## 
## Range
xrange<- x[length(x)] - x[1]
#segments(x0= c(x[1], x[length(x)]), x1= c(x[1], x[length(x)]), y0= 40, y1= 50)
baseline2<- baseline - (w * 2)
text(y= baseline2, x= mean(c(x[1], x[length(x)])), labels= formatC(paste(xrange, 'bp'), format= 'd', big.mark= ','), cex= cex_range, adj= c(0.5, 1))
text(y= baseline2, x= c(x[1], x[length(x)]), labels= '|', cex= cex_range, adj= c(0.5, 1))
w2<- strheight(formatC(paste(xrange, 'bp')), cex= cex_range)
## Sequence annotation
if(nrow(refbases) > 0){
    baseline3<- baseline2 - (w2 * 2)
    text(x= refbases$end, y= baseline3, labels= refbases$base, adj= c(0.5, 1), col= '%(col_seq)s', family= 'mono', font= 1, cex= cex_seq)
}
dev.off()
