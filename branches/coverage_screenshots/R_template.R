#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# This script template read by pycoverage.RPlot
# -----------------------------------------------------------------------------

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

change_colour<- function(x, ref_idx, up_colour= '%(col_cov)s'){
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

make_colour_df<- function(mcov_long, col_nuc){
    ## Produce a colour dataframe on the basis of mcov_long and colour vector in col_nuc
    col_df<- data.frame(
        mcov_long[, c('chrom', 'start', 'end')],
        colA= rep(ifelse(is.null(col_nuc), makeTransparent('green', 95), col_nuc[1]), nrow(mcov_long)),
        colC= rep(ifelse(is.null(col_nuc), makeTransparent('blue', 95), col_nuc[2]), nrow(mcov_long)),
        colG= rep(ifelse(is.null(col_nuc), makeTransparent('orange', 95), col_nuc[3]), nrow(mcov_long)),
        colT= rep(ifelse(is.null(col_nuc), makeTransparent('red', 95), col_nuc[4]), nrow(mcov_long)),
        colZ= rep('%(col_cov)s', nrow(mcov_long)),
        stringsAsFactors= FALSE
        )
    ## If the range is too wide, use only one colour:
    region_size<- max(mcov_long$end) - min(mcov_long$start)
    if(region_size > %(maxres)s) {
        col_df[,c('colA', 'colC', 'colG', 'colT', 'colZ')]<- '%(col_cov)s'
    }
    return(col_df)
} 

update_colour_df<- function(colour_df, refbases){
    ## Update colour dataframe according to sequence. Where read nuc == reference set
    ## colour to 'grey' (or something)
    ## If refbases has no rows there will be no difference.
    ## ------------------------------------------------------------------------------
    if('%(col_all)s' == 'False' && (max(colour_df$end - colour_df$start) == 1)){
        ## See if intervals have single base resolution. If not, you can't
        ## distinguidh nucleotides anyway:
        colour_df<- merge(colour_df, refbases, by.x= c('chrom', 'end'), by.y= c('chrom', 'pos'), all.x= TRUE)
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
        newcolours<- t(apply(colour_df, 1, change_colour, ref_idx))
        colour_df<- data.frame(newcolours, stringsAsFactors= FALSE)
        colour_df$start<- as.integer(colour_df$start)
        colour_df$end<- as.integer(colour_df$end)
    } else {
        colour_df$base<- NA ## Add this column even if you don't have refbases.
    }
    return(colour_df)
}

# ------------------------------------------------------------------------------
# Intial settings
# ------------------------------------------------------------------------------
nonbam<- '%(nonbam)s'
mpileup_grp_bed_txt<- '%(mcov)s'
cex.axis<- ifelse(%(cex_axis)s < 0, par('cex.axis'), %(cex_axis)s)
col_nuc<- c(%(col_nuc)s)
snames<- c(%(names)s) ## c() evaluates to NULL.
col_names<- c(%(col_names)s)
bg<- c(%(bg)s)
inputlist<- c(%(inputlist)s)

# ------------------------------------------------------------------------------
# DATA INPUT
#
# data_df will look like:
#    chrom   start     end            file_name  A  C  G  T  Z  feature name strand
#1    chr7 5566757 5566829         actb.cov.bed NA NA NA NA 10 coverage   10 <NA>
#2    chr7 5566757 5566829         actb.cov.bed NA NA NA NA 10 coverage   10 <NA>     
#3    chr7 5566757 5566829         actb.cov.bed NA NA NA NA 20 coverage   20 <NA>
#4    chr7 5566757 5566829         actb.cov.bed NA NA NA NA 10 coverage   10 <NA>
#5    chr7 5566757 5566829         actb.cov.bed NA NA NA NA 20 coverage   20 <NA>
#6    chr7 5566757 5566829         actb.cov.bed NA NA NA NA  0 coverage    0 <NA>
#7    chr7 5566777 5566829            genes.gtf NA NA NA NA NA     exon ACTB <NA>
#110  chr7 5566777 5566778 bam/ds051.actb.2.bam  0  1  0  0  1 coverage <NA> <NA>
#210  chr7 5566778 5566779 bam/ds051.actb.2.bam  0  0  0  4  4 coverage <NA> <NA>
#310  chr7 5566779 5566780 bam/ds051.actb.2.bam  0  1  0  3  4 coverage <NA> <NA>
#
# ------------------------------------------------------------------------------

# NON BAM FILES

if( nonbam != ''){
    data_df<- read.table(nonbam, header= TRUE, sep= '\t', stringsAsFactors= FALSE, comment.char= '')   
}

# BAM FILES

if( mpileup_grp_bed_txt != ''){
    header<- read.table(mpileup_grp_bed_txt, header= FALSE, sep= '\t', stringsAsFactors= FALSE, nrows= 1, comment.char= '')   
    mcov<- read.table(mpileup_grp_bed_txt, header= FALSE, sep= '\t', stringsAsFactors= FALSE, skip= 1, comment.char= '')
    names(mcov)<- header
    mcov2<- reshape_mcov(mcov) ## Long format
    data_df<- rbind(data_df, mcov2)
    rm(mcov2)
    rm(mcov)
}

# -----------------------------------------------------------------------------
## For each position/interval in mcov (grouped coverage) assign a colour to each base  
col_df<- make_colour_df(data_df, col_nuc)

## Reference bases
## ---------------
refbases<- read.table('%(refbases)s', header= TRUE, sep= '\t', stringsAsFactors= FALSE, comment.char= '')

## Update colour dataframe according to reference
col_df<- update_colour_df(col_df, refbases)
data_df<- unique(merge(data_df, col_df, by.x<- c('chrom', 'start', 'end'), by.y<- c('chrom', 'start', 'end'), sort= FALSE))

## GTF file
do.gtf<- FALSE
if('%(gtf)s' != ''){
    gtf<- read.table('%(gtf)s', header= TRUE, sep= '\t', stringsAsFactors= FALSE, comment.char= '')
    gtf<- gtf[which(gtf$type %%in%% c('start_codon', 'stop_codon') == FALSE),] ## Do not annotate start and stop codon
    if(nrow(gtf) > 0){
        do.gtf<- TRUE
    }
}

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

## Vector indexes to iterate thourgh names colours etc. to features
names_i<- 1
col_names_i<- 1
bg_i<- 1

 # LAYOUT ----------------------------------------------------------------------
# For layout() you need to know what files are `coverage` and what are not. Include
# Duplicate files!
cov_x_ann<- 8 ## Proportion of height for coverage vs annotation plot. E.g. 4 means coverage is 4x higher than annotation 
plot_type<- unique(data_df[, c('file_name', 'feature')])
plot_type$feature<- ifelse(plot_type$feature == 'coverage', 'coverage', 'annotation')
plot_type<- unique(plot_type)
plot_type<- merge(plot_type, data.frame(file_name= inputlist, order=1:length(inputlist), stringsAsFactors= FALSE))
plot_type<- plot_type[order(plot_type$order),]
plot_type$exp<- ifelse(plot_type$feature == 'coverage', cov_x_ann, 1)
lay.mat<- rep(plot_type$order, times= plot_type$exp)
lay.mat<- as.matrix(c(lay.mat, rep(max(lay.mat)+1, 3)), ncol= 1) ## Add one more 'short' plot for bottom annotation

pdf('%(pdffile)s', width= pwidth/2.54, height= pheight/2.54, pointsize= %(psize)s)
layout(lay.mat)
# print(as.vector(lay.mat))
# -----------------------------------------------------------------------------
for(i in 1:nrow(plot_type)){
    file_name<- plot_type$file_name[i]
    type<- plot_type$feature[i]
    pdata<- data_df[which(data_df$file_name == file_name), ]
    if(is.null(snames)){
        libname<- basename(file_name)
    } else {
        libname<- snames[names_i]
    }
    if(type == 'coverage'){
        ## This is to have base colours stacked on top of each others
        Z<- pdata$Z
        A<- pdata$A
        C<- pdata$C + A
        G<- pdata$G + C
        T<- pdata$T + G
        ## Set maximum for y-axt
        ## ---------------------
        if('%(ylim)s' == 'max'){
            ylim<- max(data_df$Z, na.rm= TRUE)
        } else if('%(ylim)s' == 'indiv'){
            ylim<- max(pdata$Z)
        } else {
            ylim<- as.numeric('%(ylim)s')
        }
        par(las= 1, mar= c(0, 4, 0, 1), bty= 'l', xaxt= 'n', yaxt= 's', mgp= c(3, 0.7, 0)) ## oma= c(%(oma)s), mar= c(%(mar)s)
        plot(x= 0, type= 'n', xlab= '', ylab= '', ylim= c(0, ylim), xlim= c(%(xlim1)s, %(xlim2)s), cex.axis= cex.axis)
        rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col= bg[bg_i], border= 'transparent')
        if('%(nogrid)s' == 'False'){
            grid(col= 'darkgrey')
        }
        rect(xleft= pdata$start, ybottom= rep(0, nrow(pdata)), xright= pdata$end, ytop= pdata$A, col= pdata$colA, border= 'transparent')
        rect(xleft= pdata$start, ybottom= A,             xright= pdata$end, ytop= C, col= pdata$colC, border= 'transparent')
        rect(xleft= pdata$start, ybottom= C,             xright= pdata$end, ytop= G, col= pdata$colG, border= 'transparent')
        rect(xleft= pdata$start, ybottom= G,             xright= pdata$end, ytop= T, col= pdata$colT, border= 'transparent')
        rect(xleft= pdata$start, ybottom= T,             xright= pdata$end, ytop= Z, col= pdata$colZ, border= 'transparent')
        mtext(side= 3, text= libname, adj= 0.02, line= -1, col= col_names[col_names_i], cex= %(cex_names)s)
    } else {
        par(xaxt= 'n', yaxt= 'n', bty= 'n', mar= c(0,4,0,1))
        plot(0, type= 'n', ylim= c(0, 100), xlim= c(%(xlim1)s, %(xlim2)s), xlab= '', ylab= '')

        thick_bottom<- 40 - 15
        thick_top<-    40 + 15
        thin_bottom<-  40 - 5
        thin_top<-     40 + 5
        rect(xleft= pdata$start,
                    xright= pdata$end,
                    ybottom= ifelse(pdata$feature == 'CDS', thick_bottom, thin_bottom),
                    ytop= ifelse(pdata$feature == 'CDS', thick_top, thin_top),
                    col= '%(col_ann)s', border= '%(col_ann)s')
        text(x= rowMeans(pdata[, c('start', 'end')]), y= thick_top + 5, labels= paste(pdata$name, pdata$strand), adj= c(0.5,0), col= '%(col_text_ann)s', cex= cex.axis * 0.9)
    }
    names_i<- ifelse(names_i < length(snames), names_i + 1, 1)
    col_names_i<- ifelse(col_names_i < length(col_names), col_names_i + 1, 1)
    bg_i<- ifelse(bg_i < length(bg), bg_i + 1, 1)
}
## x-axis labels, tickmarks and range: Note very low level
par(xaxt= 'n', yaxt= 'n', bty= 'n', mar= c(0,4,0,1), tcl= 0.5)
plot(0, type= 'n', ylim= c(0, 100), xlim= c(%(xlim1)s, %(xlim2)s), xlab= '', ylab= '')
par(xaxt= 's')
x<- axis(side= 3, labels= FALSE, tick= FALSE)
segments(x0= x, x1= x, y0= 100, y1= 90)
text(x= x, y= 85, labels= formatC(x, format= 'd', big.mark= ','), cex.axis= cex.axis, adj= c(0.5, 1))
## Range
wcex_range<- ifelse(%(cex_range)s <= 0, ifelse(nplots > 3, 0.66, 0.7), %(cex_range)s)
xrange<- x[length(x)] - x[1]
text(text= '|', x= c(x[1], x[length(x)]), y= cex= wcex_range)
text(labels= formatC(paste(xrange, 'bp'), format= 'd', big.mark= ','), cex= wcex_range)

old<- "
x<- axis(side= 1, labels= FALSE)
axis(labels= formatC(x, format= 'd', big.mark= ','), side= 1, at= x, cex.axis= cex.axis)
mtext(text= '%(plotname)s', cex= 0.95, outer= TRUE, side= 4, las= 0, line= 0, col= 'grey50')
mtext(text= '%(ylab)s', cex= 0.95, outer= TRUE, side= 2, las= 0, line= -0.2)
if(nrow(refbases) > 0){
    cex_seq<- ifelse(%(cex_seq)s <= 0, ifelse(nplots > 3, 0.66, 0.75), %(cex_seq)s)
    mtext(at= refbases$pos,
        side= 1,
        text= refbases$base,
        line= %(line_seq)s,
        cex= cex_seq,
        col= '%(col_seq)s',
        adj= 1,
        family= 'mono',
        font= 1)
}
## Text for range
wcex_range<- ifelse(%(cex_range)s <= 0, ifelse(nplots > 3, 0.66, 0.7), %(cex_range)s)
xrange<- x[length(x)] - x[1]
mtext(text= '|', at= c(x[1], x[length(x)]), line= %(line_range)s, side= 1, xpd= NA, cex= wcex_range)
mtext(text= formatC(paste(xrange, 'bp'), format= 'd', big.mark= ','), line= %(line_range)s, side= 1, xpd= NA, cex= wcex_range)
"
dev.off()