#!/usr/bin/env Rscript

VERSION<- '0.1.0'

done<- suppressWarnings(suppressMessages(require(argparse)))
if(done == FALSE){
    cat('\nPlease install the "argparse" package. Open an R session and execute:\n\n')
    cat('> install.packages("argparse")\n\n')
    quit(save= 'no', status= 1)
}

docstring<- sprintf("DESCRIPTION \\n\\
Create a text histogram of a data column. \\n\\
\\n\\
EXAMPLE \\n\\
\\n\\
Version %s", VERSION)

parser<- ArgumentParser(description= docstring, formatter_class= 'argparse.RawTextHelpFormatter')

parser$add_argument("-i", "--input", help= "File to read or - for reading stdin", required= TRUE)
parser$add_argument("-c", "--column", help= "Column index to get data for histogram", type= 'integer', required= TRUE)
parser$add_argument("-d", "--delim", help= "Column delimiter. Default tab", default= '\t')
parser$add_argument("-w", "--width", help= "Width of the highest bar in number of chars", type= 'integer', default= 50)
parser$add_argument("-b", "--breaks", help= "Number of breaks. Default is to use 'Sturges' method in R", type= 'integer', default= -1)
parser$add_argument("-V", "--verbose", help= "Print some details of the processed data", action= 'store_true')
# parser$add_argument("-H", "--header", help= "First line in input is header", action= 'store_true')
xargs<- parser$parse_args()

## Read data
## ---------
if(xargs$input == "-") {
    dat<- read.table(file("stdin"), sep= xargs$delim, header= FALSE, stringsAsFactors= FALSE, fill= TRUE)
} else {
    dat<- read.table(xargs$input, sep= xargs$delim, header= FALSE, stringsAsFactors= FALSE, fill= TRUE)
}

xdat<- suppressWarnings(as.numeric(dat[, xargs$column]))
nas<- sum(is.na(xdat))
xdat<- xdat[!is.na(xdat)]

## Some data stats
## ---------------
if(xargs$verbose){
    cat(sprintf('N. data points: %s; Removed: %s\n', length(xdat), nas))
}
if(length(xdat) == 0){
    cat('No numeric data found.\n')
    quit()
}

## Compute Histogram
## -----------------
if(xargs$breaks <= 0){
    br<- "Sturges"
} else {
    br<- xargs$breaks
}
hg<- hist(xdat, plot= FALSE, breaks= br)

## Print
## -----
# Divide the max value in n bars
bmax<- max(hg$counts)
barsPerUnit<- xargs$width / bmax
for(i in 1:length(hg$counts)){
    bar<- paste(rep('|', hg$counts[i] * barsPerUnit), collapse= '')
    xrange<- hg$mids[i]
    space<- paste(rep(' ', xargs$width - nchar(bar) + 1), collapse= '')
    cat(sprintf('%s%s%s\n', bar, space, xrange))
}
cat('+', paste(rep('-', xargs$width - 2), collapse= ''), "+", '\n', sep= '')
cat(0, paste(rep(' ', xargs$width - 2), collapse= ''), bmax, '\n', sep= '')
cat(sprintf('|: %sx\n', barsPerUnit))
quit()

