#! /usr/bin/env Rscript

library(class)

docstring<- '
DESCRIPTION
    Convert strings of colours in hex encoding to the most similar R colour name
    Input is one or more colors passed as command line args or a file with hex codes
    one per line.
    
    Output is tab separated with <hex> <color name>

USAGE
hexToRColour.R FF0000 00FF00 0000FF 070707 070708
FF0000	red
00FF00	green
0000FF	blue
070707	gray3
070708	gray3
'

argv<- commandArgs(trailingOnly= TRUE)

if (length(argv) == 0 || argv[1] %in% c('-h', '--help')){
    cat(docstring)
    quit(save= 'no')
}

prepareColorTable<- function(){
    cols<- unique(colors())
    RGB<- col2rgb(cols)
    rgbdf<- data.frame(t(RGB))
    rownames(rgbdf)<- cols
    return(unique(rgbdf))
}

hexColToDec<- function(x){
    # x<- '#88CBED'
    x<- sub('^#', '', x, perl= TRUE)
    x<- strsplit(x, split= '')[[1]]
    stopifnot(length(x) == 6)
    rgbv<- c(paste(x[1:2], collapse= ''), paste(x[3:4], collapse= ''), paste(x[5:6], collapse= ''))
    rgbdec<- strtoi(rgbv, 16)
    return(rgbdec)
}

mapColToTable<- function(cols, clrTable){
    # cols<- hexColToDec('#88CBED')
    stopifnot(length(cols) == 3)
    outcol<- as.vector(knn(train= clrTable, test= cols, cl= rownames(clrTable), k= 1))
}

clrTable<- prepareColorTable()

if(argv[1] == '-'){
    f <- file("stdin")
    open(f)
    colv<- readLines(f)
    colv<- sub('\\n+', ' ', colv, perl= TRUE)
    colv<- unlist(strsplit(colv, ' '))
    colv<- colv[which(colv != '')]
} else {
    colv<- argv
}

for(x in colv){
    xdec<- hexColToDec(x)
    xname<- mapColToTable(xdec, clrTable)
    cat(paste(c(x, xname), collapse= '\t'))
    cat('\n')
}

quit(save= 'no')