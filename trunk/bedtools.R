# Functions to call bedtools tool from R. 
# Funtions are named bedtools.<tool name>.
# This is am edit of http://zvfak.blogspot.co.uk/2011/02/calling-bedtools-from-r.html
# -----------------------------------------------------------------------------
bedtools.intersect<-function(a, b, opt.string=""){
    #create temp files
    a.file<- tempfile(fileext = ".bed")
    b.file<- tempfile(fileext = ".bed")
    out   <- tempfile(fileext = ".bed")
    ori_scipen<- options('scipen')
    options(scipen= 99) # not to use scientific notation when writing out
    
    #write bed formatted dataframes to tempfile
    write.table(a, file= a.file, quote=F, sep="\t", col.names=F, row.names=F)
    write.table(b, file= b.file, quote=F, sep="\t", col.names=F, row.names=F)
    
    # create the command string and call the command using system()
    command<- paste('bedtools intersect', "-a", a.file, "-b", b.file , opt.string, ">", out, sep=" ")
    cat(command,"\n")
    try(system(command))
    
    res=read.table(out, header=F, stringsAsFactors= FALSE)
    unlink(a.file);unlink(b.file);unlink(out)
    options(scipen= ori_scipen)
    return(res)
}