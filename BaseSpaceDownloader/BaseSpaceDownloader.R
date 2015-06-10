#! /usr/bin/env Rscript

# TODO
# * Better reporting for failed files. Make it obvious what has failed
#   and what has succeded.

VERSION<- '0.4.0'
APP_NAME= 'Get FASTQ files' # App name to get fastq files

done<- suppressWarnings(suppressMessages(require(BaseSpaceR)))
if(done == FALSE){
    cat('Please install the "BaseSpaceR" package.\n')
    cat('See http://master.bioconductor.org/packages/release/bioc/html/BaseSpaceR.html\n\n')
    quit(save= 'no', status= 1)
}

done<- suppressWarnings(suppressMessages(require(data.table)))
if(done == FALSE){
    cat('\nPlease install the "data.table" package. Open an R session and execute:\n\n')
    cat('> install.packages("data.table")\n\n')
    quit(save= 'no', status= 1)
}


getBasespaceToken<- function(app_name= APP_NAME, x= '~/.basespace_login'){
    # x:
    #   Log in file to read and get the access token
    # app_name:
    #   Name of App to get Access to.
    #
    # Login file (~/.basespace_login) is tab separated with two columns named:
    # app_name, access_token.
    # 
    # For example (| = tab):
    # -----------------------------------
    # app_name        | access_token    
    # Get FASTQ files | dd9d20...59c5b43
    # -----------------------------------
    xf<- read.table(x, sep= '\t', header= TRUE, stringsAsFactors= FALSE)
    if(app_name %in% xf$app_name == FALSE){
        stop(sprintf("\n\nApp name '%s' not found!\n\n", app_name))
    }
    token<- xf$access_token[which(xf$app_name == app_name)]
    if(length(token) > 1){
        stop(sprintf('More than one token found for App %s', app_name))
    }
    return(token)
}

getFileSizeFromFileItem<- function(fileItem){
    # Returns the size of a file as stored in a FileItem obj.
    # fileItem:
    #      FileItem obj typically returned from a Samples obj with
    #      listFiles(s)
    # Return:
    #      Named list with file name and file size.
    fitems<- fileItem@data@Items
    stopifnot(length(fitems) == 1) ## Items is a list so presulmably it could contain more than 1
                                   ## element. Stop if it has more than 1.
    fileSize<- fitems[[1]]@Size
    fileName<- fitems[[1]]@Path # This includes Path and Name
    return(list(fileName= fileName, fileSize= fileSize))
}

compareFileSizes<- function(fileItem){
    # Compare the file size in a FileItem obj (expected) vs the size of the
    # downloaded file (observed). The downloaded file is expected to be on the
    # same path as in the fileItem obj, typically in current dir or in
    # Data/Intensities/...
    # Return:
    #   Named list with file name, exp size, obs size.  
    fileNameSize<- getFileSizeFromFileItem(fileItem)
    if(!file.exists(fileNameSize$fileName)){
        print(fileNameSize$fileName)
        print(fileItem)
        stop('File not found!')
    }
    obsSize<- file.size(fileNameSize$fileName)
    expObs<- list(fileName= fileNameSize$fileName, sizeExp= fileNameSize$fileSize, sizeObs= obsSize)
    return(expObs)
}

checkFileExists<- function(fileItem){
    # Check whether the file in FileItem obj exists on disk (i.e. jhas been
    # downloaded already) and check whether file size on disk matches file size in
    # FileItem object.
    # Return:
    #   TRUE if filename and size on disk match name and size in FileItem.
    #   FALSE otherwise.
    expFileNameSize<- getFileSizeFromFileItem(fileItem)
    if(file.exists(expFileNameSize$fileName) == FALSE){
        return(FALSE)
    }
    expObsSize<- compareFileSizes(fileItem)
    if(expObsSize$sizeExp == expObsSize$sizeObs){
        return(TRUE)
    } else {
        return(FALSE)
    }
}

getFastqFromBaseSpace<- function(
    proj_id,
    accessToken,
    dest_dir= '.' ,
    regex= '.*\\.gz$',
    sample_regex= '.*',
    echo= FALSE,
    verbose= TRUE){
    # Download fastq files for project ID from BaseSpace.
    # MEMO: Fastfile names might be slightly dofferent from BaseSpace. E.g.
    # proj_id:
    #   Project ID. You can get this from extracted project URL
    # dest_dir:
    #   Destination dir. Fastq will be in Data/....
    # accessToken:
    #   Access token. If NULL, try to read it from '~/.basespace_login' where
    #   the app line 'Get FASTQ files' is expected to be found.
    # regex:
    #   perl regexp to grep only these file names.
    # echo:
    #   If TRUE, only show which files would be downloaded.
    # Returns:
    #   Vector of downloaded files
    # TODO:
    #     This function is too complex and needs to be split.
    fileDT<- data.table(filename= NA, SampleId= NA, ExperimentName= NA)[0,]
    aAuth<- AppAuth(access_token = accessToken)
    myProj <- listProjects(aAuth)
    selProj <- Projects(aAuth, id = proj_id, simplify = TRUE) 
    sampl <- listSamples(selProj, limit= 1024) # limit appears to be constrained btw 1 and 1024. What if you have more? 

    inSampleAll <- Samples(aAuth, id = Id(sampl))
    inSample<- c()
    for(s in inSampleAll){
        if(grepl(sample_regex, s@data@SampleId, perl= TRUE)){
            inSample<- c(inSample, s)
        }
    }

    for(s in inSample){
        ff <- listFiles(s)
        for(i in 1:length(ff)){
            # Iterate thorough the FileItem objects
            f<- ff[i]
            stopifnot(length(f@data@Items) == 1)
            if( grepl(regex, Name(f), perl= TRUE) ){
                infoList<- list(filename= f@data@Items[[1]]@Name, SampleId= s@data@SampleId, ExperimentName= s@data@ExperimentName)
                fileDT<- rbindlist(list(fileDT, infoList))
                if(verbose){
                    print(f)
                }
                attempt<- 1
                if(!echo){
                    if(checkFileExists(f)){
                        cat('\nFile found, skipping download\n')
                        next
                    }
                    # Try to download file up to five times.
                    while(attempt <= 5){
                        out<- tryCatch(
                            expr= {getFiles(aAuth, id= Id(f), destDir = dest_dir, verbose = TRUE)},
                            error= function(e) {
                                cat(sprintf('\nAttempt: %s failed with\n', attempt))
                                message(e)
                                cat('\n')
                                return(NA)
                            }
                        )
                        if(is.null(out) || !is.na(out)){
                            break  # If getFiles returns NULL (older versions) or not NA download was ok.                            
                        } else {
                            attempt<- attempt + 1 # Something went wrong, try again.
                        }
                    }
                    # Check expected and downloaded file sizes:
                    expObsSize<- compareFileSizes(f)
                    if(expObsSize$sizeExp != expObsSize$sizeObs){
                        print(f)
                        stop(sprintf('\nFile %s: Expected size from FileItem object is %s but downloaded file is %s\n\n', expObsSize$fileName, expObsSize$sizeExp, expObsSize$sizeObs))
                        stop()
                    }
                }
            }
        }
    }
    return(fileDT)
}


# END_OF_FUNCTIONS <- Don't change this string it is used to source in run_test.R
# ==============================================================================
# If you just want to use the functions in BaseSpaceDownloader via source(...)
# exit after having sourced.
if(interactive()){
    cat('BaseSpaceDownloader: Functions loaded\n')
    options(show.error.messages=FALSE)
    on.exit(options(show.error.messages=TRUE))
    stop()
}
# ==============================================================================

done<- suppressWarnings(suppressMessages(require(argparse)))
if(done == FALSE){
    cat('\nPlease install the "argparse" package. Open an R session and execute:\n\n')
    cat('> install.packages("argparse")\n\n')
    quit(save= 'no', status= 1)
}

docstring<- sprintf("DESCRIPTION \\n\\
Download fastq files from BaseSpace given a project ID. \\n\\
\\n\\
EXAMPLE \\n\\
BaseSpaceDownloader.R -p 18434424 \\n\\
BaseSpaceDownloader.R -p 18434424 -r \"Ldono.*\\.gz\" \\n\\
\\n\\
Version %s", VERSION)

parser<- ArgumentParser(description= docstring, formatter_class= 'argparse.RawTextHelpFormatter')

parser$add_argument("-p", "--projid", help= "Project ID. Typically obtained from project's URL", required= TRUE)
parser$add_argument("-t", "--token", help= "Access token")
parser$add_argument("-r", "--regex", help= "Regex to filter by file name. Default all files ending in .gz", default= '.*\\.gz$')
parser$add_argument("-s", "--sample_regex", help= "Regex to filter by sample. Default all samples", default= '.*')
parser$add_argument("-e", "--echo", help= "Only show which files would be downloaded", action= 'store_true')
parser$add_argument("-o", "--outdir", help= "DEPRECATED: Output dir for fetched files.\\n", default= '.')

xargs<- parser$parse_args()
xargs$outdir<- '.'

# ==============================================================================

if(length(xargs$token) == 0){
    accessToken<- getBasespaceToken(app_name= APP_NAME)
} else {
    accessToken<- xargs$token
}

fileDT<- getFastqFromBaseSpace(
    proj_id= xargs$projid,
    accessToken= accessToken,
    dest_dir= xargs$outdir ,
    regex= xargs$regex,
    sample_regex= xargs$sample_regex,
    echo= xargs$echo,
    verbose= TRUE
)
cat('<sampleTable>\n')
write.table(fileDT, file= stdout(), sep= '\t', row.names= FALSE, quote= FALSE)
cat('</sampleTable>\n')

cat(sprintf('\n%s files found with regex "%s"\n\n', nrow(fileDT), xargs$regex))
warnings()
quit(save= 'no')
