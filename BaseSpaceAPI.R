## ------------------------------------
## Functions to work with BaseSpace API
## ------------------------------------

require(BaseSpaceR)

getBasespaceToken<- function(app_name, x= '~/.basespace_login'){
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
    token<- xf$access_token[which(xf$app_name == app_name)]
    if(length(token) > 1){
        stop(sprintf('More than one token found for App %s', app_name))
    }
    return(token)
}

getFastqFromBaseSpace<- function(proj_id, dest_dir= '.', accessToken= NULL){
    # Dwonload fastq files for project ID from BaseSpace.
    # proj_id:
    #   Project ID. You can get this from extracted project URL
    # dest_dir:
    #   destination dir for fastq.
    # accessToken:
    #   Access token. If NULL, try to read it from '~/.basespace_login' where
    #   the app line 'Get FASTQ files' is expected to be found.
    stopifnot(require(BaseSpaceR))    

    if(is.null(accessToken)){
        accessToken<- getBasespaceToken(app_name= 'Get FASTQ files')
    }

    aAuth<- AppAuth(access_token = accessToken)
    myProj <- listProjects(aAuth)
    selProj <- Projects(aAuth, id = proj_id, simplify = TRUE) 
    sampl <- listSamples(selProj, limit= 1000)
    inSample <- Samples(aAuth, id = Id(sampl), simplify = TRUE)
    for(s in inSample){ 
        f <- listFiles(s, Extensions = ".gz")
        print(Name(f))
        getFiles(aAuth, id= Id(f), destDir = dest_dir, verbose = TRUE)
    }
}
