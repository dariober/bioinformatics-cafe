## ------------------------------------
## Functions to work with BaseSpace API
## ------------------------------------

getBasespaceToken<- function(app_name, x= '~/.basespace_login'){
    # x: Log in file to read and get the access token
    # app_name: Name of App to get Access to.
    xf<- read.table(x, sep= '\t', header= TRUE, stringsAsFactors= FALSE)
    token<- xf$access_token[which(xf$app_name == app_name)]
    if(length(token) > 1){
        stop(sprintf('More than one token found for App %s', app_name))
    }
    return(token)
}
