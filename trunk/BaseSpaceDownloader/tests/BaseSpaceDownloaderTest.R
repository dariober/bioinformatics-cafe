library(RUnit)

# Test project: https://basespace.illumina.com/projects/18434424/
TEST_PROJ_ID= '18434424'

test.canReadTokenCodeFromFile<- function(){
    token<- getBasespaceToken(app_name= APP_NAME, x= 'tests/basespace_login')
    checkIdentical('1234DummyTokenCode', token, "\nToken code read from file is not as expected.\n")
}
test.failIfAppNotFound<- function(){
    checkException(getBasespaceToken(app_name= 'NONSENSE', x= 'tests/basespace_login'),
        "\nShould fail if app_name is not found in login file.\n")
}

# This file and corresponding token must really exist
realToken<- getBasespaceToken(app_name= APP_NAME, x= '~/.basespace_login')

test.canGetFilesForProjID<- function(){
    fqfiles<- getFastqFromBaseSpace(proj_id= TEST_PROJ_ID, accessToken= realToken, echo= TRUE, verbose= FALSE)
    checkTrue(length(fqfiles) > 0)
    checkTrue(length(fqfiles) == 5, "\nFail might be due to files added to project\n")
    checkIdentical("Ldono-chem1_S1_L001_R1_001.fastq.gz", Name(fqfiles[[1]]), "\n\nFail might be due changes in project\n")
}

test.canGetFilesForRegex<- function(){
    fqfiles<- getFastqFromBaseSpace(regex= '.*\\.gz', proj_id= TEST_PROJ_ID, accessToken= realToken, echo= TRUE, verbose= FALSE)
    checkTrue(length(fqfiles) == 5, "\nFail might be due to files added to project\n")

    fqfiles<- getFastqFromBaseSpace(regex= 'Ldono-chem[1-2].*\\.gz', proj_id= TEST_PROJ_ID, accessToken= realToken, echo= TRUE, verbose= FALSE)
    checkTrue(length(fqfiles) == 2, fqfiles)
    
    fqfiles<- getFastqFromBaseSpace(regex= '^Ldono-chem1.*\\.gz$', proj_id= TEST_PROJ_ID, accessToken= realToken, echo= TRUE, verbose= FALSE)
    checkIdentical("Ldono-chem1_S1_L001_R1_001.fastq.gz", Name(fqfiles[[1]]), "\n\nFail might be due changes in project\n")

    fqfiles<- getFastqFromBaseSpace(regex= '.*\\.nonsense', proj_id= TEST_PROJ_ID, accessToken= realToken, echo= TRUE, verbose= FALSE)
    checkTrue(length(fqfiles) == 0, "\nThis regex should return 0 downloadable files\n")

}

#test.failIfProjNotFound<- function(){
#    checkException(
#        getFastqFromBaseSpace(proj_id= 'NONSENSE', accessToken= realToken, echo= TRUE, verbose= FALSE)
#    )
#}
