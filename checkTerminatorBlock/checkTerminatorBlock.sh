#!/bin/bash

VERSION='0.1.0'
docstring="\n
DESCRIPTION\n
Simple wrapper around picard CheckTerminatorBlock. For each file print 0 for file\n
OK or other exit code for missing terminator block (code 100).\n
\n
\n
The exit code from checkTerminatorBlock.sh is the one from last file checked.
\n
USAGE\n
checkTerminatorBlock.sh <aln.bam> [aln2.bam ...]\n
\n
REQUIRES\n
picard.jar 1.127+ on path.\n
\n
Version ${VERSION}"

# ------------------------------------------------------------------------------

bams="$*"

if [[ ${bams} = "" || ${bams} = "-h" || ${bams} == "--help" ]]
then
    echo -e $docstring
    exit 1
fi

# See http://stackoverflow.com/questions/5947742/how-to-change-the-output-color-of-echo-in-linux
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

path=`echo $PATH | sed 's/:/ /g'`

for d in $path
do
    picard=`find $d -name 'picard.jar' 2> /dev/null`
    if [[ $picard != '' ]]
    then
        break
    fi
done

if [[ $picard == '' ]]
then
    echo -e "\npicard.jar not found on PATH\n"
    exit 1
fi

#picard=`which picard.jar`
#if [[ $? != 0 ]]
#then
#    echo -e "\npicard.jar not found on PATH or non-executable\n"
#    exit 1
#fi

for bam in $bams
do
    if [[ ! -f ${bam} ]]
    then
        printf "\n${RED}File '${bam}' not found.${NC}\n\n"
        exit 1
    fi
    result=$((java -Xmx200m -jar $picard CheckTerminatorBlock I=$bam) 2>&1) # Redirect stderr to stdout then stdout to variable
    x=$?
    result=`echo $result`
    if [[ ${result} == *" HAS_TERMINATOR_BLOCK "* && ${x} == 0 ]]
    then
        col=${GREEN}
    elif [[ ${result} == *" DEFECTIVE "* ]]
    then
        col=${RED}
    else
        echo -e "\n$result\n\nUnexpected output for file $bam\n"
        exit 1
    fi
    printf "${col}$bam\t$x${NC}\n"
done

exit $x
