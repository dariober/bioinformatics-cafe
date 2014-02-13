#!/usr/bin/env Rscript

tryCatch({
  library(zinba)  
  }, error= function(e) {
    print("Install zinba package")
    quit(save= 'no')
  }
)

tryCatch({
  library(argparse)  
  }, error= function(e) {
    print("Install argparse package")
    quit(save= 'no')
  }
)

parser<-ArgumentParser(description= 'Run ZINBA')

parser$add_argument('--seq', help= 'Path to experimental reads', required= TRUE)
parser$add_argument('--align', help= 'Path to alignability directory', required= TRUE)
parser$add_argument('--twoBit', help= 'Path to genome build in .2bit format', required= TRUE)
parser$add_argument('--outfile', help= 'Prefix for outputted files', required= FALSE)

parser$add_argument('--refinepeaks', action='store_true', help= 'Refine peaks?')
parser$add_argument('--input', help= 'path to mapped input reads if available (default is "none")', required= FALSE)
parser$add_argument('--filetype', help= 'Either "bed", "bowtie", or "tagAlign"', required= FALSE, default= 'bed')
parser$add_argument('--threshold', help= 'FDR threshold, default is 0.05', required= FALSE, default= 0.05, type= 'double')
parser$add_argument('--numProc', help= 'Number of CPUs to use, must be less than max available   (default 1)', default= 1, type= 'integer')
parser$add_argument('--extension', help= 'Average fragment library length (size selected)', default= 200, type= 'integer')

parser$add_argument('--athresh', help= 'Number of hits per read allowed during mapping process', default= 4, type= 'integer')

args <- parser$parse_args()


quit(save= 'no')

if (args$refinepeaks == TRUE){
  refinepeaks<- 1
} else {
  refinepeaks<- 0
}

# ------------------------------------------------------------------------------
# INPUT 
# ------------------------------------------------------------------------------

## FAIRE-Seq pull down file. BED format
#pulldown_bed<- 'rhh005_5_6pctOC_3F.merged.clean.bed' 

## INPUT BED file
#input_bed<- 'rhh006_5_6pctIP_3I.merged.clean.bed'

## Alignability directory
#aln_dir<- '/lustre/sblab/berald01/reference_data/zinba_alignability/map36_hg19'

## 2bit refernce genome
#twoBitFile<- '/lustre/sblab/berald01/reference_data/genomes/hg19.2bit'

## Prefix for output files
#outprefix<- 'rhh005_5_6pctOC_3F.merged'

# ------------------------------------------------------------------------------
# ZINBA
# ------------------------------------------------------------------------------

zinba(
  refinepeaks= refinepeaks, #refine peaks? 1 for yes, 0 for no
  seq= args$seq,   #path to mapped experimental reads
  input= args$input ,  #path to mapped input reads if available (default is "none")
  filetype=  args$filetype,                                          #either 'bed', 'bowtie', or 'tagAlign'
  threshold= args$threshold,                                          #FDR threshold, default is 0.05
  align= args$align,                                            #path to alignability directory
  numProc= args$numProc,                                            #number of CPUs to use, must be less than max available   (default 1)
  twoBit= args$twoBit,                               #path to genome build in .2bit format
  outfile= args$outprefix,                           #prefix for outputted files
  extension= args$extension
)

#####################
# OPTIONAL PARAMETERS
#####################
#  basecountfile=, #path to basecount file if refinepeaks is 1
#  broad=, #broad setting, TRUE or FALSE (default)
#  printFullOut=, #print original data with enrichment estimates, 1 for yes (more space required), 0 for no (default)
#  interaction=, #whether or not to considering interaction during model selection, TRUE (default) or FALSE
#  mode=, #either "peaks" for peak calling (default) or "CNV" for calling likely amplified CNV regions for reads in "seq" (input reads are best)
#  FDR= #either TRUE (default) or FALSE. If false, then uses posterior probability to threshold peaks using 1-threshold

quit(save= 'no')

## ---------------------------------------------------
## Alignmability directory Needs to be done only once:
## ---------------------------------------------------

generateAlignability(
  mapdir= args$align,                   #mappability directory from unpacked mappability files
  outdir= args$align,                   #directory for processed files, used later in analysis
  athresh= args$athresh,                       #number of hits per read allowed during mapping process
  extension= args$extension,                   #average fragment library length
  twoBitFile= args$twoBit            #path to downloaded genome build file in .2bit format
)

