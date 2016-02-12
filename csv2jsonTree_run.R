###############################################################################
# csv2jsonTree_run.R, command line R tools
# Author: Martial M Sankar
# copyright: 2015-2016, SIB Swiss Institute of Bioinformatics
###############################################################################

source('csv2jsonTree_function.R')
args <- commandArgs(trailingOnly = TRUE)

inputcsv <- args[1] # provide input csv as first args,
outjson <- args[2] # provide path and filename for output json.
  
res <- csv2jsonTree(fn= inputcsv,
                    outfn = outjson, sub = NULL)