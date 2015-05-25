#!/usr/bin/Rscript

nSlaves <- 1

#DEBUG <- TRUE
DEBUG <- FALSE

# Prepend EBI Metabolights User library 
#.libPaths( c( "/nfs/public/rw/homes/tc_cm01/metabolights/software/R/x86_64-redhat-linux-gnu-library/3.1", .libPaths() ))

# Load libraries
#library(Risa)
library(xcms)
library(IPO)
library(mzR) 

msfile <- "/nfs/production/metabolights/steffen/metabolights/IPB-2014-01/msdata.real/LTI225-09-2neg_1-E,4_01_24328.mzData"
# Debug output
args <- commandArgs(trailingOnly = TRUE)
msfile <- args[1]

directory <- "/vol/R/BioC/devel/mtbls2/inst/mzData"
accession <- "MTBLS2"

resultsprefix <- "/vol/R/BioC/devel/mtbls2/inst/IPOresults"
outdir <- paste(resultsprefix, accession, sep="/")

cat(directory, "\n")
cat(outdir, "\n")

## Create output directory
dir.create(outdir, showWarnings = FALSE)

## Obtain absolute file paths
#msfile <- paste(directory, msfile, sep="/")
    
peakpickingParameters <- getDefaultXcmsSetStartingParams('centWave')

if (DEBUG) {
   ## reduce what gets optimised
   peakpickingParameters$min_peakwidth <- 20 
   peakpickingParameters$max_peakwidth <- 20 

   peakpickingParameters$prefilter <- 3
   peakpickingParameters$value_of_prefilter <- 100 
   peakpickingParameters$mzdiff <- 0

   peakpickingParameters$ppm <- c(30,60) 
   peakpickingParameters$snthresh <- c(10,20) 
} else {
  ## more elaborate optimisation

   peakpickingParameters$prefilter <- c(3, 10)
   peakpickingParameters$value_of_prefilter <- c(100, 1000)
   peakpickingParameters$ppm <- c(25,50) 
   peakpickingParameters$snthresh <- c(5,25) 

}

iporesult <- optimizeXcmsSet(files=msfile,
                             params=peakpickingParameters, 
			     nSlaves=nSlaves,
                             subdir=paste(outdir, basename(msfile), sep="/"))

## Write information to file
parameters <- t(sapply(iporesult$best_settings$parameters, function(x) x[1]))

write.csv(cbind(msfile, parameters),
          file=paste(outdir, "/", basename(msfile), "/ipoparameters.csv", sep=""))



