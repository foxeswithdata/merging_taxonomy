rm(list = ls())

library(neonUtilities)
library(tidyverse)
library(neonMicrobe)
library(ShortRead)
library(Biostrings)
library(dada2)
library(dplyr)
library(lubridate)
library(ggplot2)

source("load_meta_data_neonMicrobe.R")
source("qcMetadata.R")

### DATA from Laura filter out to get all the relevant sites
datasets <- read.csv("data/MRR_2024_dataset_descript.csv")
datasets <- datasets %>%
  filter(MRR_use == "Y")


### Following NeonMicrobe download data instructions (will find link but this is their download data vignette)
setBaseDirectory(dir=getwd())
# makeDataDirectories()

## BELOW: chose which sites to download metadata for

sites_to_DL <- datasets$Site
# sites_to_DL <- sites_to_DL[!(sites_to_DL %in% c("HARV"))]
# sites_to_DL <- "HARV" ## Test one site for now
# sites_to_DL <- c("GUAN", "HARV", "OSBS", "UNDE")


meta_16s <- downloadSequenceMetadata(sites = sites_to_DL, 
                                     startYrMo = "2016-01", endYrMo = "2021-12", 
                                     targetGene = "16S", outDir = "data/sequenceMeta",
                                     include.provisional = TRUE)

# below are fixes data for qcMetadata input ## Something was not working here with the NeonMicrobe file
class(meta_16s$setDate)
columns <- lapply(meta_16s, function(x){
  
  if(all(class(x) %in% c("POSIXct", "POSIXt"))){
    return(TRUE)
  }
  else{
    return(FALSE)
  }
})

## Make sure all "NA" or "<NA> or na inputs are empty strings 
for (i in unlist(which(columns == FALSE))){
  test <- meta_16s[i]
  test[is.na(test)] <- ""
  test[test == "<NA>"] <- ""
  test[test == "NA"] <- ""
  meta_16s[i] <- test
}



## Remove any files without a raw Data File Path
meta_16s <- meta_16s %>%
  filter(rawDataFilePath != "")

## Now quality control - this is returning to the NeonMicrobe instrictions
meta_16s_qc <- qcMetadata(meta_16s, pairedReads = "Y", rmFlagged = "Y", outDir = "data/sequenceMeta")
rm(meta_16s)

meta_16s_qc <- meta_16s_qc %>%
  filter(rawDataFilePath != "")
meta_16s_qc$date_ymd <- as.Date(format(as.Date(meta_16s_qc$collectDate, format="%Y-%m-%d %H:%M:%S"), "%Y-%m-%d"))