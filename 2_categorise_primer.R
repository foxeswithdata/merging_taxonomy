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

source("primer_2_grep.R")
source("load_primers.R")

#below file has all the 15 sites from 2016 to 2021. this file may need to be changed based if downloaded on other computers
meta_16_qc <- read.csv(file = "data/sequenceMeta/mmg_metadata_16SrRNA_QCd_20240616.csv")
meta_16_qc <- meta_16_qc %>%
  filter(siteID %in% sites_to_DL)


# Process dates, add additional information
meta_16_qc$date_ymd <- as.Date(format(as.Date(meta_16_qc$collectDate, format="%Y-%m-%d %H:%M:%S"), "%Y-%m-%d"))
meta_16_qc$yr <- year(meta_16_qc$date_ymd)
meta_16_qc$mnth <- month(meta_16_qc$date_ymd)
meta_16_qc$sampleCodeID <- paste(meta_16_qc$siteID, meta_16_qc$yr, meta_16_qc$mnth, sep="_") 

meta_16_qc$processed_date_ymd <- as.Date(format(as.Date(meta_16_qc$processedDate.seq, format="%Y-%m-%d"), "%Y-%m-%d"))
meta_16_qc$processed_yr <- year(meta_16_qc$processed_date_ymd)
meta_16_qc$processed_mnth <- month(meta_16_qc$processed_date_ymd)

## Below is code to generate primer summaries according to the testProtocolVersion.seq information
test_sum <- meta_16_qc |>
  select(processed_yr, processed_mnth, laboratoryName, sequencingFacilityID, yr, mnth, domainID, siteID, testProtocolVersion.seq) |>
  dplyr::summarise(n = dplyr::n(), .by = c(processed_yr, processed_mnth, laboratoryName, sequencingFacilityID, yr, mnth, domainID, siteID, testProtocolVersion.seq)) |>
  dplyr::arrange(processed_yr)

test_sum$testProtocolVersion.seq <- as.factor(test_sum$testProtocolVersion.seq) 
unique(test_sum$testProtocolVersion.seq)
test_sum$testProtocolVersion.seq <- factor(test_sum$testProtocolVersion.seq, levels = c("BMI_markerGeneSequencingSOP_v1", 
                                                                                        "BMI_markerGeneSequencingSOP_v2", 
                                                                                        "BMI_markerGenes_16Sv4v5SOP_v1.1", 
                                                                                        "BMI_markerGenes_16Sv4v5SOP_v1.2", 
                                                                                        "GMCF_markerGenes_16SSOP_v1a"))
test_sum$testProtocolVersion.seq_num <- as.numeric(test_sum$testProtocolVersion.seq)

write.csv(test_sum, file= "dplr_protocol_version_by_site_by_date.csv")

### Further to organise by sequencerRunID

test_sum <- meta_16_qc |>
  select(processed_yr, processed_mnth, laboratoryName, sequencingFacilityID, sequencerRunID, testProtocolVersion.seq) |>
  dplyr::summarise(n = dplyr::n(), .by = c(processed_yr, processed_mnth, laboratoryName, sequencingFacilityID, sequencerRunID, testProtocolVersion.seq)) |>
  dplyr::arrange(processed_yr)

test_sum$testProtocolVersion.seq <- as.factor(test_sum$testProtocolVersion.seq) 
unique(test_sum$testProtocolVersion.seq)
test_sum$testProtocolVersion.seq <- factor(test_sum$testProtocolVersion.seq, levels = c("BMI_markerGeneSequencingSOP_v1", 
                                                                                        "BMI_markerGeneSequencingSOP_v2", 
                                                                                        "BMI_markerGenes_16Sv4v5SOP_v1.1", 
                                                                                        "BMI_markerGenes_16Sv4v5SOP_v1.2", 
                                                                                        "GMCF_markerGenes_16SSOP_v1a"))
test_sum$testProtocolVersion.seq_num <- as.numeric(test_sum$testProtocolVersion.seq)

write.csv(test_sum, file= "dplr_protocol_version_by_sequencerRunID.csv")

#### Categorise METADATA files by primer

meta_16_qc$primer_cat <- ifelse(meta_16_qc$testProtocolVersion.seq %in% c("BMI_markerGeneSequencingSOP_v2", "BMI_markerGeneSequencingSOP_v1"), "PRIMER1", "PRIMER2")
meta_16_qc$primer_cat[meta_16_qc$uid.rawFiles %in% unique_runs] <- primers_2019

meta_16_qc_primer1 <- meta_16_qc %>%
  filter(meta_16_qc$primer_cat == "PRIMER1")
meta_16_qc_primer2 <- meta_16_qc %>%
  filter(meta_16_qc$primer_cat == "PRIMER2")

meta_16_qc_primer1 <- unique(meta_16_qc_primer1)
meta_16_qc_primer2 <- unique(meta_16_qc_primer2)

write.csv(meta_16_qc_primer1, "outputs/meta_primer1_all.csv", row.names = FALSE)
write.csv(meta_16_qc_primer2, "outputs/meta_primer2_all.csv", row.names = FALSE)




######## LEGACY 
### Below is the code to test primers by actual sequence in the file 


## Step1 chose files to check
meta_16_qc <- meta_16_qc %>%
  filter(yr %in% c(2018, 2019))
unique_runs <- unique(meta_16_qc$uid.rawFiles)

# Step2 check for all the files chosen
primers_2019 <- unlist(lapply(unique_runs, function(seqRunID){
  meta_sub <- meta_16_qc[meta_16_qc$uid.rawFiles == seqRunID,]
  j = 1
  while(!file.exists(paste0("data/raw_sequence/16S/",meta_sub$rawDataFileName[j])) && j <= nrow(meta_sub)){
    j = j + 1
  }
  file_location = paste0("data/raw_sequence/16S/",meta_sub$rawDataFileName[j])
  fastq <- readFastq(file_location)
  detected = FALSE
  i = 1
  outString = "ERROR"
  while(detected == FALSE && i <= length(fastq)){
    read <- as.character(fastq[i]@sread)
    read <- substr(read, 1, 30)
    if(stringr::str_detect(read, primer_16S_fwd_primer1_grep)){
      detected = TRUE
      outString = "PRIMER1"
    }
    else if(stringr::str_detect(read, primer_16S_rev_primer1_grep)){
      detected = TRUE
      outString = "PRIMER1"
    }
    else if(stringr::str_detect(read, primer_16S_fwd_primer2_grep)){
      detected = TRUE
      outString = "PRIMER2"
    }
    else if(stringr::str_detect(read, primer_16S_rev_primer2_grep)){
      detected = TRUE
      outString = "PRIMER2"
    }
    
    else{
      detected = FALSE
      i = i + 1
    }
  }
  return (outString)
}))

## Some analysis code, can be changed later

meta_16_qc$primers <- primers_2019
meta_16_qc$direction <- ifelse(grepl("_R1", meta_16_qc$rawDataFileName), "forward", "reverse")
meta_16_qc_primer1 <- meta_16_qc %>%
  filter(meta_16_qc$primer == "PRIMER1")
meta_16_qc_sum <- meta_16_qc |>
  select(sampleCodeID, yr, mnth, siteID, direction, primers, testProtocolVersion.seq) |>
  dplyr::summarise(n = dplyr::n(), .by = c(sampleCodeID, siteID, yr, mnth, direction, primers, testProtocolVersion.seq))


