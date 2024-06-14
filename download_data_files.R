rm(list = ls())

library(neonUtilities)
library(tidyverse)
library(neonMicrobe)
library(ShortRead)
library(Biostrings)
library(dada2)
library(dplyr)
library(ggplot2)

### NeonMicrobe Vignette for reference: https://people.ucsc.edu/~claraqin/analyze-neon-greatplains-16s.R

### DATA from Laura filter out to get all the relevant sites
datasets <- read.csv("data/MRR_2024_dataset_descript.csv")
datasets <- datasets %>%
  filter(MRR_use == "Y")


### Following NeonMicrobe download data instructions (will find link but this is their download data vignette)
setBaseDirectory(dir=getwd())
sites_to_DL <- datasets$Site
sites_to_DL <- "HARV" ## Test one site for now

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
  message(i)
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
meta_16s_qc <- meta_16s_qc %>%
  filter(rawDataFilePath != "")

## This didn't work wo I wrote all the file names in toDL.csv and wrote a script to download all files
## This script is "data/raw_sequence/downloadAll.sh 
## run from terminal and make sure files go into data/raw_sequence/16S when finished
downloadRawSequenceData(meta_16s_qc, outDir = ".", overwrite = TRUE, verbose = TRUE)

meta_16s_qc_toDL <- select(meta_16s_qc, rawDataFilePath)
write.csv(meta_16s_qc_toDL, quote = FALSE, row.names = FALSE,  "data/raw_sequence//toDL.csv")


## Now use the Harvard sites only - metadata input file might need to be adjusted for other sites
meta_16_qc_HARV <- read.csv(file = "data/sequenceMeta/mmg_metadata_16SrRNA_QCd_20240614.csv")
# meta_16_qc <- read.csv(file = "data/sequenceMeta/mmg_metadata_16SrRNA_QCd_20240603.csv")
meta_16_qc_HARV$date_ymd <- as.Date(format(as.Date(meta_16_qc_HARV$collectDate, format="%Y-%m-%d %H:%M:%S"), "%Y-%m-%d"))
# meta_16_qc$date_ymd <- as.Date(format(as.Date(meta_16_qc$collectDate, format="%Y-%m-%d %H:%M:%S"), "%Y-%m-%d"))


meta_16_qc_HARV_pre_2019 <- meta_16_qc_HARV %>%
  filter(meta_16_qc_HARV$date_ymd < ymd("2020-01-01"))
meta_16_qc_HARV_post_2019 <- meta_16_qc_HARV %>%
  filter(meta_16_qc_HARV$date_ymd >= ymd("2020-01-01"))


meta_16s_qc_toRemove <- meta_16_qc_HARV %>%
  filter(date_ymd < ymd("2016-01-01"))


# unique_runs <- unique(meta_16_qc$sequencerRunID)

# fl_nm
# fl_nm <- meta_16_qc$rawDataFileName[grepl("BMI_Plate6", meta_16_qc$rawDataFileName, ignore.case = TRUE)]



NEONMICROBE_DIR_SEQUENCE()
####### FOR PrE 2019 DATASETS ---
# from: https://data.neonscience.org/documents/-/document_library_display/kV4WWrbEEM2s/view_file/3370805?_110_INSTANCE_kV4WWrbEEM2s_redirect=https%3A%2F%2Fdata.neonscience.org%2Fdocuments%2F-%2Fdocument_library_display%2FkV4WWrbEEM2s%2Fview%2F2237401%3F_110_INSTANCE_kV4WWrbEEM2s_keywords%3D%26_110_INSTANCE_kV4WWrbEEM2s_topLink%3Dhome%26_110_INSTANCE_kV4WWrbEEM2s_advancedSearch%3Dfalse%26_110_INSTANCE_kV4WWrbEEM2s_delta2%3D20%26_110_INSTANCE_kV4WWrbEEM2s_cur2%3D2%26p_r_p_564233524_resetCur%3Dfalse%26_110_INSTANCE_kV4WWrbEEM2s_andOperator%3Dtrue%26_110_INSTANCE_kV4WWrbEEM2s_delta1%3D20
# Table 1
# Accessed June 14 2024
primer_16S_fwd_pre_2019 = "CCTACGGGNBGCASCAG"
primer_16S_rev_pre_2019 = "GGACTACNVGGGTATCTAATCC"
trim_trackReads <- trimPrimers16S(
  meta_16_qc_HARV_pre_2019$rawDataFileName, in_subdir = "raw", out_subdir = "1_trimmed/pre_2019", meta = meta_16_qc, 
  multithread = FALSE,# set multithread = FALSE on Windows computers though
  primer_16S_fwd =primer_16S_fwd_pre_2019, 
  primer_16S_rev = primer_16S_rev_pre_2019
)
####### FOR POST 2019 DATASETS --- 
# from: 
# "Improved Bacterial 16S rRNA Gene (V4 and V4-5) and Fungal Internal T
# Transcribed Spacer Marker Gene Primers for Microbial Community Surveys
# Waletrs et al 2015
# (Table 1 rows 3 for FWD and 5 for REV)
# 515f Modiﬁed GTGYCAGCMGCCGCGGTAA Parada et al. (4)
# 926r CCGYCAATTYMTTTRAGTTT Parada et al. (4)
primer_16S_fwd_post_2019 = "GTGYCAGCMGCCGCGGTAA"
primer_16S_rev_post_2019 = "CCGYCAATTYMTTTRAGTTT"
trim_trackReads <- trimPrimers16S(
  meta_16_qc_HARV_post_2019$rawDataFileName, in_subdir = "raw", out_subdir = "1_trimmed/post_2019", meta = meta_16_qc, 
  multithread = FALSE,# set multithread = FALSE on Windows computers though
  primer_16S_fwd = primer_16S_fwd_post_2019, 
  primer_16S_rev = primer_16S_rev_post_2019
)

### COMBINDED TRIM 
# Reference email from Laura Super/Hugh XXX in 2024
# Use fwd primer post 2019 and rev primer for pre 2019
trim_trackReads <- trimPrimers16S(
  meta_16_qc_HARV$rawDataFileName, in_subdir = "raw", out_subdir = "1_trimmed/combined", meta = meta_16_qc, 
  multithread = FALSE,# set multithread = FALSE on Windows computers though
  primer_16S_fwd = primer_16S_fwd_post_2019, 
  primer_16S_rev = primer_16S_rev_pre_2019
)



filter_trackReads <- qualityFilter16S(
  meta_16_qc_HARV_pre_2019$rawDataFileName, in_subdir = "1_trimmed/pre_2019", out_subdir = "2_filtered/pre_2019",
  meta = meta_16_qc_HARV_pre_2019, truncLen = 220, maxEE = 8, multithread = FALSE
)
filter_trackReads <- qualityFilter16S(
  meta_16_qc_HARV_post_2019$rawDataFileName, in_subdir = "1_trimmed/post_2019", out_subdir = "2_filtered/pre_2019",
  meta = meta_16_qc_HARV_post_2019, truncLen = 220, maxEE = 8, multithread = FALSE
)
filter_trackReads <- qualityFilter16S(
  meta_16_qc_HARV$rawDataFileName, in_subdir = "1_trimmed/combined", out_subdir = "2_filtered/combined",
  meta = meta_16_qc_HARV, truncLen = 220, maxEE = 8, multithread = FALSE
)










unique_runs <- unique(meta_16_qc_HARV$sequencerRunID)

for(i in 1:length(unique_runs)) {
  meta_thisrun <- meta_16_qc_HARV[which(meta_16_qc_HARV$sequencerRunID==unique_runs[i]),]
  fl_nm_thisrun <- meta_thisrun$rawDataFileName
  dada_out <- runDada16S(
    fl_nm_thisrun, in_subdir = "2_filtered/combined", meta = meta_16_qc_HARV,
    out_seqtab = file.path(NEONMICROBE_DIR_OUTPUTS(), "HARV/combined",
                           paste0("HARV_asv_", unique_runs[i], ".Rds")),
    out_track = file.path(NEONMICROBE_DIR_OUTPUTS(), "HARV/combined",
                          paste0("HARV_track_", unique_runs[i], ".csv")),
    verbose = FALSE,
    multithread = FALSE
  )
}

seqtab_joined <- mergeSequenceTables(
  tables = file.path(NEONMICROBE_DIR_OUTPUTS(), "HARV/combined",
                     paste0("HARV_asv_combined_", unique_runs, ".Rds"))
)

seqtab_collapse_filename <- file.path(NEONMICROBE_DIR_OUTPUTS(), "HARV/combined",
                                      "NEON_16S_seqtab_HARV_COMBINED_COLLAPSED.Rds")
t0 <- Sys.time()
seqtab_collapse <- collapseNoMismatch(seqtab_joined)
saveRDS(seqtab_collapse, seqtab_collapse_filename)
t1 <- Sys.time()
t1 - t0 # Took 6.93 hours on socs-stats.ucsc.edu

seqtab_collapse <- readRDS(seqtab_collapse_filename)