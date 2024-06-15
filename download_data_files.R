rm(list = ls())

library(neonUtilities)
library(tidyverse)
library(neonMicrobe)
library(ShortRead)
library(Biostrings)
library(dada2)
library(dplyr)
library(ggplot2)

source("load_meta_data_neonMicrobe.R")
source("primer_2_grep.R")
source("write_qiime2_manifest.R")

### NeonMicrobe Vignette for reference: https://people.ucsc.edu/~claraqin/analyze-neon-greatplains-16s.R

### DATA from Laura filter out to get all the relevant sites
datasets <- read.csv("data/MRR_2024_dataset_descript.csv")
datasets <- datasets %>%
  filter(MRR_use == "Y")


### Following NeonMicrobe download data instructions (will find link but this is their download data vignette)
setBaseDirectory(dir=getwd())
makeDataDirectories()
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

primer_16S_fwd_primer1 = "CCTACGGGNBGCASCAG"
primer_16S_fwd_primer1_grep = primer_to_grep(primer_16S_fwd_primer1)
primer_16S_rev_primer1 = "GGACTACNVGGGTATCTAATCC"
primer_16S_rev_primer1_grep = primer_to_grep(primer_16S_rev_primer1)
primer_16S_fwd_primer2 = "GTGYCAGCMGCCGCGGTAA"
primer_16S_fwd_primer2_grep = primer_to_grep(primer_16S_fwd_primer2)
primer_16S_rev_primer2 = "CCGYCAATTYMTTTRAGTTT"
primer_16S_rev_primer2_grep = primer_to_grep(primer_16S_rev_primer2)



unique_runs <- unique(meta_16_qc_HARV$sequencerRunID[meta_16_qc_HARV$date_ymd >= ymd("2019-01-01") &
                        meta_16_qc_HARV$date_ymd <= ymd("2019-12-31")])

primers_2019 <- unlist(lapply(unique_runs, function(seqRunID){
  meta_sub <- meta_16_qc_HARV[meta_16_qc_HARV$sequencerRunID == seqRunID,]
  file_location = paste0("data/raw_sequence/16S/",meta_sub$rawDataFileName[1])
  fastq <- readFastq(file_location)
  detected = FALSE
  i = 1
  outString = "ERROR"
  while(detected == FALSE && i <= length(fastq)){
    read <- as.character(fastq[i]@sread)
    if(stringr::str_detect(read, primer_16S_fwd_primer1_grep)){
      detected = TRUE
      outString = "PRIMER1"
    }
    else if(stringr::str_detect(read, primer_16S_fwd_primer2_grep)){
      detected = TRUE
      outString = "PRIMER2"
    }
    else if(stringr::str_detect(read, primer_16S_rev_primer1_grep)){
      detected = TRUE
      outString = "PRIMER1"
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

meta_16_qc_HARV_primer1 <- meta_16_qc_HARV %>%
  filter(meta_16_qc_HARV$date_ymd < ymd("2019-01-01")) %>%
  filter(date_ymd >= ymd("2016-01-01"))
meta_16_qc_HARV_primer2 <- meta_16_qc_HARV %>%
  filter(meta_16_qc_HARV$date_ymd >= ymd("2020-01-01"))

meta_16_qc_HARV_2019_primer1 <- meta_16_qc_HARV %>%
  filter(sequencerRunID %in% unique_runs[primers_2019 == "PRIMER1"])
meta_16_qc_HARV_2019_primer2 <-meta_16_qc_HARV %>%
  filter(sequencerRunID %in% unique_runs[primers_2019 == "PRIMER2"])

meta_16_qc_HARV_primer1 <- rbind(meta_16_qc_HARV_primer1, meta_16_qc_HARV_2019_primer1)
meta_16_qc_HARV_primer2 <- rbind(meta_16_qc_HARV_primer1, meta_16_qc_HARV_2019_primer2)

qiime2_HARV_primer1 <- write_qiime2_manifest_file(meta_16_qc_HARV_primer2,
                                                   forwrdPrimer = primer_16S_fwd_primer1,
                                                   reversePrimer = primer_16S_rev_primer1,
                                                   runDir = "R1",
                                                   Subdirectory = "/home/laura_s/project_guy/2022_microbiome/2022_qiime2/NEON_Qiime2_2023_Biodiversity_Conservation_NST/BART/16S")
qiime2_HARV_primer2 <- write_qiime2_manifest_file(meta_16_qc_HARV_primer2,
                                                    forwrdPrimer = primer_16S_fwd_primer2,
                                                    reversePrimer = primer_16S_rev_primer2,
                                                    runDir = "R1",
                                                    Subdirectory = "/home/laura_s/project_guy/2022_microbiome/2022_qiime2/NEON_Qiime2_2023_Biodiversity_Conservation_NST/BART/16S")

write.csv(meta_16_qc_HARV_primer1, "outputs/meta_primer1_HARV.csv")
write.csv(meta_16_qc_HARV_primer2, "outputs/meta_primer2_HARV.csv")
write.csv(qiime2_HARV_primer1, "outputs/qiime2_manifest_primer1_HARV.csv")
write.csv(qiime2_HARV_primer2, "outputs/qiime2_manifest_primer2_HARV.csv")


# unique_runs <- unique(meta_16_qc$sequencerRunID)

# fl_nm
# fl_nm <- meta_16_qc$rawDataFileName[grepl("BMI_Plate6", meta_16_qc$rawDataFileName, ignore.case = TRUE)]


####### FOR PrE 2019 DATASETS ---
# from: https://data.neonscience.org/documents/-/document_library_display/kV4WWrbEEM2s/view_file/3370805?_110_INSTANCE_kV4WWrbEEM2s_redirect=https%3A%2F%2Fdata.neonscience.org%2Fdocuments%2F-%2Fdocument_library_display%2FkV4WWrbEEM2s%2Fview%2F2237401%3F_110_INSTANCE_kV4WWrbEEM2s_keywords%3D%26_110_INSTANCE_kV4WWrbEEM2s_topLink%3Dhome%26_110_INSTANCE_kV4WWrbEEM2s_advancedSearch%3Dfalse%26_110_INSTANCE_kV4WWrbEEM2s_delta2%3D20%26_110_INSTANCE_kV4WWrbEEM2s_cur2%3D2%26p_r_p_564233524_resetCur%3Dfalse%26_110_INSTANCE_kV4WWrbEEM2s_andOperator%3Dtrue%26_110_INSTANCE_kV4WWrbEEM2s_delta1%3D20
# Table 1
# Accessed June 14 2024
primer_16S_fwd_primer1 = "CCTACGGGNBGCASCAG"
primer_16S_rev_primer1 = "GGACTACNVGGGTATCTAATCC"
trim_trackReads <- trimPrimers16S(
  meta_16_qc_HARV_primer1$rawDataFileName, 
  in_subdir = "raw", 
  out_subdir = "1_trimmed/primer1", meta = meta_16_qc_HARV_primer1, 
  multithread = FALSE,# set multithread = FALSE on Windows computers though
  primer_16S_fwd =primer_16S_fwd_primer1, 
  primer_16S_rev = primer_16S_rev_primer1
)

####### FOR POST 2019 DATASETS --- 
# from: 
# "Improved Bacterial 16S rRNA Gene (V4 and V4-5) and Fungal Internal T
# Transcribed Spacer Marker Gene Primers for Microbial Community Surveys
# Waletrs et al 2015
# (Table 1 rows 3 for FWD and 5 for REV)
# 515f Modiﬁed GTGYCAGCMGCCGCGGTAA Parada et al. (4)
# 926r CCGYCAATTYMTTTRAGTTT Parada et al. (4)
primer_16S_fwd_primer2 = "GTGYCAGCMGCCGCGGTAA"
primer_16S_rev_primer2 = "CCGYCAATTYMTTTRAGTTT"
trim_trackReads <- trimPrimers16S(
  meta_16_qc_HARV_primer2$rawDataFileName, in_subdir = "raw", out_subdir = "1_trimmed/primer2", 
  meta = meta_16_qc_HARV_primer2, 
  multithread = FALSE,# set multithread = FALSE on Windows computers though
  primer_16S_fwd = primer_16S_fwd_primer2, 
  primer_16S_rev = primer_16S_rev_primer2
)




filter_trackReads <- qualityFilter16S(
  meta_16_qc_HARV_primer1$rawDataFileName, in_subdir = "1_trimmed/primer1", 
  out_subdir = "2_filtered/primer1",
  meta = meta_16_qc_HARV_primer1, truncLen = 220, maxEE = 8, multithread = FALSE
)
filter_trackReads <- qualityFilter16S(
  meta_16_qc_HARV_primer2$rawDataFileName, in_subdir = "1_trimmed/primer2", 
  out_subdir = "2_filtered/primer2",
  meta = meta_16_qc_HARV_primer2, truncLen = 220, maxEE = 8, multithread = FALSE
)






###### 



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