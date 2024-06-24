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
source("primer_2_grep.R")
source("write_qiime2_manifest.R")
source("qcMetadata.R")
### NeonMicrobe Vignette for reference: https://people.ucsc.edu/~claraqin/analyze-neon-greatplains-16s.R

### DATA from Laura filter out to get all the relevant sites
datasets <- read.csv("data/MRR_2024_dataset_descript.csv")
datasets <- datasets %>%
  filter(MRR_use == "Y")


### Following NeonMicrobe download data instructions (will find link but this is their download data vignette)
setBaseDirectory(dir=getwd())
# makeDataDirectories()
sites_to_DL <- datasets$Site
sites_to_DL <- sites_to_DL[!(sites_to_DL %in% c("HARV"))]
sites_to_DL <- "HARV" ## Test one site for now
sites_to_DL <- c("GUAN", "HARV", "OSBS", "UNDE")
meta_16s <- downloadSequenceMetadata(sites = sites_to_DL, 
                                     startYrMo = "2016-01", endYrMo = "2021-12", 
                                     targetGene = "16S", outDir = "data/sequenceMeta",
                                     include.provisional = TRUE)

meta_16s <- downloadSequenceMetadata(sites = sites_to_DL, 
                                     startYrMo = "2019-01", endYrMo = "2019-12", 
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

## This didn't work wo I wrote all the file names in toDL.csv and wrote a script to download all files
## This script is "data/raw_sequence/downloadAll.sh 
## run from terminal and make sure files go into data/raw_sequence/16S when finished
downloadRawSequenceData(meta_16s_qc, outDir = ".", overwrite = TRUE, verbose = TRUE)

meta_16s_qc$date_ymd <- as.Date(format(as.Date(meta_16s_qc$collectDate, format="%Y-%m-%d %H:%M:%S"), "%Y-%m-%d"))
meta_16s_qc_toDL <- meta_16s_qc %>%
  filter(date_ymd >= ymd("2019-01-01") & date_ymd < ymd("2020-01-01")) %>%
  select(rawDataFilePath)

write.csv(meta_16s_qc_toDL, quote = FALSE, row.names = FALSE,  "data/raw_sequence/16S/toDL.csv")


## Now use the Harvard sites only - metadata input file might need to be adjusted for other sites
meta_16_qc_HARV <- read.csv(file = "data/sequenceMeta/mmg_metadata_16SrRNA_QCd_20240614.csv")
meta_16_qc <- read.csv(file = "data/sequenceMeta/mmg_metadata_16SrRNA_QCd_20240616.csv")
meta_16_qc <- meta_16_qc %>%
  filter(siteID %in% sites_to_DL)







meta_16_qc_HARV <- filter(meta_16_qc, siteID == "HARV")

# meta_16_qc <- meta_16s_qc
meta_16_qc_HARV$date_ymd <- as.Date(format(as.Date(meta_16_qc_HARV$collectDate, format="%Y-%m-%d %H:%M:%S"), "%Y-%m-%d"))
meta_16_qc_HARV$yr <- year(meta_16_qc_HARV$date_ymd)
meta_16_qc_HARV$mnth <- month(meta_16_qc_HARV$date_ymd)
meta_16_qc_HARV$sampleCodeID <- paste(meta_16_qc_HARV$siteID, meta_16_qc_HARV$yr, meta_16_qc_HARV$mnth, sep="_") 



meta_16_qc$date_ymd <- as.Date(format(as.Date(meta_16_qc$collectDate, format="%Y-%m-%d %H:%M:%S"), "%Y-%m-%d"))
meta_16_qc$yr <- year(meta_16_qc$date_ymd)
meta_16_qc$mnth <- month(meta_16_qc$date_ymd)
meta_16_qc$sampleCodeID <- paste(meta_16_qc$siteID, meta_16_qc$yr, meta_16_qc$mnth, sep="_") 

meta_16_qc$processed_date_ymd <- as.Date(format(as.Date(meta_16_qc$processedDate.seq, format="%Y-%m-%d"), "%Y-%m-%d"))
meta_16_qc$processed_yr <- year(meta_16_qc$processed_date_ymd)
meta_16_qc$processed_mnth <- month(meta_16_qc$processed_date_ymd)
# meta_16_qc$sampleCodeID <- paste(meta_16_qc$siteID, meta_16_qc$yr, meta_16_qc$mnth, sep="_") 



test_sum <- meta_16_qc |>
  select(processed_yr, processed_mnth, laboratoryName, sequencingFacilityID, yr, mnth, domainID, siteID, testProtocolVersion.seq) |>
  dplyr::summarise(n = dplyr::n(), .by = c(processed_yr, processed_mnth, laboratoryName, sequencingFacilityID, yr, mnth, domainID, siteID, testProtocolVersion.seq)) |>
  dplyr::arrange(processed_yr)

test_sum$testProtocolVersion.seq <- as.factor(test_sum$testProtocolVersion.seq) 
test_sum$testProtocolVersion.seq <- factor(test_sum$testProtocolVersion.seq, levels = c("BMI_markerGeneSequencingSOP_v1", 
                                                                 "BMI_markerGeneSequencingSOP_v2", 
                                                                 "BMI_markerGenes_16Sv4v5SOP_v1.1", 
                                                                 "BMI_markerGenes_16Sv4v5SOP_v1.2", 
                                                                 "GMCF_markerGenes_16SSOP_v1a"))
test_sum$testProtocolVersion.seq_num <- as.numeric(test_sum$testProtocolVersion.seq)

write.csv(test_sum, file= "dplr_protocol_version_by_site_by_date.csv")


View(test_sum)

meta_16s_qc_toDL <- meta_16_qc %>%
  filter(date_ymd >= ymd("2017-01-01") & date_ymd < ymd("2018-01-01")) %>%
  select(rawDataFilePath)

write.csv(meta_16s_qc_toDL, quote = FALSE, row.names = FALSE,  "data/raw_sequence/16S/toDL.csv")



dirs_to_create <- paste(meta_16_qc$siteID, meta_16_qc$yr, meta_16_qc$mnth, "16S", sep="/") 
write.csv(dirs_to_create, "dirs_to_create.csv", row.names = FALSE, quote = FALSE, col.names =  FALSE)

primer_16S_fwd_primer1 = "CCTACGGGNBGCASCAG"
primer_16S_fwd_primer1_grep = primer_to_grep(primer_16S_fwd_primer1)
primer_16S_rev_primer1 = "GACTACNVGGGTATCTAATCC"
primer_16S_rev_primer1_grep = primer_to_grep(primer_16S_rev_primer1)

primer_16S_fwd_primer2 = "GTGYCAGCMGCCGCGGTAA"
primer_16S_fwd_primer2_grep = primer_to_grep(primer_16S_fwd_primer2)
primer_16S_rev_primer2 = "CCGYCAATTYMTTTRAGTTT"
primer_16S_rev_primer2_grep = primer_to_grep(primer_16S_rev_primer2)



# unique_runs <- unique(meta_16_qc_HARV$uid.rawFiles[meta_16_qc_HARV$date_ymd >= ymd("2018-01-01") &
#                                                      meta_16_qc_HARV$date_ymd <= ymd("2019-12-31")])
# unique_runs <- unique(meta_16_qc_HARV$uid.rawFiles[meta_16_qc_HARV$sampleCodeID %in% c("HARV_2018_10", "HARV_2019_7", "HARV_2020_7")])
unique_runs <- unique(meta_16_qc$uid.rawFiles[meta_16_qc$date_ymd >= ymd("2019-01-01") &
                                                  meta_16_qc$date_ymd <= ymd("2019-12-31")])

meta_16s_qc_toDL <- meta_16_qc %>%
  filter(date_ymd >= ymd("2018-01-01") & date_ymd < ymd("2019-01-01")) %>%
  select(rawDataFilePath)
unique_runs <- unique(meta_16_qc$uid.rawFiles[meta_16_qc$rawDataFilePath %in% meta_16s_qc_toDL$rawDataFilePath])
#
# unique_runs <- unique(meta_16_qc$sequencerRunID)
# 
# nrow(meta_16_qc %>%
#   dplyr::mutate(year = year(date_ymd)) %>%
#   filter(year(date_ymd) %in% c(2019)))
# 
# %>%
#   select(date_ymd, dnaSampleID, year, siteID, sequencerRunID) %>%
#   group_by(date_ymd,year, siteID, sequencerRunID)%>%
#   dplyr::count(dnaSampleID)

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


# harv_test <-meta_16_qc[meta_16_qc$uid.rawFiles %in% unique_runs,]

harv_test <-meta_16_qc_HARV[meta_16_qc_HARV$sampleCodeID %in% c("HARV_2018_10", "HARV_2019_7", "HARV_2020_7"),]
length(primers_2019)
harv_test$primers <- primers_2019
harv_test$direction <- ifelse(grepl("_R1", harv_test$rawDataFileName), "forward", "reverse")
harv_test_primer1 <- harv_test %>%
  filter(harv_test$primer == "PRIMER1")
harv_test_sum <- harv_test |>
  select(sampleCodeID, yr, mnth, siteID, direction, primers, testProtocolVersion.seq) |>
  dplyr::summarise(n = dplyr::n(), .by = c(sampleCodeID, siteID, yr, mnth, direction, primers, testProtocolVersion.seq))

write.csv()
harv_test_sum_filename <- harv_test$rawDataFileName[harv_test$primers == "ERROR"]

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
meta_16_qc_HARV_primer2 <- rbind(meta_16_qc_HARV_primer2, meta_16_qc_HARV_2019_primer2)

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

meta_16_qc$primer_cat <- ifelse(meta_16_qc$testProtocolVersion.seq %in% c("BMI_markerGeneSequencingSOP_v2", "BMI_markerGeneSequencingSOP_v1"), "PRIMER1", "PRIMER2")
meta_16_qc$primer_cat[meta_16_qc$uid.rawFiles %in% unique_runs] <- primers_2019

file_list_to_remove <- list.files("data/raw_sequence/16S/", pattern = ".fastq.gz")
toremove <- which((file_list_to_remove %in% meta_16_qc$rawDataFileName))
file.remove(paste0("data/raw_sequence/16S/", file_list_to_remove[toremove]))

meta_16_qc_primer1 <- meta_16_qc %>%
  filter(meta_16_qc$primer_cat == "PRIMER1")
meta_16_qc_primer2 <- meta_16_qc %>%
  filter(meta_16_qc$primer_cat == "PRIMER2")

# meta_16_qc_2019_primer1 <- meta_16_qc %>%
#   filter(sampleCodeID %in% unique_runs[primers_2019 == "PRIMER1"])
# meta_16_qc_2019_primer2 <-meta_16_qc %>%
#   filter(sampleCodeID %in% unique_runs[primers_2019 == "PRIMER2"])

# meta_16_qc_primer1 <- rbind(meta_16_qc_primer1, meta_16_qc_2019_primer1)
meta_16_qc_primer1 <- unique(meta_16_qc_primer1)
# meta_16_qc_primer2 <- rbind(meta_16_qc_primer2, meta_16_qc_2019_primer2)
meta_16_qc_primer2 <- unique(meta_16_qc_primer2)



qiime2_primer1 <- write_qiime2_manifest_file(meta_16_qc_primer1,
                                                  forwrdPrimer = primer_16S_fwd_primer1,
                                                  reversePrimer = primer_16S_rev_primer1,
                                                  runDir = "R1",
                                                  Subdirectory = "/home/laura_s/project_guy/2024_MRR")
qiime2_primer2 <- write_qiime2_manifest_file(meta_16_qc_primer2,
                                                  forwrdPrimer = primer_16S_fwd_primer2,
                                                  reversePrimer = primer_16S_rev_primer2,
                                                  runDir = "R1",
                                                  Subdirectory = "/home/laura_s/project_guy/2024_MRR")

write.csv(meta_16_qc_primer1, "outputs/meta_primer1_all.csv", row.names = FALSE)
meta_16_qc_primer1 <- read.csv("outputs/meta_primer1_all.csv")
write.csv(meta_16_qc_primer2, "outputs/meta_primer2_all.csv", row.names = FALSE)
meta_16_qc_primer2 <- read.csv("outputs/meta_primer2_all.csv")
write.csv(qiime2_primer1, "outputs/qiime2_manifest_primer1_all.csv", row.names = FALSE)
write.csv(qiime2_primer2, "outputs/qiime2_manifest_primer2_all.csv", row.names = FALSE)

meta_16_qc_primer1 <- meta_16_qc_primer1 %>%
  filter(meta_16_qc_primer1$siteID %in% sites_to_DL) %>%
  filter(meta_16_qc_primer1$yr < 2018)


unique_runs <- unique(meta_16_qc_primer1$sampleCodeID)
for (run_ID in unique_runs[1:length(unique_runs)]){
  meta_sub <- meta_16_qc_primer1[meta_16_qc_primer1$sampleCodeID == run_ID,]
  qiime2_primer1_sub <- write_qiime2_manifest_file(meta_sub,
                                                   forwrdPrimer = primer_16S_fwd_primer1,
                                                   reversePrimer = primer_16S_rev_primer1,
                                                   runDir = "R1",
                                                   Subdirectory = "/home/laura_s/project_guy/2024_MRR")
  qiime2_primer1_sub_short_out <- write_qiime2_manifest_file_short(meta_sub, Subdirectory = "/home/laura_s/project_guy/2024_MRR")
  qiime2_primer1_sub_short <- qiime2_primer1_sub_short_out$qiime2_manifest
  qiime2_primer1_sub_wget <- qiime2_primer1_sub_short_out$qiime2_wget
  qiime2_primer1_sub_errors <- qiime2_primer1_sub_short_out$errors
  write.csv(meta_sub, paste0("outputs/primer1/meta/", run_ID, "_meta_primer1_all.csv"), row.names = FALSE)
  write.csv(qiime2_primer1_sub, paste0("outputs/primer1/qiime_full/", run_ID, "_qiime2_manifest_primer1.csv"), row.names = FALSE)
  write.csv(qiime2_primer1_sub_short, paste0("outputs/primer1/qiime_simple/", run_ID, "_qiime2_manifest_primer1_short.csv"), row.names = FALSE)
  write.csv(qiime2_primer1_sub_wget, paste0("outputs/primer1/qiime_wget/", run_ID, "_qiime2_manifest_primer1_wget.csv"), row.names = FALSE)
  if(nrow(qiime2_primer1_sub_errors) > 0){
    write.csv(qiime2_primer1_sub_errors, paste0("outputs/primer1/qiime_simple/", run_ID, "_error_reads_qiime2_manifest_primer1_short.csv"), row.names = FALSE)
  }
}



unique_runs <- unique(meta_16_qc_primer2$sampleCodeID)
for (run_ID in unique_runs){
  meta_sub <- meta_16_qc_primer2[meta_16_qc_primer2$sampleCodeID == run_ID,]
  qiime2_primer2_sub <- write_qiime2_manifest_file(meta_sub,
                                                   forwrdPrimer = primer_16S_fwd_primer2,
                                                   reversePrimer = primer_16S_rev_primer2,
                                                   runDir = "R1",
                                                   Subdirectory = "/home/laura_s/project_guy/2024_MRR")
  qiime2_primer2_sub_short_out <- write_qiime2_manifest_file_short(meta_sub, Subdirectory = "/home/laura_s/project_guy/2024_MRR")
  qiime2_primer2_sub_short <- qiime2_primer2_sub_short_out$qiime2_manifest
  qiime2_primer2_sub_wget <- qiime2_primer2_sub_short_out$qiime2_wget
  qiime2_primer2_sub_errors <- qiime2_primer2_sub_short_out$errors
  write.csv(meta_sub, paste0("outputs/primer2/meta/", run_ID, "_meta_primer2_all.csv"), row.names = FALSE)
  write.csv(qiime2_primer2_sub, paste0("outputs/primer2/qiime_full/", run_ID, "_qiime2_manifest_primer2.csv"), row.names = FALSE)
  write.csv(qiime2_primer2_sub_short, paste0("outputs/primer2/qiime_simple/", run_ID, "_qiime2_manifest_primer2_short.csv"), row.names = FALSE)
  write.csv(qiime2_primer2_sub_wget, paste0("outputs/primer2/qiime_wget/", run_ID, "_qiime2_manifest_primer2_wget.csv"), row.names = FALSE)
  if(nrow(qiime2_primer1_sub_errors) > 0){
    write.csv(qiime2_primer2_sub_errors, paste0("outputs/primer2/qiime_simple/", run_ID, "_error_reads_qiime2_manifest_primer2_short.csv"), row.names = FALSE)
  }
}




# unique_runs <- unique(meta_16_qc$sequencerRunID)

# fl_nm
# fl_nm <- meta_16_qc$rawDataFileName[grepl("BMI_Plate6", meta_16_qc$rawDataFileName, ignore.case = TRUE)]


####### FOR PrE 2019 DATASETS ---
# from: https://data.neonscience.org/documents/-/document_library_display/kV4WWrbEEM2s/view_file/3370805?_110_INSTANCE_kV4WWrbEEM2s_redirect=https%3A%2F%2Fdata.neonscience.org%2Fdocuments%2F-%2Fdocument_library_display%2FkV4WWrbEEM2s%2Fview%2F2237401%3F_110_INSTANCE_kV4WWrbEEM2s_keywords%3D%26_110_INSTANCE_kV4WWrbEEM2s_topLink%3Dhome%26_110_INSTANCE_kV4WWrbEEM2s_advancedSearch%3Dfalse%26_110_INSTANCE_kV4WWrbEEM2s_delta2%3D20%26_110_INSTANCE_kV4WWrbEEM2s_cur2%3D2%26p_r_p_564233524_resetCur%3Dfalse%26_110_INSTANCE_kV4WWrbEEM2s_andOperator%3Dtrue%26_110_INSTANCE_kV4WWrbEEM2s_delta1%3D20
# Table 1
# Accessed June 14 2024

files_trimmed <- list.files("outputs/mid_process/16S/1_trimmed/primer1/")
files_to_remove <- which(!(files_trimmed %in% meta_16_qc_HARV_primer1$rawDataFileName))
files_to_trim <- which(!(meta_16_qc_HARV_primer1$rawDataFileName %in% files_trimmed))

file.remove(paste0("outputs/mid_process/16S/1_trimmed/primer1/", files_trimmed[files_to_remove]))


View(meta_16_qc_HARV_primer1[files_to_trim,])

primer_16S_fwd_primer1 = "CCTACGGGNBGCASCAG"
primer_16S_rev_primer1 = "GGACTACNVGGGTATCTAATCC"
trim_trackReads <- trimPrimers16S(
  meta_16_qc_HARV_primer1$rawDataFileName[files_to_trim], 
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

files_trimmed <- list.files("outputs/mid_process/16S/1_trimmed/primer2/")
files_to_remove <- which(!(files_trimmed %in% meta_16_qc_HARV_primer2$rawDataFileName))
files_to_trim <- which(!(meta_16_qc_HARV_primer2$rawDataFileName %in% files_trimmed))

file.remove(paste0("outputs/mid_process/16S/1_trimmed/primer2/", files_trimmed[files_to_remove]))


primer_16S_fwd_primer2 = "GTGYCAGCMGCCGCGGTAA"
primer_16S_rev_primer2 = "CCGYCAATTYMTTTRAGTTT"
trim_trackReads <- trimPrimers16S(
  meta_16_qc_HARV_primer2$rawDataFileName[files_to_trim], in_subdir = "raw", out_subdir = "1_trimmed/primer2", 
  meta = meta_16_qc_HARV_primer2, 
  multithread = FALSE,# set multithread = FALSE on Windows computers though
  primer_16S_fwd = primer_16S_fwd_primer2, 
  primer_16S_rev = primer_16S_rev_primer2
)

files_trimmed <- list.files("outputs/mid_process/16S/2_filtered/primer1/")
files_to_remove <- which(!(files_trimmed %in% meta_16_qc_HARV_primer1$rawDataFileName))
files_to_filter <- which(!(meta_16_qc_HARV_primer1$rawDataFileName %in% files_trimmed))

file.remove(paste0("outputs/mid_process/16S/2_filtered/primer1/", files_trimmed[files_to_remove]))


filter_trackReads <- qualityFilter16S(
  meta_16_qc_HARV_primer1$rawDataFileName[files_to_filter], in_subdir = "1_trimmed/primer1", 
  out_subdir = "2_filtered/primer1",
  meta = meta_16_qc_HARV_primer1, truncLen = 220, maxEE = 8, multithread = FALSE
)



files_trimmed <- list.files("outputs/mid_process/16S/2_filtered/primer2/")
files_to_remove <- which(!(files_trimmed %in% meta_16_qc_HARV_primer2$rawDataFileName))
files_to_filter <- which(!(meta_16_qc_HARV_primer2$rawDataFileName %in% files_trimmed))

file.remove(paste0("outputs/mid_process/16S/2_filtered/primer2/", files_trimmed[files_to_remove]))


filter_trackReads <- qualityFilter16S(
  meta_16_qc_HARV_primer2$rawDataFileName[files_to_filter], in_subdir = "1_trimmed/primer2", 
  out_subdir = "2_filtered/primer2",
  meta = meta_16_qc_HARV_primer2, truncLen = 220, maxEE = 8, multithread = FALSE
)



meta_16_qc_primer1 <- read.csv("outputs/meta_primer1_all.csv")
meta_16_qc_HARV_primer1 <- meta_16_qc_primer1 %>%
  filter(siteID == "HARV")
###### 

unique_runs <- unique(meta_16_qc_HARV_primer1$sampleCodeID)

for(i in 1:length(unique_runs)) {
  meta_thisrun <- meta_16_qc_HARV_primer1[which(meta_16_qc_HARV_primer1$sampleCodeID==unique_runs[i]),]
  fl_nm_thisrun <- meta_thisrun$rawDataFileName
  dada_out <- runDada16S(
    fl_nm_thisrun, in_subdir = "2_filtered/primer1/", meta = meta_16_qc_HARV_primer1,
    out_seqtab = file.path(NEONMICROBE_DIR_OUTPUTS(), "HARV/primer1",
                           paste0("HARV_asv_", unique_runs[i], ".Rds")),
    out_track = file.path(NEONMICROBE_DIR_OUTPUTS(), "HARV/primer1",
                          paste0("HARV_track_", unique_runs[i], ".csv")),
    verbose = FALSE,
    multithread = FALSE
  )
}

seqtab_joined <- mergeSequenceTables(
  tables = file.path(NEONMICROBE_DIR_OUTPUTS(), "HARV/primer1",
                     paste0("HARV_asv_", unique_runs, ".Rds"))
)

seqtab_collapse_filename <- file.path(NEONMICROBE_DIR_OUTPUTS(), "HARV/primer1",
                                      "NEON_16S_seqtab_HARV_COMBINED_COLLAPSED.Rds")
t0 <- Sys.time()
seqtab_collapse <- collapseNoMismatch(seqtab_joined)
saveRDS(seqtab_collapse, seqtab_collapse_filename)
t1 <- Sys.time()
t1 - t0 # Took 6.93 hours on socs-stats.ucsc.edu

seqtab_collapse <- readRDS(seqtab_collapse_filename)





meta_16_qc_primer2 <- read.csv("outputs/meta_primer2_all.csv")
meta_16_qc_HARV_primer2 <- meta_16_qc_primer2 %>%
  filter(siteID == "HARV")
unique_runs <- unique(meta_16_qc_HARV_primer2$sampleCodeID)

## Errors in K22DL
for(i in 1:length(unique_runs)) {
  meta_thisrun <- meta_16_qc_HARV_primer2[which(meta_16_qc_HARV_primer2$sampleCodeID==unique_runs[i]),]
  fl_nm_thisrun <- meta_thisrun$rawDataFileName
  dada_out <- runDada16S(
    fl_nm_thisrun, in_subdir = "2_filtered/primer2/", meta = meta_16_qc_HARV_primer2,
    out_seqtab = file.path(NEONMICROBE_DIR_OUTPUTS(), "HARV/primer2",
                           paste0("HARV_asv_", unique_runs[i], ".Rds")),
    out_track = file.path(NEONMICROBE_DIR_OUTPUTS(), "HARV/primer2",
                          paste0("HARV_track_", unique_runs[i], ".csv")),
    verbose = FALSE,
    multithread = FALSE
  )
}



seqtab_joined <- mergeSequenceTables(
  tables = file.path(NEONMICROBE_DIR_OUTPUTS(), "HARV/primer2",
                     paste0("HARV_asv_", unique_runs, ".Rds"))
)


seqtab_collapse_filename <- file.path(NEONMICROBE_DIR_OUTPUTS(), "HARV/primer2",
                                      "NEON_16S_seqtab_HARV_COMBINED_COLLAPSED.Rds")
t0 <- Sys.time()
seqtab_collapse <- collapseNoMismatch(seqtab_joined)
saveRDS(seqtab_collapse, seqtab_collapse_filename)
t1 <- Sys.time()
t1 - t0 # Took 6.93 hours on socs-stats.ucsc.edu

seqtab_collapse <- readRDS(seqtab_collapse_filename)



## currently here


