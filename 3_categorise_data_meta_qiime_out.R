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

source("write_qiime2_manifest.R")
source("load_primers.R")


### Read by primers
meta_16_qc_primer1 <- read.csv("outputs/meta_primer1_all.csv")
meta_16_qc_primer2 <- read.csv("outputs/meta_primer2_all.csv")


### Run for each primer 
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