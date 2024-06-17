write_qiime2_manifest_file <- function(meta, forwrdPrimer, reversePrimer, runDir, Subdirectory){
  
  meta <- select(meta, c("dnaSampleID","rawDataFileName", "rawDataFilePath",  
                         "domainID", "siteID", "plotID", "internalLabID.seq",
                         "instrument_model", "illuminaAdapterKit", "illuminaIndex1", "illuminaIndex2",
                         "targetGene", "sampleMaterial", "nucleicAcidQuantMethod", "yr", "mnth"))
  meta <- unique(meta)
  meta$direction <- ifelse(grepl("_R1", meta$rawDataFileName), "forward", "reverse")
  
  meta <- meta %>% pivot_wider(id_cols = c("dnaSampleID", "domainID", "siteID", "plotID", "internalLabID.seq",
                               "instrument_model", "illuminaAdapterKit", "illuminaIndex1", "illuminaIndex2",
                               "targetGene", "sampleMaterial", "nucleicAcidQuantMethod", "yr", "mnth"),
                               values_from = c("rawDataFileName", "rawDataFilePath"),
                               names_from = "direction")
  
  processedSeqFileNameLocation = stringr::str_locate(meta$rawDataFileName_forward,"_R1")
  processedSeqFileNameLocation <- processedSeqFileNameLocation[,1]
  processedSeqFileName = stringr::str_sub(meta$rawDataFileName_forward, 1, processedSeqFileNameLocation - 1)
  meta$processedSeqFileName <- paste(processedSeqFileName, "fastq", sep = ".")
  
  
  qiime2_manifest <- data.frame("sample-id" = meta$dnaSampleID,
                               "forward-absolute-filepath" = paste(Subdirectory, meta$siteID, meta$yr, meta$mnth, "16S", meta$rawDataFileName_forward, sep = "/"),
                               "reverse-absolute-filepath" = paste(Subdirectory, meta$siteID, meta$yr, meta$mnth, "16S", meta$rawDataFileName_reverse, sep = "/"),
                               "wget_forward-absolute-filepath" = meta$rawDataFilePath_forward,
                               "wget_reverse-absolute-filepath" = meta$rawDataFilePath_reverse,
                               Subdirectory = Subdirectory,
                               processedSeqFileName = meta$processedSeqFileName,
                               domainID = meta$domainID,
                               siteID = meta$siteID,
                               plotID = meta$plotID,
                               internalLabID.seq = meta$internalLabID.seq,
                               instrument_model = meta$instrument_model,
                               illuminaAdapterKit = meta$illuminaAdapterKit,
                               IlluminaIndex1 = meta$illuminaIndex1,
                               IlluminaIndex2 = meta$illuminaIndex2,
                               targetGene = meta$targetGene,
                               sampleMaterial = meta$sampleMaterial,
                               nucleicAcidQuantMethod = meta$nucleicAcidQuantMethod,
                               forwardPrimer = forwrdPrimer,
                               reversePrimer = reversePrimer,
                               runDir = runDir
                               ) 
  
  return (qiime2_manifest)
}
