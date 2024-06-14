#' Download Sequence Metadata
#'
#' Loads soil marker gene sequencing metadata for specified target gene, site(s) and date(s),
#' with an option to download output by providing a valid output directory. This function uses
#' \code{\link[neonUtilities]{loadByProduct}} to conduct the downloads.
#'
#' Function by Lee F. Stanish and Clara Qin (2020). Currently available for testing only.
#'
#' @param sites Either the string 'all', meaning all available sites, or a character vector of 4-letter NEON site codes, e.g. c('ONAQ','RMNP'). Defaults to PRESET_SITES parameter in params.R.
#' @param startYrMo,endYrMo Either NA (default), meaning all available dates, or a character vector in the form YYYY-MM, e.g. 2017-01.
#' @param targetGene '16S' or 'ITS'.
#' @param sequencingRuns Either the string 'all', meaning all available sequencing runs, or a character vector of NEON sequencing run IDs, e.g. c('C25G9', 'B69PP').
#' @param dpID NEON data product of interest. Default is soil marker gene sequences, and currently code only works for marker genes data products.
#' @param outDir Directory where a copy of the downloaded data will be saved. By default (NULL), this is file.path(NEONMICROBE_DIR_SEQMETA(), "raw_metadata"). If no copy should be saved, set outDir=FALSE.
#'
#' @return Data frame containing joined records from across the NEON soil marker gene sequence metadata, subsetted according to function arguments.
#' @export
#'
#' @examples
#' \dontrun {
#' meta_16s <- downloadSequenceMetadata(
#'   startYrMo = "2017-07", endYrMo = "2017-07",
#'   sites = c("KONZ", "CPER", "NOGP"),
#'   targetGene = "16S"
#' )
#' }
downloadSequenceMetadata <- function(sites='all', startYrMo=NA, endYrMo=NA, targetGene= "all",
                                     sequencingRuns = "", dpID = "DP1.10108.001", outDir = NULL,
                                     include.provisional = FALSE) {
  # author: Lee Stanish
  # date: 2020-08-13
  # function loads soil marker gene sequencing metadata for target gene, site(s) and date(s)
  # option to download output by providing a valid output directory
  # sites: character vector of valid site ID's, or 'all' for all sites
  # targetGene: '16S',  'ITS', 'all'
  # startYrMo: start date, format YYYY-MM
  # endYrMo: end date, format YYYY-MM
  # dpID: NEON data product of interest. Default is soil marker gene sequences, and currently code only works for this dpID
  # outDir: directory for outputs. Defaults to output directory in parameters file
  
  # library(neonUtilities)
  # library(plyr)
  # library(dplyr)
  
  # check valid data values entered
  ## validate dpID ##
  if(!grepl("DP1", dpID) | !grepl('\\.001', dpID) | !grepl('10108|20280|20282', dpID)) {
    message("Invalid Data Product ID: must follow convention 'DP1.[5-digit value].001' and must be a marker genes data product ID")
    return(NULL)
  } else {
    dpID <- dpID
  }
  
  # validate target gene
  if(!grepl("16S|ITS|all", targetGene)) {
    message("Invalid targetGene: must be either '16S', 'ITS', 'all' ")
    return(NULL)
  } else {
    targetGene <- targetGene
  }
  
  # validate site(s)
  terrSiteList <- c("all","HARV","SCBI","OSBS","GUAN","UNDE","KONZ","ORNL","TALL","WOOD","CPER","CLBJ","YELL","NIWO",
                    "SRER","ONAQ","WREF","SJER","TOOL","BONA","PUUM","BART","BLAN","SERC","SCBI","DSNY","JERC","LAJA",
                    "TREE","STEI","KONA","UKFS","MLBS","GRSM","LENO","DELA","NOGP","DCFS","STER","RMNP","OAES","MOAB",
                    "JORN","ABBY","TEAK","SOAP","BARR","DEJU","HEAL")
  if(!any(sites %in% terrSiteList)){
    message("Invalid site(s): must be a valid NEON site or 'all'")
    return(NULL)
  } else {
    sites <- sites
  }
  
  # validate output directory
  if(is.null(outDir)) {
    outDir <- file.path(NEONMICROBE_DIR_SEQMETA(), "raw_metadata")
  }
  if(!identical(outDir, FALSE)) {
    if(!dir.exists(outDir) ) {
      message("Output directory does not exist")
      return(NULL)
    }
  }
  
  
  message("loading metadata...")
  mmgL1 <- neonUtilities::loadByProduct(dpID, sites, package = 'expanded', check.size = F, startdate = startYrMo, enddate = endYrMo, include.provisional = include.provisional) # output is a list of each metadata file
  
  
  # for target data product and targetGene: extract lists into data.frames
  if(grepl("10108", dpID)) {
    seq16S <- mmgL1$mmg_soilMarkerGeneSequencing_16S
    seq16S$targetGene <-"16S rRNA"
    seqITS <- mmgL1$mmg_soilMarkerGeneSequencing_ITS
    seqITS$targetGene <- "ITS"
    pcr16S <- mmgL1$mmg_soilPcrAmplification_16S
    pcrITS <- mmgL1$mmg_soilPcrAmplification_ITS
    raw <- mmgL1$mmg_soilRawDataFiles
    dna <- mmgL1$mmg_soilDnaExtraction
    seq <- rbind(seq16S, seqITS)
    pcr <- rbind(pcr16S, pcrITS)
    varfile <- mmgL1$variables_10108
    
    if(targetGene=="16S") {
      message("filtering to 16S data")
      seq <- seq16S
      pcr <- pcr16S
    }
    if(targetGene=="ITS") {
      message("filtering to ITS data")
      seq <- seqITS
      pcr <- pcrITS
    }
  }
  
  if(grepl("20280", dpID)) {
    seq16S <- mmgL1$mmg_benthicMarkerGeneSequencing_16S
    seq16S$targetGene <-"16S rRNA"
    seqITS <- mmgL1$mmg_benthicMarkerGeneSequencing_ITS
    seqITS$targetGene <- "ITS"
    pcr16S <- mmgL1$mmg_benthicPcrAmplification_16S
    pcrITS <- mmgL1$mmg_benthicPcrAmplification_ITS
    raw <- mmgL1$mmg_benthicRawDataFiles
    dna <- mmgL1$mmg_benthicDnaExtraction
    seq <- rbind(seq16S, seqITS)
    pcr <- rbind(pcr16S, pcrITS)
    varfile <- mmgL1$variables_20280
    
    if(targetGene=="16S") {
      message("filtering to 16S data")
      seq <- seq16S
      pcr <- pcr16S
    }
    if(targetGene=="ITS") {
      message("filtering to ITS data")
      seq <- seqITS
      pcr <- pcrITS
    }
  }
  
  if(grepl("20282", dpID)) {
    seq16S <- mmgL1$mmg_swMarkerGeneSequencing_16S
    seq16S$targetGene <-"16S rRNA"
    seqITS <- mmgL1$mmg_swMarkerGeneSequencing_ITS
    seqITS$targetGene <- "ITS"
    pcr16S <- mmgL1$mmg_swPcrAmplification_16S
    pcrITS <- mmgL1$mmg_swPcrAmplification_ITS
    raw <- mmgL1$mmg_swRawDataFiles
    dna <- mmgL1$mmg_swDnaExtraction
    seq <- rbind(seq16S, seqITS)
    pcr <- rbind(pcr16S, pcrITS)
    varfile <- mmgL1$variables_20282
    
    if(targetGene=="16S") {
      message("filtering to 16S data")
      seq <- seq16S
      pcr <- pcr16S
    }
    if(targetGene=="ITS") {
      message("filtering to ITS data")
      seq <- seqITS
      pcr <- pcrITS
    }
  }
  
  # remove unnecessary/redundant columns from tables
  raw <- dplyr::select(raw, -domainID, -siteID, -namedLocation, -laboratoryName, -sequencingFacilityID, -collectDate, -dnaSampleCode)
  dna <- dplyr::select(dna, -domainID, -siteID, -namedLocation, -laboratoryName, -collectDate)
  pcr <- dplyr::select(pcr, -domainID, -siteID, -namedLocation, -laboratoryName, -collectDate)
  
  # convert factors to characters (bug in output of loadByProduct)
  i <- sapply(seq, is.factor)
  seq[i] <- lapply(seq[i], as.character)
  j <- sapply(raw, is.factor)
  raw[j] <- lapply(raw[j], as.character)
  j <- sapply(dna, is.factor)
  dna[j] <- lapply(dna[j], as.character)
  
  
  # If specified, filter by sequencing run ID
  if(sequencingRuns[1] != "") {
    raw <- raw[which(raw$sequencerRunID %in% sequencingRuns), ]
    # Validate sequencing run ID argument
    if(nrow(raw) == 0) {
      warning("After filtering by specified sequencing run ID(s), no records remain. Double-check your sequencing run  ID(s).")
      return(NULL)
    }
  }
  
  # Join sequencing metadata with raw data files metadata
  if(targetGene=="16S") {
    if(any(grepl("ITS", raw$rawDataFileName))) {
      rawCleaned <- raw[-grep("ITS", raw$rawDataFileName), ]
    } else {
      rawCleaned <- raw
    }
    joinedTarget <- dplyr::left_join(rawCleaned, seq, by=c('dnaSampleID', 'sequencerRunID'))
    out <- joinedTarget[!is.na(joinedTarget$uid.y), ]
  }
  if(targetGene=="ITS") {
    if(any(grepl("16S", raw$rawDataFileName))) {
      rawCleaned <- raw[-grep("16S", raw$rawDataFileName), ]
    } else {
      rawCleaned <- raw
    }
    joinedTarget <- dplyr::left_join(rawCleaned, seq, by=c('dnaSampleID', 'sequencerRunID'))
    out <- joinedTarget[!is.na(joinedTarget$uid.y), ]
  }
  if(targetGene=="all") {
    joinedTarget <- dplyr::left_join(raw, seq, by=c('dnaSampleID', 'sequencerRunID'))
    out <- joinedTarget[!is.na(joinedTarget$uid.y), ]
    message(paste0(length(grep("16S", out$rawDataFileName)), " 16S records and ", length(grep("ITS", out$rawDataFileName)), " ITS records found."))
  }
  
  # clean up redundant column names
  names(out) <- gsub("\\.x", ".rawFiles", names(out))
  names(out) <- gsub("\\.y", ".seq", names(out))
  
  # join with DNA extraction metadata
  outDNA <- dplyr::left_join(out, dna, by=c('plotID', 'dnaSampleID'))
  # clean up redundant column names
  names(outDNA) <- gsub("\\.x", ".seq", names(outDNA))
  names(outDNA) <- gsub("\\.y", ".dna", names(outDNA))
  names(outDNA)[names(outDNA)=="uid"] <- 'uid.dna'
  names(outDNA)[names(outDNA)=="remarks"] <- 'remarks.dna'
  names(outDNA)[names(outDNA)=="dataQF"] <- 'dataQF.dna'
  names(outDNA)[names(outDNA)=="processedBy"] <- 'processedBy.seq'
  names(outDNA)[names(outDNA)=="processedDate"] <- 'processedDate.dna'
  names(outDNA)[names(outDNA)=="publicationDate"] <- 'publicationDate.dna'
  names(outDNA)[names(outDNA)=="dnaProcessedBy"] <- 'processedBy.dna'
  
  # join with PCR amplification metadata
  outPCR <- dplyr::left_join(outDNA, pcr, by=c('plotID', 'dnaSampleID', 'targetGene'))
  names(outPCR)[names(outPCR)=="uid"] <- "uid.pcr"
  names(outPCR)[names(outPCR)=="processedDate"] <- "processedDate.pcr"
  names(outPCR)[names(outPCR)=="testProtocolVersion"] <- "testProtocolVersion.pcr"
  names(outPCR)[names(outPCR)=="qaqcStatus"] <- "qaqcStatus.pcr"
  names(outPCR)[names(outPCR)=="processedBy"] <- "processedBy.pcr"
  names(outPCR)[names(outPCR)=="remarks"] <- "remarks.pcr"
  names(outPCR)[names(outPCR)=="dataQF"] <- "dataQF.pcr"
  names(outPCR)[names(outPCR)=="publicationDate"] <- "publicationDate.pcr"
  names(outPCR)[names(outPCR)=="internalLabID.y"] <- "internalLabID.pcr"
  
  # download local copy to output dir path
  # unless user provides outDir=FALSE
  if(!identical(outDir, FALSE)) {
    if(targetGene != "all") {
      write.csv(outPCR, paste0(outDir, "/mmg_soilMetadata_", targetGene, "_", sub(" ", "_", gsub(":", "", Sys.time())), ".csv"),
                row.names=FALSE)
    } else {
      out16S <- outPCR[grep("16S", outPCR$targetGene), ]
      outITS <- outPCR[grep("ITS", outPCR$targetGene), ]
      write.csv(out16S, paste0(outDir, "/mmg_soilMetadata_16S_", sub(" ", "_", gsub(":", "", Sys.time())), ".csv"),
                row.names=FALSE)
      write.csv(outITS, paste0(outDir, "/mmg_soilMetadata_ITS_", sub(" ", "_", gsub(":", "", Sys.time())), ".csv"),
                row.names=FALSE)
    }
    message(paste0("metadata downloaded to: ", outDir) )
    
    # download variables file (required for zipsByUri)
    write.csv(varfile, paste0(outDir, "/mmg_variables.csv") )
    message(paste0("variables file downloaded to: ", outDir) )
  }
  
  
  return(outPCR)
}
