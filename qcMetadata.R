qcMetadata <- function(metadata, outDir=NULL, pairedReads="Y", rmDupes=TRUE, rmFlagged="N", verbose=FALSE) {
  # library(plyr)
  options(stringsAsFactors = FALSE)
  
  # metadata <- readSequenceMetadata(metadata)
  
  # validate pairedReads
  if(!(pairedReads %in% c("Y", "N")) ) {
    stop("value for argument pairedReads invalid. Must be 'Y' or 'N'.")
  }
  
  # validate rmFlagged
  if(!(rmFlagged %in% c("Y", "N")) ) {
    stop("value for argument rmFlagged invalid. Must be 'Y' or 'N'.")
  }
  
  # get targetGene and confirm that only one targetGene is in input data set
  targetGene <- unique(metadata$targetGene)
  if(length(targetGene) > 1) {
    stop("more than one targetGene in input data set. Only one targetGene can be QCed at a time.")
  }
  
  # validate output folder for QCed metadata
  if(is.null(outDir)) {
    outDir <- file.path(NEONMICROBE_DIR_SEQMETA(), "qc_metadata")
  }
  if(!dir.exists(outDir) ) {
    message("Output directory does not exist")
    return(NULL)
    # dir.create(outDir, recursive=TRUE)
  }
  
  # Print size of dataset
  message(paste("Input dataset contains", nrow(metadata), "rows.") )
  
  # Remove flagged records, if rmFlagged="Y"
  if(rmFlagged=="Y") {
    message("Removing records with existing quality flag...")
    # Define flag values for removal based on LOV values #
    flagVals <- "Fail|legacyData"
    # Remove NA values #
    # metadata[is.na(metadata)] <- ""
    flagFields <- grep("qaqcStatus|dataQF", names(metadata))
    ind <- vector()
    for(i in flagFields) {
      flagged <- grep(flagVals, metadata[,i])
      ind <- c(ind,flagged)
    }
    if(length(ind)==0) {
      message("No flagged records found.")
    } else {
      ind <- unique(ind)
      numDupes <- length(ind)
      message(paste0(length(ind), " flagged records found. Removing flagged record(s).") )
      metadata <- metadata[-ind, ]
    }
    Sys.sleep(3)
  }
  
  # check for and remove duplicate sequence file names
  message("QC check for duplicate sequence file names...")
  # Add pause #
  Sys.sleep(1)
  dupeSeqIDs <- as.character(metadata$rawDataFileName[duplicated(metadata$rawDataFileName)] )
  if(length(dupeSeqIDs)==0) {
    message("QC check Pass. No duplicate sequence file names.")
  } else {
    numDupes <- length(dupeSeqIDs)
    message(paste0("QC check Fail. ", numDupes, " duplicate sequence file names found. Removing duplicated file(s).") )
    if(verbose) message(paste0("Removing duplicated row: ", which(duplicated(metadata$rawDataFileName)) ) )
    metadata <- metadata[!duplicated(metadata$rawDataFileName), ]
  }
  
  
  # check for and flag duplicate dnaSampleIDs
  message("QC checking duplicate dnaSampleIDs...")
  # Add pause #
  Sys.sleep(2)
  metadata$runDir <- ""
  metadata$runDir[grep("R1", metadata$rawDataFileDescription)] <- "R1"
  metadata$runDir[grep("R2", metadata$rawDataFileDescription)] <- "R2"
  dnaIDsPerRunDir <- paste(metadata$dnaSampleID, metadata$sequencerRunID, metadata$runDir, sep="-")
  metadata$duplicateDnaSampleIDFlag <- "0"
  if(any(duplicated(dnaIDsPerRunDir)) ) {
    ind <- which(duplicated(dnaIDsPerRunDir))
    metadata$duplicateDnaSampleIDFlag[ind] <- "1"
    if(rmDupes==TRUE) {
      metadata <- metadata[-ind,]
      message("Duplicated dnaSampleID record(s) found and removed.")
    } else {
      message("Duplicated dnaSampleID(s) found. Flagging affected record(s).")
    }
  } else {
    message("No duplicated dnaSampleID records found.")
  }
  
  # Subset qced data to un-flagged records
  metaFlagged <- metadata[metadata$duplicateDnaSampleIDFlag=="1", ]
  metaNotFlagged <- metadata[metadata$duplicateDnaSampleIDFlag=="0", ]
  
  metaNotFlagged$dnaSampleID <- as.character(metaNotFlagged$dnaSampleID)
  
  # Check existence of R1 (and R2 based on user input) #
  dnaSampTab <- data.frame(with(metaNotFlagged, table(dnaSampleID, runDir)) )
  
  # convert factors to characters
  dfType <- sapply(dnaSampTab, class)
  colsToFix <- names(dnaSampTab[which(dfType=='factor')])
  dnaSampTab[colsToFix] <- sapply(dnaSampTab[colsToFix], as.character)
  
  # Handle sequence data with missing run direction
  missingR1 <- dnaSampTab$dnaSampleID[which(dnaSampTab$Freq[dnaSampTab$runDir=="R1"]==0)]
  missingR2 <- dnaSampTab$dnaSampleID[which(dnaSampTab$Freq[dnaSampTab$runDir=="R2"]==0)]
  metaNotFlagged$runDirFlag <- "0"
  message("Check for missing forward or reverse read...")
  # Add pause #
  Sys.sleep(2)
  # Handle sequence data with missing R2 data
  if(length(missingR2)>0) {
    message(paste0("Reverse read missing from ", length(missingR2), " records") )
    # If specified, remove R1 file
    if(pairedReads=="Y") {
      message("Removing R1 files lacking a matching R2 file (default action when pairedReads='Y')" )
      metaNotFlagged <- metaNotFlagged[-intersect(which(metaNotFlagged$runDir=="R1"), which(metaNotFlagged$dnaSampleID %in% missingR2) ), ]
    } else {
      metaNotFlagged$runDirFlag[intersect(which(metaNotFlagged$runDir=="R1"), which(metaNotFlagged$dnaSampleID %in% missingR2) )] <- "1"
      message("Flagging R1 files lacking a matching R2 file (default action when pairedReads='N')" )
    }
  }
  # Handle sequence data with missing R1 data
  if(length(missingR1)>0) {
    message(paste0("Forward read missing from ", length(missingR1), " record(s). Removing R2 file for affected sequence data set(s).") )
    # remove R2 file
    metaNotFlagged <- metaNotFlagged[-intersect(which(metaNotFlagged$runDir=="R2"), which(metaNotFlagged$dnaSampleID %in% missingR1) ), ]
  }
  
  # Recombine original flagged records and remaining records post-initial flagging.
  out <- suppressMessages(plyr::join(metaNotFlagged, metaFlagged))
  
  write.csv(out, file.path(outDir, paste0("mmg_metadata_", gsub("\\ ", "", targetGene), "_QCd_", gsub("-", "", Sys.Date()), '.csv')), row.names = FALSE)
  cat(paste("Output QCed file contains", nrow(out), "rows. File saved to the following directory:", outDir))
  cat("\nNOTE: Always review output before proceeding with analysis.")
  return(out)
}


readSequenceMetadata <- function(metadata) {
  metadata_load_err <- FALSE
  if(class(metadata) == "data.frame") return(metadata)
  if(class(metadata) == "character") {
    if(file.exists(metadata)) {
      return(read.csv(metadata))
    }
  }
  stop("'metadata' must be the data.frame output from downloadSequenceMetadata() or ",
       "the filepath to a local csv copy of the output from downloadSequenceMetadata()")
}



downloadRawSequenceData <- function(metadata, outDir = NULL, overwrite=FALSE,
                                    ignore_tar_files=TRUE, checkSize=TRUE, verbose=FALSE) {
  
  # library(utils)
  options(stringsAsFactors = FALSE)
  
  metadata <- readSequenceMetadata(metadata)
  
  if(is.null(outDir)) {
    outDir <- NEONMICROBE_DIR_SEQUENCE()
  }
  
  if(!dir.exists(outDir)) {
    warning("Specified output directory does not exist.")
    return(invisible(NULL))
  }
  
  if(ignore_tar_files) {
    tar_ind <- grep('\\.tar\\.gz', metadata$rawDataFileName)
    if(length(tar_ind) > 0) {
      metadata <- metadata[-tar_ind, ]
      if(nrow(metadata) > 0) {
        message(length(tar_ind), " row(s) in metadata associated with batch-level sequence data
                were ignored prior to downloading raw sequence data.")
      } else {
        stop("No rows remain in metadata after removing rows associated with batch-level sequence
             data. Consider setting ignore_tar_files = FALSE")
      }
    }
  }
  
  # Get unique file names
  metadata.u <- metadata[!duplicated(metadata$rawDataFilePath), ]
  message(paste("There are", nrow(metadata.u), "unique raw sequence files to download."))
  
  #Loop to check existence and cumulative size of files
  cat("checking file sizes...\n")
  fileSize <- 0
  n_goodfiles <- 0
  idxrem <- vector()
  for(i in 1:nrow(metadata.u)) {
    # get file metadata
    response <- httr::HEAD(metadata.u$rawDataFilePath[i])
    # check for file found
    if(is.null(httr::headers(response)[["Content-Length"]])) {
      cat(paste('No file found for url ', metadata.u$rawDataFilePath[i], '. Skipping\n', sep=''))
      idxrem <- c(idxrem, i)
    } else {
      # grab file size
      fileSize <- fileSize + as.numeric(httr::headers(response)[["Content-Length"]])
      n_goodfiles <- n_goodfiles + 1
    }
  }
  #  Sum up file sizes and convert bytes to MB
  totalFileSize <- fileSize/1e6
  
  # Remove missing files from download list
  if(length(idxrem)>0) {
    metadata.u <- metadata.u[-idxrem, ]
  }
  
  if(checkSize==TRUE) {
    resp <- readline(paste("Continuing will download",nrow(metadata.u), "files totaling approximately",
                           totalFileSize, "MB. Do you want to proceed? y/n: ", sep=" "))
    if(!(resp %in% c("y","Y"))) stop("Stopping")
  }else{
    cat("Downloading", n_goodfiles, "files totaling approximately", totalFileSize," MB.\n")
  }
  
  # Get target gene to know which directory within outDir to sort into
  targetGene <- metadata.u$targetGene
  targetGene <- sub("16S rRNA", "16S", targetGene)
  if(!all(targetGene %in% c("ITS", "16S"))) {
    warning("Target gene is missing from at least one fastq file. These files will not be sorted ",
            "into a subdirectory within the specified outDir.")
  }
  
  # Download files
  download_success <- list()
  
  progressbar <- txtProgressBar(min = 0, max = nrow(metadata.u), style = 3)
  
  for(i in 1:nrow(metadata.u)) {
    if(is.na(targetGene[i])) {
      destfile <- file.path(outDir, metadata.u$rawDataFileName[i])
    } else {
      destfile <- file.path(outDir, targetGene[i], metadata.u$rawDataFileName[i])
    }
    if(file.exists(destfile) & !identical(overwrite, TRUE)) {
      if(verbose) message(destfile, " already exists and 'overwrite' is FALSE. Skipping.")
      download_success[[i]] <- "already_existed"
      next
    } else {
      tryCatch({
        download_success[[i]] <- utils::download.file(
          url = metadata.u$rawDataFilePath[i],
          destfile = destfile,
          quiet = !verbose)
      }, error = function(e) { # Occasionally an error arises because _fastq should be replaced by .fastq
        tryCatch({
          revised_url <- sub("_fastq", ".fastq", as.character(metadata.u$rawDataFilePath[i]))
          download_success[[i]] <- download.file(
            url = revised_url,
            destfile = destfile,
            method="curl",
            quiet = !verbose)
        }, error = function(f) {
          message("Could not download from URL: ", metadata.u$rawDataFilePath[i])
          download_success[[i]] <- 2
        })
      })
    }
    if(verbose) message("Finished downloading ", destfile, ".\n")
    setTxtProgressBar(progressbar, i)
  }
  close(progressbar)
  
  message("Finished download raw sequence files to ", outDir)
  return(invisible(download_success))
}





