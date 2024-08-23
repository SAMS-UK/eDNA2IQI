#' Remove unspecified or ambiguous bases from the sequence reads.
#'
#' @param folder A string. Location of raw data files.
#' @param raw_files A string. List of files from same location.
#'
#' @return A list.

step1.2_remove_unspecified <- function(folder, raw_files) {
#  browser()

  messageColour("Function starting: 'remove_unspecified'
           Removing unspecified bases from sequence reads \n\n", "message")

  # Set file to store details of cleaning and filtering
  dataCleanfile <- paste(folder, "outputData/dataCleaningDetails.txt",sep="")
  sink(dataCleanfile, append = TRUE)
  cat(as.character(Sys.time()))
  cat("\n\nFunction 1.2: remove_unspecified \n\n")
  sink()


  # Set file name of raw data in input folder
  AA=unlist(raw_files)

  fwd_input <- paste(folder, AA[grep("_R1_001.fastq.gz",AA,ignore.case = TRUE)],sep="")
  rev_input <- paste(folder, AA[grep("_R2_001.fastq.gz",AA,ignore.case = TRUE)],sep="")


  # Set file name of raw data in intermediate folder
  path_no_unspec <- paste(folder, "outputData/intermediate/no_Ns",sep="")


  if (!dir.exists(path_no_unspec)) {dir.create(path_no_unspec)}
  fwd_no_unspec <- paste(path_no_unspec, gsub("^.*/", "", fwd_input), sep = "/")
  rev_no_unspec <- paste(path_no_unspec, gsub("^.*/", "", rev_input), sep = "/")


  # Skip step and if files already present
  if ((sum(!file.exists(fwd_no_unspec)) +
       sum(!file.exists(rev_no_unspec))) != 0) {

  # Remove Ns from the sequences and save as .gz files in intermediate folder

  dada2::filterAndTrim(fwd_input, fwd_no_unspec, rev_input, rev_no_unspec,
                       maxN = 0, multithread = FALSE, verbose = FALSE)


  # DATA FLAG: flag if <15000 reads per sample before denoising

  names <- c()

  for (x in seq_along(fwd_input)) {
    fw_in_L <- length(ShortRead::id(ShortRead::readFastq(fwd_input[x])))
    if (fw_in_L < 15000) {
      names <- append(names, sub(".*//", "", sub("_R1.*", "", fwd_input[x])))
    }
  }

  if (length(names) > 0) {
    messageColour("Low raw read count for following samples.\n", "warning")
messageColour("Confidence in prediction may be low. Consider investigation of
raw read count and denoising statistics. \n", "warning")
for (y in names){
  messageColour(y, "warningSamples")
  messageColour("\n", "warningSamples")
}
cat("\n")

  sink(dataCleanfile, append = TRUE)
  cat(as.character(Sys.time()))
  cat("\nLow raw read count for following samples.\n")
  cat("Confidence in prediction may be low. Consider investigation of raw read")
  cat(" count and denoising statistics.\n")
  for (y in names){
    cat(y, "\n")
  }
  cat("\n\n")
  sink()
  }

  # Add details of Ns cleaning to file
  fw_in_L <- c()#create empty vector
  fw_out_L <- c()
  sink(dataCleanfile, append = TRUE)
  cat(as.character(Sys.time()))
  for (x in seq_along(fwd_input)) {
    fw_in_L[x] <- length(ShortRead::id(ShortRead::readFastq(fwd_input[x])))
    fw_out_L[x] <- length(ShortRead::id(ShortRead::readFastq(fwd_no_unspec[x])))

    cat("\nFiles ", sub(folder, "", sub(".fastq.gz", "", fwd_input[x])), " & ",
        sub(folder, "", sub(".fastq.gz", "", rev_input[x])), " - read in ",
        fw_in_L[x], " paired-sequences, output ", fw_out_L[x], " (",
        round(100 * (fw_out_L[x] / fw_in_L[x]), 1),
        "%) filtered paired-sequences ", sep = "")
  }
  cat("\n\n")
  sink()

  ####DATA FLAG: flag if data loss >70% ####
  if (length(names[(fw_out_L / fw_in_L) < 0.3]) > 0) {
messageColour("High loss of reads when removing unspecified or ambigous bases in
following samples, data of low quality. Treat IQI predictions with caution.",
              "warning")
    for (y in names[(fw_out_L / fw_in_L) < 0.3]) {
      messageColour(y, "warningSamples")
      messageColour("\n", "warningSamples")
    }
    cat("\n")

    sink(dataCleanfile, append = TRUE)
    cat(as.character(Sys.time()))
    cat("\nHigh loss of reads when removing unspecified or ambigous bases in ")
    cat("following samples, data of low quality. \n")
    cat("Data loss >70% of reads. Treat IQI predictions with caution.\n")
    for (y in names[(fw_out_L / fw_in_L) < 0.3]) {
      cat(y, "\n")
    }
    cat("\n\n")
    sink()
  }

  ####Check all files created ####
  fwd_actual <- stringr::str_subset(dir(path_no_unspec), "_R1_001.fastq.gz")
  rev_actual <- stringr::str_subset(dir(path_no_unspec), "_R2_001.fastq.gz")
  #use grep and then ignore.case.##
  fwd_expected <- stringr::str_subset(raw_files, "_R1_001.fastq.gz")
  rev_expected <- stringr::str_subset(raw_files, "_R2_001.fastq.gz")

  fwd_actual <- intersect(fwd_actual, fwd_expected)
  rev_actual <- intersect(rev_actual, rev_expected)

      if (identical(fwd_actual, fwd_expected) &&
          identical(rev_actual, rev_expected)) {
      cat("All samples in batch processed.\n")

      sink(dataCleanfile, append = TRUE)
      cat(as.character(Sys.time()))
      cat("\nAll samples in batch processed.\n\n")
      sink()
      } else {
      cat("The following samples were not processed:\n")
      cat(fwd_expected[!(fwd_expected %in% fwd_actual)], "\n")
      cat(rev_expected[!(rev_expected %in% rev_actual)], "\n\n")

      sink(dataCleanfile, append = TRUE)
      cat(as.character(Sys.time()))
      cat("\nThe following samples were not processed:\n")
      cat(fwd_expected[!(fwd_expected %in% fwd_actual)], "\n")
      cat(rev_expected[!(rev_expected %in% rev_actual)], "\n\n")
      sink()
    }

if (!identical(sub("_R1_.*", "", fwd_actual), sub("_R2_.*", "", rev_actual))) {
    fwd_actual <- fwd_actual[(sub("_R1_.*", "", fwd_actual)
                              %in% sub("_R2_.*", "", rev_actual))]
    rev_actual <- rev_actual[(sub("_R2_.*", "", rev_actual)
                              %in% sub("_R1_.*", "", fwd_actual))]
  }

  fwd_no_unspec <- file.path(path_no_unspec, fwd_actual)
  rev_no_unspec <- file.path(path_no_unspec, rev_actual)
  }

  messageColour("Function finished: 'remove_unspecified' \n\n", "message")

  return(list(fwd_no_unspec, rev_no_unspec))
}

