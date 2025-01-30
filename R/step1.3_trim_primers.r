#' Step1.3_trim_primers
#' Remove the primer sequences from the sequence reads
#'
#' @param folder A string. Location of raw data files.
#' @param fwd_no_unspec A character vector. Contains the file location of
#' forward reads of the input samples where Ns were removed.
#' @param rev_no_unspec A character vector. Contains the file location of
#' reverse reads of the input samples where Ns were removed.
#' @param parameteroptions A dataframe.
#' @param auto_download Boolean to override console prompt to automatically
#'   download external dependencies. This can be useful for automating or
#'   testing this function.
#'
#' @return A list.
#'
step1.3_trim_primers <- function(folder, fwd_no_unspec, rev_no_unspec,
                         parameteroptions, auto_download = FALSE) {

  messageColour("Function starting: 'trim_primers'
  Removing primer sequences from sequence reads \n\n", "message")

  #Set file to store details of cleaning and filtering
  dataCleanfile <- paste(folder, "outputData/dataCleaningDetails.txt",sep="")
  sink(dataCleanfile, append = TRUE)
  cat(as.character(Sys.time()))
  cat("\n\nFunction 1.3: trim_primers \n\n")
  sink()

  #Change folder of no_Ns to trimmed
  path_trimmed <- paste(folder, "outputData/intermediate/", "trimmed",sep="")

  if (!dir.exists(path_trimmed)) {
    dir.create(path_trimmed)
  }
  
  #create filename for trimmed forward
  fwd_trimmed <- sub("no_Ns", "trimmed", fwd_no_unspec)
  #and reverse
  rev_trimmed <- sub("no_Ns", "trimmed", rev_no_unspec)

  # Skip step and if files already present
  if ((sum(!file.exists(fwd_trimmed)) + sum(!file.exists(rev_trimmed))) != 0) {

  #Load the cutadapt tool if required
  downloadExternal(auto_download = auto_download)
    if (!file.exists(file.path(find.package("eDNA2IQI"), "extdata", "cutadapt_v1.exe"))) {
    # If the file does not exist, stop execution with an error message
    stop()
  }

  cutadapt_exe <- file.path(find.package("eDNA2IQI"), "extdata/cutadapt_v1.exe")


  #Specify primer sequences
  p=extractParameters(parameteroptions,"global","primers")

  primers <- c(sub(",.*", "", p), sub(".*,", "", p))
  fwd_primer <- primers[1]
  rev_primer <- primers[2]

  #Specify reverse complement primer sequences, call dada2
  fwd_primer_rc <- dada2::rc(fwd_primer)
  rev_primer_rc <- dada2::rc(rev_primer)

  #Add flags to primers for cutadapt
  fwd_flags <- paste("-g", fwd_primer, "-a", rev_primer_rc)
  rev_flags <- paste("-G", rev_primer, "-A", fwd_primer_rc)

  #fwd_trimmed is the complete path
  #fwd_trimmed_wo_folder is the path minus the folder
    #output trimmed file name and location
  (fwd_trimmed_wo_folder <- sub(folder, "", fwd_trimmed))
  (rev_trimmed_wo_folder <- sub(folder, "", rev_trimmed))
    #input files for cutadapt, following trimming.
  (fwd_no_unspec_wo_folder <- sub(folder, "", fwd_no_unspec))
  (rev_no_unspec_wo_folder <- sub(folder, "", rev_no_unspec))

  #checks if 1st character a "/", if TRUE, removes it
  if (substr(fwd_trimmed_wo_folder[1], 1, 1) == "/") {
    fwd_trimmed_wo_folder <- sub("/", "", fwd_trimmed_wo_folder)
    rev_trimmed_wo_folder <- sub("/", "", rev_trimmed_wo_folder)
    fwd_no_unspec_wo_folder <- sub("/", "", fwd_no_unspec_wo_folder)
    rev_no_unspec_wo_folder <- sub("/", "", rev_no_unspec_wo_folder)
  }

  setwd(folder)
  #Set variables
    n=extractParameters(parameteroptions,"cutadapt","-n","Char2Vect")
    m=extractParameters(parameteroptions,"cutadapt","-m","Char2Vect")
    j=extractParameters(parameteroptions,"cutadapt","-j","Char2Vect")

  #Remove primer sequences off the sequences using cutadapt.exe
  #method changes with operating system 

    #detect operating system
op_sys <- Sys.info()[1]

  for (i in seq_along(fwd_no_unspec)) {
    #linux and mac
    if (op_sys =="Linux" || op_sys =="Darwin") {
    system2("wine", args = c(cutadapt_exe, fwd_flags, rev_flags, "--discard-untrimmed",
                                   "-n", n, "-m", m, "-j", j, "-o",
                      fwd_trimmed_wo_folder[i], "-p", rev_trimmed_wo_folder[i],
                      fwd_no_unspec_wo_folder[i], rev_no_unspec_wo_folder[i]
                                   ))
    }
    #windows
    else if (op_sys =="Windows") {
    system2(cutadapt_exe, args = c(fwd_flags, rev_flags, "--discard-untrimmed",
                                            "-n", n, "-m", m, "-j", j, "-o",
                                            fwd_trimmed_wo_folder[i], "-p", rev_trimmed_wo_folder[i],
                                            fwd_no_unspec_wo_folder[i], rev_no_unspec_wo_folder[i],
                                            "--quiet"))

               }
}

  # Calculations for data flags
  names <- c()
  for (x in seq_along(fwd_no_unspec)) {
    fw_in_L <- length(ShortRead::id(ShortRead::readFastq(fwd_no_unspec[x])))
    fw_out_L <- length(ShortRead::id(ShortRead::readFastq(fwd_trimmed[x])))
    if (fw_out_L / fw_in_L < 0.05) {
     names <- append(names, sub(folder, "", sub("_R1.*", "", fwd_no_unspec[x])))
    }
  }

  if (length(names) > 0) {
  messageColour("Primer sequence not present in following samples.\n",
                "warning")
  messageColour("Has it been pre-filtered, and the right primers been used?\n",
                "warning")
    for (y in names){
      messageColour(y, "warningSamples")
      messageColour("\n", "warningSamples")
    }
    cat("\n")

    sink(dataCleanfile, append = TRUE)
    cat(as.character(Sys.time()))
    cat("\nPrimer sequence not present in following samples.\n")
    cat("Has it been pre-filtered, and the right primers been used? \n")
    for (y in names){
      cat(y, "\n")
    }
    cat("\n\n")
    sink()
  }

  # Write details of cleaning
  fw_in_L <- c()
  fw_out_L <- c()
  sink(dataCleanfile, append = TRUE)
  cat(as.character(Sys.time()))
  for (x in seq_along(fwd_no_unspec)) {
    fw_in_L <- length(ShortRead::id(ShortRead::readFastq(fwd_no_unspec[x])))
    fw_out_L <- length(ShortRead::id(ShortRead::readFastq(fwd_trimmed[x])))
  cat("\nFiles ", sub(".*/no_Ns/", "", sub(".fastq.gz", "", fwd_no_unspec[x])),
        " & ", sub(".*/no_Ns/", "", sub(".fastq.gz", "", rev_no_unspec[x])),
        " - ", fw_in_L - fw_out_L, " reads removed from ", fw_in_L, "(",
        round(100 * (fw_in_L - fw_out_L) / fw_in_L, 1),
        "%), as no primers present", sep = "")

  }
  cat("\n\n")
  sink()


    # Check all files created
    fwd_actual <- stringr::str_subset(dir(path_trimmed), "_R1_001.fastq.gz")
    rev_actual <- stringr::str_subset(dir(path_trimmed), "_R2_001.fastq.gz")
    fwd_expected <- sub(".*/trimmed/", "", fwd_trimmed_wo_folder)
    rev_expected <- sub(".*/trimmed/", "", rev_trimmed_wo_folder)

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

    fwd_trimmed <- file.path(path_trimmed, fwd_actual)
    rev_trimmed <- file.path(path_trimmed, rev_actual)

  }

  messageColour("Function finished: 'trim_primers' \n\n", "message")

  return(list(fwd_trimmed, rev_trimmed))
}
