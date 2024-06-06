#' Extract the initial line from fastq.gz files and report data frame with
#' unique instrument_run pairs and sample names.
#'
#' @param folder A string. Location where the raw fastq.gz files can be found.
#'
#' @importFrom magrittr %>%
#'
#' @return A data frame.
Step1.1_extract_identifier <- function(folder) {
  
  messageColour("Function starting: 'Step1.1_extract_identifier'
   Extracting unique instrument_run pairs and sample names from fastq files \n\n",
                "message")

  columns <- c("Instrument", "RunID", "FlowCellID", "Lane", "Tile", "X", "Y",
               "ReadNum", "FilterFlag", "Indexes", "SampleNumber")

  raw_files <- sort(fs::dir_ls(folder, glob = "*.gz"))

  run_metadata <- vroom::vroom(raw_files, col_names = columns, delim = ":",
                               n_max = 1, show_col_types = FALSE)

  raw_files_names <- sort(list.files(folder, pattern = ".gz",
                                     full.names = FALSE))

  run_metadata <- cbind(run_metadata, raw_files_names)


  # Make a dataframe with the unique Instrument_Run ID with the sample names
  raw_files <- data.frame(Unique = paste(run_metadata$Instrument,
                                         run_metadata$RunID, sep = "_"),
                         Filename = run_metadata$raw_files_names)

  #Determine what rows belong to one batch
  samples_row <- raw_files %>%
    dplyr::group_by(Unique) %>%
    dplyr::group_rows()

  # Make a list where each element a list of filenames per batch
  samples <- list()
  for (i in seq_along(samples_row)) {
    samples <- append(samples, list(raw_files$Filename[samples_row[[i]]]))
  }

  messageColour("Function finished: 'Step1.1_extract_identifier' \n\n", "message")

  return(samples)
}
