#' step0.5_capitaliseR1_R2
#'
#' Reads in raw sequence data, capitalises any "r1" and "r2" in the file metadata components of
#' an Illumina Miseq output filename.
#' This addresses software and filesharing bugs that lose case sensitivity, resulting in the
#' conversion of filenames to lower case.
#' The read orientation, indicated by R1 or R2 is essential for the bioinformatic process, and making it the
#' correct case saves further unnessecarily flexible code.
#'
#' @param folder, string naming the location of the raw FASTQ files to be processed
#'
#' @return NULL, code changes the capitalisation on the input files
#'
#' @section Example usage: step0.5_capitaliseR1_R2(folder)
step0.5_capitaliseR1_R2 <- function(folder) {

  # Set the path to the directory containing the files
  directory_path=folder

  # List all files in the directory
  files <- list.files(directory_path, full.names = TRUE)
  # fwd 
  # Iterate over each file and rename if needed
  for (file in files) {
    # Extract the file name without path
    file_name <- basename(file)
      # Check if the file name contains the pattern "_R1_00"
    if (grepl("_r1_00", file_name)) {
      # Replace lowercase 'r' with uppercase 'R'
      new_file_name <- gsub("_r1_00", "_R1_00", file_name)
      # Construct the new file path
      new_file_path <- file.path(dirname(file), new_file_name)
      # Rename the file
      file.rename(file, new_file_path)
      cat("Renamed:", file, "to", new_file_path, "\n")
    }
  }
  # rev
  # Iterate over each file and rename if needed
  for (file in files) {
    # Extract the file name without path
    file_name <- basename(file)
      # Check if the file name contains the pattern "_R1_00"
    if (grepl("_r2_00", file_name)) {
      # Replace lowercase 'r' with uppercase 'R'
      new_file_name <- gsub("_r2_00", "_R2_00", file_name)
      # Construct the new file path
      new_file_path <- file.path(dirname(file), new_file_name)
      # Rename the file
      file.rename(file, new_file_path)
      cat("Renamed:", file, "to", new_file_path, "\n")
    }
  }

}
