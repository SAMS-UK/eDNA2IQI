#' step0.5_capitaliseR1_R2
#'
#' Reads in the input data, capitalises any "r1" and "r2" in the file metadata componenets of 
#' an Illumina Miseq output filename
#' This addresses software and filesharing bugs that lose case sensitivity, resulting in the 
#' conversion of filenames to lower case. 
#' The read orientation, indicated by R1 or R2 is essential for the bioinformatic process, and making it the 
#' correct case saves further unnessecarily flexible code.
#'
#' @param folder, string naming the location of the raw FASTQ files to be processed
#'
#' @return NULL, code changes the capitalisation on the input files
#' @export
#'
#' @examples step0.5_capitaliseR1_R2(folder)
step0.5_capitaliseR1_R2 <- function(folder) {

  # Set the path to the directory containing the files
  directory_path=folder

  # List all files in the directory
  files <- list.files(directory_path, full.names = TRUE)

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

#' B_fastparameteroptions
#' Function author: TomW
#' 
#' updates parameteroptions to 'fast' options, for debugging
#'
#' @param parameterOptions
#' parameterOptions is loaded as part of the package
#'
#' @return dataframe, with fast parameter options for code testing
#' 
#'
#' @examples B_fastparameteroptions(parameterOptions)
B_fastparameteroptions <- function(parameterOptions){
#create debugging parameters
#turn off quality check if it is working
parameterOptions[3,3] = "FALSE"

# change to current optimal truncation
parameterOptions[24,3] = "c(272,186)"

#speed up error learning
parameterOptions[27,3] = "1.00E+03"

#single thread DADA2
#optional
#parameterOptions[49,3] = "1"

#make everything verbose
parameterOptions[26,3] = "TRUE"
parameterOptions[34,3] = "TRUE"
parameterOptions[35,3] = "TRUE"
parameterOptions[43,3] = "TRUE"
parameterOptions[48,3] = "TRUE"
parameterOptions[50,3] = "TRUE"
return(parameterOptions)
}


