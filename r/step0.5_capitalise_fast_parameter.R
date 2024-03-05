#' A_CapitaliseR1_R2
#' Reads in the input data, capitalises any r1 and r2 in string L001_R1 etc
#' convenient function, to address MS Explorer converting files to lower case during transfer.
#' @param fastqfolder string naming the location of the files to be processed
#'
#' @return NULL, code changes the capitalisation on the input files
#' @export
#'
#' @examples CapitaliseR1_R2(path_to_fastq_folder)
Step0.5_capitaliseR1_R2 <- function(fastqfolder) {

  # Set the path to the directory containing the files
  directory_path=fastqfolder

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
#' Last edit date: 14 Jan 24
#' updates parameteroptions to 'fast' options, for debugging
#'
#' @param parameterOptions
#' parameterOptions is loaded as part of the package
#'
#' @return dataframe, with fast parameter options for code testing
#' @export
#'
#' @examples fastparameteroptions(parameterOptions)
B_fastparameteroptions=function(parameterOptions){
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


