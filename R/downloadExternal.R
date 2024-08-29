#' Download the cutadapt.exe file and SILVA 138.1 taxa
#' reference database from the SAMS THREDDS server, or
#' from source.
#'
#'
#' @export
#' @return NULL
#'
#' @section Example usage: eDNA2IQI::downloadExternal()
#'

downloadExternal <- function() {


####taxa database

  if (!file.exists(file.path(find.package("eDNA2IQI"), "extdata", "referenceDatabase_full.fa.gz"))) {

    # Prompt the user for permission to download the file
    user_input <- readline(prompt = "The reference database file that needs to be stored locally is missing. Do you want to attempt to download it? (y/n): ")

    if (tolower(user_input) == "y") {

      tryCatch({
        # Attempt to download the file from SAMS Thredds server
        url <- "https://thredds.sams.ac.uk/thredds/fileServer/Full_database/referenceDatabase_full.fa.gz"
        filePath <- file.path(find.package("eDNA2IQI"), "extdata", "referenceDatabase_full.fa.gz")
        options(timeout = 300)
        messageColour("Updating taxa reference database from SAMS Thredds server \n\n", "warnMessage")
        utils::download.file(url, filePath, quiet = TRUE)
        message("SAMS Thredds server download successful")
      }, error = function(e) {
        # Handle error if SAMS Thredds download fails
        message("\nSAMS Thredds server download unsuccessful \n\n", conditionMessage(e))
        message("\n\n Attempting download from Zenodo, may be unstable in the future")

        tryCatch({
          # Attempt to download the file from Zenodo server
          url <- "https://zenodo.org/records/4587955/files/silva_nr99_v138.1_train_set.fa.gz?download=1"
          filePath <- file.path(find.package("eDNA2IQI"), "extdata", "referenceDatabase_full.fa.gz")
          options(timeout = 300)
          messageColour("Updating taxa reference database from Zenodo server \n\n", "warnMessage")
          utils::download.file(url, filePath, quiet = TRUE)
          message("Zenodo server download successful")
        }, error = function(e) {
          # Handle error if Zenodo download fails
          message("\nZenodo server download unsuccessful \n\n", conditionMessage(e))
          messageColour("\n
                    Automatic downloads have failed \n
                    Try manual download from: \n
                    https://thredds.sams.ac.uk/thredds/fileServer/Full_database/referenceDatabase_full.fa.gz, \n
                    or \n
                    https://zenodo.org/records/4587955/files/silva_nr99_v138.1_train_set.fa.gz \n\n
                    and save file as `referenceDatabase_full.fa.gz` \n
                    in the /extdata/ directory in the eDNA2IQI R library files", "message")
          message("\n")
        })
      })
    } else {
      stop(messageColour("Reference database needs to be downloaded to continue with taxa allocation. \n
          Proceed with automatic download, or manual download from: \n
                    https://thredds.sams.ac.uk/thredds/fileServer/Full_database/referenceDatabase_full.fa.gz, \n
                    or \n
                    https://zenodo.org/records/4587955/files/silva_nr99_v138.1_train_set.fa.gz \n\n
                    and save file as `referenceDatabase_full.fa.gz` \n
                    in the /extdata/ directory in the eDNA2IQI R library files \n", "message"))
    }


    if (!file.exists(file.path(find.package("eDNA2IQI"), "extdata", "referenceDatabase_full.fa.gz"))) {
      # If the file does not exist, stop execution with an error message
      message("\n Reference taxa database not successfully downloaded \n")
    }

  }

else{
    message("\nReference taxa database already downloaded")
}

##cutadapt

  if (!file.exists(file.path(find.package("eDNA2IQI"), "extdata", "cutadapt_v1.exe"))) {

    # Prompt the user for permission to download the file
    user_input <- readline(prompt = "The Cutadapt .exe file that needs to be stored locally is missing. Do you want to attempt to download it? (y/n): ")

    if (tolower(user_input) == "y") {

      tryCatch({
        # Attempt to download the file from SAMS Thredds server
        url <- "https://thredds.sams.ac.uk/thredds/fileServer/cutadapt_executable/cutadapt_v1.exe"
        filePath <- file.path(find.package("eDNA2IQI"), "extdata", "cutadapt_v1.exe")
        options(timeout = 300)
        messageColour("Updating Cutadapt .exe file from SAMS Thredds server \n\n", "warnMessage")
        utils::download.file(url, filePath, quiet = TRUE)
        message("SAMS Thredds server download successful")
      }, error = function(e) {
        # Handle error if SAMS Thredds download fails
        message("\nSAMS Thredds server download unsuccessful \n\n", conditionMessage(e))
        message("\n\n Attempting download from Github, may be unstable in the future")

        tryCatch({
          # Attempt to download the file from Github
          url <- "https://github.com/marcelm/cutadapt/releases/download/v4.0/cutadapt.exe"
          filePath <- file.path(find.package("eDNA2IQI"), "extdata", "referenceDatabase_full.fa.gz")
          options(timeout = 300)
          messageColour("Updating Cutadapt .exe file from Github \n\n", "warnMessage")
          utils::download.file(url, filePath, quiet = TRUE)
          message("Github download successful")
        }, error = function(e) {
          # Handle error if Github download fails
          message("\nGithub server download unsuccessful \n\n", conditionMessage(e))
          messageColour("\n
                    Automatic downloads have failed \n
                    Try manual download from: \n
                    https://thredds.sams.ac.uk/thredds/fileServer/cutadapt_executable/cutadapt_v1.exe, \n
                    or \n
                    https://github.com/marcelm/cutadapt/releases/download/v4.0/cutadapt.exe \n\n
                    and save file as `cutadapt_v1.exe` \n
                    in the /extdata/ directory in the eDNA2IQI R library files", "message")
          message("\n")
        })
      })
    } else {
      stop(messageColour("Reference database needs to be downloaded to continue with taxa allocation. \n
          Proceed with automatic download, or manual download from: \n
                    https://thredds.sams.ac.uk/thredds/fileServer/cutadapt_executable/cutadapt_v1.exe, \n
                    or \n
                    https://github.com/marcelm/cutadapt/releases/download/v4.0/cutadapt.exe \n\n
                    and save file as `cutadapt_v1.exe` \n
                    in the /extdata/ directory in the eDNA2IQI R library files \n", "message"))
    }


    if (!file.exists(file.path(find.package("eDNA2IQI"), "extdata", "cutadapt_v1.exe"))) {
      # If the file does not exist, stop execution with an error message
      message("\n Cutadapt .exe file not successfully downloaded \n")
    }

  }

  else{
    message("\nCutadapt .exe file already downloaded")
  }


##insert downloading of trained models here if too big for github download
}
