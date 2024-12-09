#' Download the cutadapt.exe file and SILVA 138.1 taxa
#' reference database from the SAMS THREDDS server, or
#' from source.
#' @param auto_download Boolean to override console prompt to automatically
#'   download external dependencies. This can be useful for automating or
#'   testing this function.
#'
#' @export
#' @return NULL
#'
#' @section Example usage: eDNA2IQI::downloadExternal()
#'

downloadExternal <- function(auto_download = FALSE) {


  # File paths
  filePathDB <- file.path(find.package("eDNA2IQI"), "extdata", "referenceDatabase_full.fa.gz")
  filePathCutadapt <- file.path(find.package("eDNA2IQI"), "extdata", "cutadapt_v1.exe")
  filePathModel <- file.path(find.package("eDNA2IQI"), "extdata", "RFModel_rarefy_5000_taxalevel_family_BMB_Run1_4b.rda")
  
  # Download reference database
  if (!file.exists(filePathDB)) {
    if(auto_download == FALSE) {
    # Prompt user for permission to download
    user_input <- readline(prompt = "The reference database file that needs to be stored locally is missing. Do you want to attempt to download it? (y/n): ")
    } else {
      user_input <- "y"
    }

    if (tolower(user_input) == "y") {

      tryCatch({
        # Attempt to download from SAMS Thredds server
        url <- "https://thredds.sams.ac.uk/thredds/fileServer/Full_database/referenceDatabase_full.fa.gz"
        options(timeout = 300)
        messageColour("Updating taxa reference database from SAMS Thredds server \n\n", "warnMessage")
        utils::download.file(url, filePathDB, quiet = TRUE)

        # Check if file was successfully downloaded
        if (file.exists(filePathDB)) {
          message("Download of the reference taxa database from the SAMS Thredds server successful")
        }

      }, error = function(e) {
        # Handle error if SAMS Thredds download fails
        message("\nSAMS Thredds server download unsuccessful \n\n", conditionMessage(e))
        message("\n\n Attempting download from Zenodo, may be unstable in the future")

        # Nested tryCatch for downloading from Zenodo if the first attempt fails
        tryCatch({
          url <- "https://zenodo.org/records/4587955/files/silva_nr99_v138.1_train_set.fa.gz"
          options(timeout = 300)
          messageColour("Updating taxa reference database from Zenodo server \n\n", "warnMessage")
          utils::download.file(url, filePathDB, quiet = TRUE)

          # Check if file was successfully downloaded
          if (file.exists(filePathDB)) {
            message("Zenodo server download successful")
          }

        }, error = function(e) {
          # Handle error if Zenodo download fails
          message("\nZenodo server download unsuccessful \n\n", conditionMessage(e))
          messageColour(paste(
            "\nAutomatic downloads have failed \n",
            "Try manual download from: \n",
            "https://thredds.sams.ac.uk/thredds/fileServer/Full_database/referenceDatabase_full.fa.gz, \n",
            "or \n",
            "https://zenodo.org/records/4587955/files/silva_nr99_v138.1_train_set.fa.gz \n\n",
            "and save file as `referenceDatabase_full.fa.gz` \n",
            "in the /extdata/ directory in the eDNA2IQI R library files"
          ), "message")
          message("\n")
        })
      })

    } else {
      message(messageColour(paste(
        "Reference database needs to be downloaded to continue with taxa allocation. \n",
        "Proceed with automatic download, or manual download from: \n",
        "https://thredds.sams.ac.uk/thredds/fileServer/Full_database/referenceDatabase_full.fa.gz, \n",
        "or \n",
        "https://zenodo.org/records/4587955/files/silva_nr99_v138.1_train_set.fa.gz \n\n",
        "and save file as `referenceDatabase_full.fa.gz` \n",
        "in the /extdata/ directory in the eDNA2IQI R library files \n"
      ), "message"))
    }

    if (!file.exists(filePathDB)) {
      # If the file does not exist, notify the user
      message("\nReference taxa database not successfully downloaded \n")
    }

  } else {
    message("\nReference taxa database already downloaded")
  }

  # Download Cutadapt
  if (!file.exists(filePathCutadapt)) {
    if(auto_download == FALSE) {
    # Prompt user for permission to download
    user_input <- readline(prompt = "The Cutadapt .exe file that needs to be stored locally is missing. Do you want to attempt to download it? (y/n): ")
    } else {
      user_input <- "y"
}
    if (tolower(user_input) == "y") {

      tryCatch({
        # Attempt to download from SAMS Thredds server
        url <- "https://thredds.sams.ac.uk/thredds/fileServer/cutadapt_executable/cutadapt_v1.exe"
        options(timeout = 300)
        messageColour("Updating Cutadapt .exe file from SAMS Thredds server \n\n", "warnMessage")
        utils::download.file(url, filePathCutadapt, quiet = TRUE, mode = "wb")

        # Check if file was successfully downloaded
        if (file.exists(filePathCutadapt)) {
          message("Cutadapt .exe downlaod from SAMS Thredds server download successful")
        }

      }, error = function(e) {
        # Handle error if SAMS Thredds download fails
        message("\nSAMS Thredds server download unsuccessful \n\n", conditionMessage(e))
        message("\n\n Attempting download from Github, may be unstable in the future")

        # Nested tryCatch for downloading from Github if the first attempt fails
        tryCatch({
          url <- "https://github.com/marcelm/cutadapt/releases/download/v4.0/cutadapt.exe"
          options(timeout = 300)
          messageColour("Updating Cutadapt .exe file from Github \n\n", "warnMessage")
          utils::download.file(url, filePathCutadapt, quiet = TRUE, mode = "wb")

          # Check if file was successfully downloaded
          if (file.exists(filePathCutadapt)) {
            message("Cutadapt.exe download from Github successful")
          }

        }, error = function(e) {
          # Handle error if Github download fails
          message("\nCutadapt.exe download from Github unsuccessful \n\n", conditionMessage(e))
          messageColour(paste(
            "\nAutomatic downloads have failed \n",
            "Try manual download from: \n",
            "https://thredds.sams.ac.uk/thredds/fileServer/cutadapt_executable/cutadapt_v1.exe, \n",
            "or \n",
            "https://github.com/marcelm/cutadapt/releases/download/v4.0/cutadapt.exe \n\n",
            "and save file as `cutadapt_v1.exe` \n",
            "in the /extdata/ directory in the eDNA2IQI R library files"
          ), "message")
          message("\n")
        })
      })

    } else {
      message(messageColour(paste(
        "Cutadapt.exe needs to be downloaded to continue with data denoising. \n",
        "Proceed with automatic download, or manual download from: \n",
        "https://thredds.sams.ac.uk/thredds/fileServer/cutadapt_executable/cutadapt_v1.exe, \n",
        "or \n",
        "https://github.com/marcelm/cutadapt/releases/download/v4.0/cutadapt.exe \n\n",
        "and save file as `cutadapt_v1.exe` \n",
        "in the /extdata/ directory in the eDNA2IQI R library files \n"
      ), "message"))
    }

    if (!file.exists(filePathCutadapt)) {
      # If the file does not exist, notify the user
      message("\nCutadapt.exe file not successfully downloaded \n")
    }

  } else {
    message("\nCutadapt.exe file already downloaded")
  }


  # Download RF_model from Thredds
  
  
 if (!file.exists(filePathModel)) {
    if (auto_download == FALSE) {
      # Prompt user for permission to download
      user_input <- readline(prompt = "The Random Forest model that needs to be stored locally is missing. Do you want to attempt to download it? (y/n): ")
    } else {
      user_input <- "y"
    }

    if (tolower(user_input) == "y") {
      tryCatch({
        # Attempt to download from SAMS Thredds server
        url <- "https://thredds.sams.ac.uk/thredds/fileServer/RF_Model/RFModel_rarefy_5000_taxalevel_family_BMB_Run1_4b.rda"
        options(timeout = 300)
        messageColour("Updating Random Forest Model file from SAMS Thredds server \n\n", "warnMessage")
        utils::download.file(url, filePathModel, quiet = TRUE, mode = "wb")

        # Check if file was successfully downloaded
        if (file.exists(filePathModel)) {
          message("Random Forest Model download from SAMS Thredds server successful")
        }

      }, error = function(e) {
        # Handle error if SAMS Thredds download fails
        message("\nSAMS Thredds server download unsuccessful \n\n", conditionMessage(e))
        message("\n\n Try again soon, or contact SAMS or SEPA for assistance")
      })

    } else {
      message(messageColour(paste(
        "The Random Forest Models need to be downloaded to predict IQIs. \n",
        "Proceed with automatic download, or manual download from: \n",
        "https://thredds.sams.ac.uk/thredds/fileServer/RF_Model/RFModel_rarefy_5000_taxalevel_family_BMB_Run1_4b.rda, \n",
        "and save file as `RFModel_rarefy_5000_taxalevel_family_BMB_Run1_4b.rda` \n",
        "in the /extdata/ directory in the eDNA2IQI R library files \n"
       ), "message"))
    }
  }
  else {
    message("\nRF model file already downloaded")
  }
  
}
