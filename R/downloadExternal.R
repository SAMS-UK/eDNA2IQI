#' Download the random forest models,  cutadapt.exe file and SILVA 138.1 taxa
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


# Handle thredds 1st time user reg, and store token for subsequent downloads
# Download with token; prefers libcurl headers; falls back to curl package

.flask_base      <- "http://valetudo:8447" # this will be https://thredds.sams.ac.uk eventually
.register_url  <- paste0(.flask_base, "/api/register")  
.api  <- paste0(.flask_base, "/api/download") # GET with Bearer token is fine


# load token creds
flaskcli_load_cfg <- function(cfgpath) {
  p <- path.expand(cfgpath)
  if (file.exists(p)) fromJSON(p, simplifyVector = TRUE) else NULL
}

# save token creds
flaskcli_save_cfg <- function(cfg, cfgpath) {
  pdir <- dirname(path.expand(cfgpath))
  if (!dir.exists(pdir)) dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
  writeLines(jsonlite::toJSON(cfg, auto_unbox = TRUE, pretty = TRUE), path.expand(cfgpath))
}



# Get or create a token (prompts once)
# and registration if no token present!
flaskcli_token <- function(fileq, cfgpath) {
  cfg <- flaskcli_load_cfg(cfgpath)
  # Saved credentials is found so return token
  if (!is.null(cfg) && !is.null(cfg$token) && nzchar(cfg$token)) return(cfg$token)
  
  # No saved credentials so prompt for registration 

  # ---- helpers ----
.ask_menu <- function(title, choices, other_label = "Other") {
  full <- c(choices, other_label)
  repeat {
    cat("\n", title, "\n", sep = "")
    sel <- utils::menu(full, title = NULL, graphics = FALSE)
    if (sel == 0L) { cat("Please choose a number.\n"); next }
    if (sel <= length(choices)) return(list(value = choices[sel], is_other = FALSE))
    # "Other" chosen
    txt <- readline(paste0(other_label, " (type your answer): "))
    txt <- trimws(txt)
    if (nzchar(txt)) return(list(value = txt, is_other = TRUE))
    cat("Please enter a non-empty value.\n")
  }
}

# ------------------

cat("First-time registration:\n")
file <- fileq
name <- readline("YOUR NAME: ")



repeat {
  email <- readline("EMAIL: ")
  res <- sanitizeEmail(email)
  if (res$valid) {
    email <- res$value  # take the sanitized version
    break
  }
  cat("ehmmm — invalid email:", res$reason, "\n")
}


# Menus
sector_choices <- c("industry", "consultancy", "academic", "government/regulator")
env_choices    <- c("marine", "freshwater", "brackish", "terrestrial")
app_choices    <- c("aquaculture", "mining")

sector <- .ask_menu("Sector:", sector_choices)$value
application_env <- .ask_menu("Application environment:", env_choices)$value
application <- .ask_menu("Application:", app_choices)$value

# Free-text (still asked explicitly)
nationality <- readline("YOUR NATIONALITY: ")

cat("\nSubmitting registration to: ", .register_url, "\n", sep = "")

req <- httr2::request(.register_url) |>
  httr2::req_method("POST") |>
  httr2::req_body_json(list(
    name = name,
    email = email,
    sector = sector,
    application_env = application_env,
    application = application,
    nationality = nationality,
    file = file
  ))

resp <- httr2::req_perform(req)
httr2::resp_check_status(resp)
token <- httr2::resp_body_json(resp)$token
if (!nzchar(token)) stop("No token returned by server.", call. = FALSE)

flaskcli_save_cfg(list(
  token = token,
  name = name,
  email = email,
  sector = sector,
  application_env = application_env,
  application = application,
  nationality = nationality
), cfgpath)

}


# ensure users arent entering rubbish
sanitizeEmail <- function(email) {
  # handle NULL/NA/non-char
  if (is.null(email) || length(email) == 0L || is.na(email)) {
    return(list(valid = FALSE, value = NA_character_, reason = "missing"))
  }
  e <- as.character(email[[1]])
  
  # trim & strip wrappers
  e <- trimws(e)
  e <- sub("^<\\s*(.+?)\\s*>$", "\\1", e)                 # remove <> if present
  e <- gsub("[\u00A0\u200B\u200C\u200D\uFEFF]", "", e)     # strip NBSP/zero-width
  
  # must contain exactly one "@"
  parts <- strsplit(e, "@", fixed = TRUE)[[1]]
  if (length(parts) != 2) {
    return(list(valid = FALSE, value = e, reason = "must contain one @"))
  }
  
  local  <- parts[1]
  domain <- parts[2]
  
  # normalise dots & case
  local  <- gsub("\\.+", ".", local)
  domain <- tolower(gsub("\\.+", ".", domain))
  
  # no leading/trailing dots
  local  <- gsub("^\\.|\\.$", "", local)
  domain <- gsub("^\\.|\\.$", "", domain)
  
  # length limits (RFC-ish)
  if (nchar(local) == 0L || nchar(local) > 64L) {
    return(list(valid = FALSE, value = NA_character_, reason = "local-part length"))
  }
  if ((nchar(local) + 1 + nchar(domain)) > 254L) {
    return(list(valid = FALSE, value = NA_character_, reason = "address too long"))
  }
  
  # local-part allowed chars (dots already checked for edges/dupes)
  if (!grepl("^[A-Za-z0-9.!#$%&'*+/=?^_`{|}~-]+$", local)) {
    return(list(valid = FALSE, value = NA_character_, reason = "invalid local-part chars"))
  }
  
  # domain must have at least one dot and valid labels
  labels <- strsplit(domain, ".", fixed = TRUE)[[1]]
  if (length(labels) < 2L) {
    return(list(valid = FALSE, value = NA_character_, reason = "domain needs a dot"))
  }
  # label syntax & lengths
  if (any(nchar(labels) < 1L | nchar(labels) > 63L)) {
    return(list(valid = FALSE, value = NA_character_, reason = "domain label length"))
  }
  if (any(!grepl("^[A-Za-z0-9-]+$", labels))) {
    return(list(valid = FALSE, value = NA_character_, reason = "invalid domain chars"))
  }
  if (any(grepl("^-|-$", labels))) {
    return(list(valid = FALSE, value = NA_character_, reason = "domain label starts/ends with -"))
  }
  # TLD: letters only, min 2
  tld <- labels[length(labels)]
  if (!grepl("^[A-Za-z]{2,63}$", tld)) {
    return(list(valid = FALSE, value = NA_character_, reason = "invalid TLD"))
  }
  
  sanitized <- paste0(local, "@", paste(labels, collapse = "."))
  list(valid = TRUE, value = sanitized, reason = NA_character_)
}



flask_download <- function(file, destpath,cfgpath) {
  fileq <- utils::URLencode(file, reserved = TRUE)
  token <- flaskcli_token(fileq, cfgpath)
  url   <- paste0(.api, "?file=", fileq)
  print(url)
  
  utils::download.file(
    url      = url,
    destfile = destpath,
    method   = "libcurl",
    mode     = "wb",
    quiet    = FALSE,
    headers  = c(Authorization = paste("ednasams25", token))
  )
}





downloadExternal <- function(auto_download = FALSE) {


  # File paths
  filePathDB <- file.path(find.package("eDNA2IQI"), "extdata", "referenceDatabase_full.fa.gz")
  cfgpath <- file.path(find.package("eDNA2IQI"), "extdata", "thredds_creds.json") # needed to store credentials for thredds api after successful registration
  filePathCutadapt <- file.path(find.package("eDNA2IQI"), "extdata", "cutadapt_v1.exe")
  filePathModel <- file.path(find.package("eDNA2IQI"), "extdata", "RFModel_rarefy_5000_taxalevel_family_BMB_Run1_5_multiple.rda")

  # 1. Download reference database----
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
        threddsfile <- "/thredds/fileServer/Full_database/referenceDatabase_full.fa.gz"
        options(timeout = 300)
        messageColour("Updating taxa reference database from SAMS THREDDS server \n\n", "warnMessage")
        flask_download(threddsfile, filePathDB, cfgpath)

        # Check if file was successfully downloaded
        if (file.exists(filePathDB)) {
          message("Download of the reference taxa database from the SAMS THREDDS server successful")
        }

      }, error = function(e) {
        # Handle error if SAMS Thredds download fails
        message("\nSAMS THREDDS server download unsuccessful \n\n", conditionMessage(e))
        message("\n\n Attempting download from Zenodo repository. External links may change or become unavailable in the future, contact support to alert them if this occurs.")

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
    #message on denied download
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

  # 2. Download Cutadapt----
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
        messageColour("Updating Cutadapt .exe file from SAMS THREDDS server \n\n", "warnMessage")
        utils::download.file(url, filePathCutadapt, quiet = TRUE, mode = "wb")

        # Check if file was successfully downloaded
        if (file.exists(filePathCutadapt)) {
          message("Cutadapt .exe download from SAMS Thredds server download successful")
        }

      }, error = function(e) {
        # Handle error if SAMS Thredds download fails
        message("\nSAMS THREDDS server download unsuccessful \n\n", conditionMessage(e))
        message("\n\n Attempting download from the Cutadapt Github. External links may change or become unavailable in the future, contact support to alert them if this occurs.")

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
    #message on denied download
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


  #3. Download RF_model from SAMS THREDDs server----
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
        url <- "https://thredds.sams.ac.uk/thredds/fileServer/RF_Model/RFModel_rarefy_5000_taxalevel_family_BMB_Run1_5_multiple.rda"
        options(timeout = 300)
        messageColour("Updating Random Forest Model file from SAMS THREDDS server \n\n", "warnMessage")
        utils::download.file(url, filePathModel, quiet = TRUE, mode = "wb")

        # Check if file was successfully downloaded
        if (file.exists(filePathModel)) {
          message("Random Forest Model download from SAMS Thredds server successful")
        }

      }, error = function(e) {
        # Handle error if SAMS Thredds download fails
        message("\nSAMS THREDDS server download unsuccessful \n\n", conditionMessage(e))
        message("\n\n Try again soon, or contact SAMS or SEPA for assistance")
      })
    #message on denied download
    } else {
      message(messageColour(paste(
        "The Random Forest Models need to be downloaded to predict IQIs. \n",
        "Proceed with automatic download, or manual download from: \n",
        "https://thredds.sams.ac.uk/thredds/fileServer/RF_Model/RFModel_rarefy_5000_taxalevel_family_BMB_Run1_5_multiple.rda, \n",
        "and save file as `RFModel_rarefy_5000_taxalevel_family_BMB_Run1_5_multiple.rda` \n",
        "in the /extdata/ directory in the eDNA2IQI R library files \n"
       ), "message"))
    }
  }
  else {
    message("\nRF model file already downloaded")
  }

}
