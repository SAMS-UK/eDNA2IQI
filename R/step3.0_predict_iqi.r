#' Generates the predicted IQI values based on the optimized random forest
#' model.
#'
#' @param folder A string. Location of raw fastq data files.
#' @param S16_reads A data frame. Contains the bacterial taxa read counts which
#' should be used to predict the IQI values.
#'
#' @importFrom utils glob2rx
#'
#' @return A data frame. As inputted, with added 'predicted_IQI' column
#' @export
#' @section Example usage: myPREDs <- step3.0_predict_iqi(folder = "C:/full/file/path/to/fastq/S16_reads", S16_reads = myTAXA)
#'

step3.0_predict_iqi <- function(folder, S16_reads) {
  messageColour("Function starting: 'predict_iqi'
  Predicting IQI using optimized random forest model  \n\n", "message")
  step3_start_time <- as.character(Sys.time())
  message("Step 3 Start Date/Time:",step3_start_time)

####Check Folder name and add / at the end if its not there
if (substr(folder, nchar(folder), nchar(folder)) != "/") {
  folder <- paste0(folder, "/")
}

  # Find trained random forest, extract info on taxa level and rarefaction rate
  potenRfs <- dir(file.path(find.package("eDNA2IQI"), "extdata"),
                  full.names = TRUE)
  RF <- subset(potenRfs, grepl(glob2rx("*/RF*.Rdata"), potenRfs,ignore.case = TRUE))
  load(RF);class(RF)
  RF=tolower(RF)
  rarefaction_rate <- as.numeric(stringr::str_extract(RF, "(?<=rarefy_)[0-9]+"))
  taxalevel <- as.character(stringr::str_extract(RF, "(?<=_taxalevel_)[^\\.]+"))

  # Check for samples  with fewer reads than rarefaction rate
  too_low <- which(rowSums(S16_reads) < rarefaction_rate)
  read_totals <- rowSums(S16_reads)
  readcounts_df <- data.frame("Denoised_Reads" = read_totals, row.names = rownames(S16_reads))

  # Store samples with low reads in a separate data frame
  excluded_samples <- data.frame(SampleID = names(too_low),
                                 Predicted_IQI = NA,
                                 stringsAsFactors = FALSE)
  row.names(excluded_samples) <- excluded_samples$SampleID
  excluded_samples$SampleID <- NULL

    if (!rlang::is_empty(too_low)) {
    messageColour("\nRead count under ", "warning")
    messageColour(rarefaction_rate, "warning")
    messageColour(" for the following samples. IQI was not calculated.\n", "warning")

    messageColour(paste(names(too_low), collapse = "\n"), "warningSamples")

    #remove samples under rarefaction
    S16_reads <- S16_reads[!rownames(S16_reads)%in%names(too_low),]

    }

   if(nrow(S16_reads)==0){
     message("\nNo samples exceeded rarefaction limit. Exiting function.")
     return(NULL)
   }

  # Rarefy data
  #output is a matrix
	set.seed(123)
  rare_data <- suppressWarnings(vegan::rrarefy(S16_reads, rarefaction_rate))
  rare_data <- as.data.frame(rare_data)
	for_barplot <- rare_data

  # write rarefied dataframe
	  utils::write.csv(rare_data, file = file.path(folder,
               "/outputData/rarefied_taxa_allocated_reads_", taxalevel, ".csv", fsep = ""),row.names = TRUE)



  # Clean sample names
  rownames(rare_data)=sub("_R1.*", "", rownames(rare_data),ignore.case = TRUE)


  if(ncol(rare_data)<2){message("\nOnly one or fewer taxa match RF taxa, breaking");return(NULL)}

  rare_data <- rare_data[,c( colnames(rare_data) %in% RF_final_reduced$coefnames)]

  # Check if taxa used for RF predictions are present in any sample, if not add
  # them as zero values

  if (FALSE %in% (RF_final_reduced$coefnames %in% colnames(rare_data))) {
    name <- RF_final_reduced$coefnames[!(RF_final_reduced$coefnames
                                         %in% colnames(rare_data))]


    messageColour("\n\nThe following taxa were used to train the randomForest but
are not present in the testing data: \n\n", "warning")
    for (x in name) {
      messageColour(x, "warningSamples")
      cat("\n")
    }
    messageColour("\nTreat IQI predictions with caution. \n\n", "warning")

    # Add columns with zero values for missing taxa
    rare_data[name] <- 0
  }

# Check if taxa used for RF predictions are present in individual samples

  name <- RF_final_reduced$coefnames[!(RF_final_reduced$coefnames
                                         %in% colnames(rare_data))]
  for (i in seq_along(rownames(rare_data))) {
    missingS <- colnames(rare_data)[(rare_data[i, ] == 0) &
                                      (!(colnames(rare_data) %in% name))]
    if (length(missingS) > 0) {
      messageColour("The following taxa were used to train the randomForest but are not present in sample ", "warning")
      messageColour(rownames(rare_data)[i], "warning")
      messageColour(": \n", "warning")
      for (x in missingS) {
        messageColour(x, "warningSamples")
        cat("\n")
      }
      messageColour("Treat IQI predictions with caution. \n\n", "warning")
    }


  # Make predictions with the trained random forest
  rare_data$Predicted_IQI <- stats::predict(RF_final_reduced, rare_data)

  rare_data=rare_data[, c("Predicted_IQI", setdiff(names(rare_data), "Predicted_IQI"))]


  #  Make sure excluded_samples has all columns in rare_data
  all_columns <- names(rare_data)

  # Add missing columns to excluded_samples with NA values
  for (col in all_columns) {
    if (!col %in% names(excluded_samples)) {
      excluded_samples[[col]] <- NA
    }
  }

  # Ensure column order matches between rare_data and excluded_samples
  excluded_samples <- excluded_samples[, all_columns]

  #clean names
  rownames(excluded_samples)=sub("_R1.*", "", rownames(excluded_samples),ignore.case = TRUE)
  rownames(readcounts_df)=sub("_R1.*", "", rownames(readcounts_df),ignore.case = TRUE)

   # Combine the data frames
  final_data <- rbind(rare_data,excluded_samples)

  final_data$Denoised_Reads <- readcounts_df$Denoised_Reads[match(rownames(final_data), rownames(readcounts_df))]
  final_data$eDNA2IQI_Version <- packageVersion("eDNA2IQI")
  final_data$Minimum_read_requirement <- rarefaction_rate
  #move reads and version to first col
  final_data <- final_data[, c("eDNA2IQI_Version", "Denoised_Reads", "Minimum_read_requirement", setdiff(names(final_data), c("eDNA2IQI_Version", "Denoised_Reads", "Minimum_read_requirement")))]


  #write file
  utils::write.csv(final_data, file = file.path(folder,
               "/outputData/predicted_IQIs.csv", fsep = ""),row.names = TRUE)

  # Make rarefied data taxaplot

S16_readsB=for_barplot
S16_readsB$SampleID=rownames(S16_readsB)
A=grep("SampleID",colnames(S16_readsB))

#Add a column for SampleID for merging with excluded samples
#rare_data$SampleID <- rownames(rare_data)

#move SampleID to first column
S16_readsB=S16_readsB[,c(A,1:(A-1))]

drawBarplot(S16_readsB, folder)
message("Step 3 Start Date/Time:",step3_start_time)
message("Step 3 End Date/Time:",as.character(Sys.time()))
return(rare_data)
 }
}

