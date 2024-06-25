#' Generates the predicted IQI values based on the optimized random forest
#' model.
#'
#' @param S16_reads A data frame. Contains the bacterial taxa read counts which
#' should be used to predict the IQI values.
#'
#' @importFrom utils glob2rx
#'
#' @return A data frame. As inputted, with added 'predicted_IQI' column
#' @export
#'


step3.0_predict_iqi <- function(S16_reads) {
  messageColour("Function starting: 'predict_iqi'
  Predicting IQI using optimized random forest model  \n\n", "message")
  

  # Find trained random forest, extract info on taxa level and rarefaction rate
  potenRfs <- dir(file.path(find.package("eDNA2IQI"), "extdata"),
                  full.names = TRUE)
  RF <- subset(potenRfs, grepl(glob2rx("*/RF*.Rdata"), potenRfs,ignore.case = TRUE))
  load(RF);class(RF)
  RF=tolower(RF)
  (rarefaction_rate <- as.numeric(stringr::str_extract_all(RF,
                                            "(?<=rf_rarefy_).+(?=_taxalevel)")))


  # Check for samples  with fewer reads than rarefaction rate
  too_low <- which(rowSums(S16_reads) < rarefaction_rate)


    if (!rlang::is_empty(too_low)) {
    messageColour("IQI can not be predicted for the following samples:\n",
                  "warning")
    messageColour(names(too_low),"warningSamples")
    messageColour("\nRead count under ", "warning")
    messageColour(rarefaction_rate, "warning")
    S16_reads <- S16_reads[!rownames(S16_reads)%in%names(too_low),]
    
    }

   

   if(nrow(S16_reads)==0){message("\nNo samples exceeded rarefaction limit, breaking");return(NULL)}
  # Rarefy data
  #output is a matrix
    rare_data <- vegan::rrarefy(S16_reads, rarefaction_rate)
    rare_data <- as.data.frame(rare_data)

  # Clean sample names
  rownames(rare_data)=sub("_R1.*", "", rownames(rare_data),ignore.case = TRUE)

  
  if(ncol(rare_data)<2){message("\nOnly one or fewer taxa match RF taxa, breaking");return(NULL)}

  rare_data <- rare_data[,c( colnames(rare_data) %in% RF_final_reduced$coefnames)]

  # Check if taxa used for RF predictions are present in any sample, if not add
  # them as zero values

  if (FALSE %in% (RF_final_reduced$coefnames %in% colnames(rare_data))) {
    name <- RF_final_reduced$coefnames[!(RF_final_reduced$coefnames
                                         %in% colnames(rare_data))]
  

    messageColour("The following taxa were used to train the randomForest but
are not present in the testing data: \n\n", "warning")
    for (x in name) {
      messageColour(x, "warning")
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
        messageColour(x, "warning")
        cat("\n")
      }
      messageColour("Treat IQI predictions with caution. \n\n", "warning")
    }
  }

  # Make predictions with the trained random forest
  rare_data$predicted_IQI <- stats::predict(RF_final_reduced, rare_data)
  
  rare_data=rare_data[, c("predicted_IQI", setdiff(names(rare_data), "predicted_IQI"))]
  
  return(rare_data)
}
