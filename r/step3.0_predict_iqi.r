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
#makes for easier debugging to include the names of the objects returned by the relevant functions as
#arguments in the subsequent functions.  Easier to read and quicker to debug.
#don't use 'data'

Step3.0_predict_iqi <- function(S16_reads) {
  messageColour("Function starting: 'predict_iqi'
  Predicting IQI using optimized random forest model  \n\n", "message")
  #S16_reads=TaxAssigned1
  #glob2rx is a utils::function
  #it is used to convert wildcard expressions (globs) to regular expressions.

  # Find trained random forest, extract info on taxa level and rarefaction rate
  potenRfs <- dir(file.path(find.package("eDNA2IQI"), "extdata"),
                  full.names = TRUE)
  RF <- subset(potenRfs, grepl(glob2rx("*/RF*.Rdata"), potenRfs,ignore.case = TRUE))
  load(RF);class(RF)
  RF=tolower(RF)
  (rarefaction_rate <- as.numeric(stringr::str_extract_all(RF,
                                            "(?<=rf_rarefy_).+(?=_taxalevel)")))

  # Load trained random forest
  #S16_reads=S16_reads2

  # Check for samples  with fewer reads than rarefaction rate
  too_low <- which(rowSums(S16_reads) < rarefaction_rate)

  #too_low <- which(rowSums(S16_reads) < 10)
    if (!rlang::is_empty(too_low)) {
    messageColour("IQI can not be predicted for the following samples:\n",
                  "warning")
    messageColour(names(too_low),"warningSamples")
    messageColour("\nRead count under ", "warning")
    messageColour(rarefaction_rate, "warning")
    S16_reads <- S16_reads[!rownames(S16_reads)%in%names(too_low),]
    #identical(names(S16_reads2),names(S16_reads))#=TRUE, sense check.
    }

    #S16BU=S16_reads#code_testing
    #S16_reads=S16BU#code_testing
    #if(nrow(S16_reads==0)){message("TW. There are no samples >rarefaction limit so breaking");return(NULL)}
  #browser()

   if(nrow(S16_reads)==0){message("\nNo samples exceeded rarefaction limit, breaking");return(NULL)}
  # Rarefy data
  #S16_reads=Output_Full_Sepa
    rare_data <- vegan::rrarefy(S16_reads, rarefaction_rate)#output is a matrix
    rare_data <- as.data.frame(rare_data)

  # Clean sample names
  rownames(rare_data)=sub("_R1.*", "", rownames(rare_data),ignore.case = TRUE)

  # Restrict data to bacterial reads that are used by the random forest
  # object class changes to numeric if only one matching taxa (i.e., one column)
  # so check and break if necessary (only applies to debugging data)
  if(ncol(rare_data)<2){message("\nOnly one or fewer taxa match RF taxa, breaking");return(NULL)}

  rare_data <- rare_data[,c( colnames(rare_data) %in% RF_final_reduced$coefnames)]

  # Check if taxa used for RF predictions are present in any sample, if not add
  # them as zero values

  if (FALSE %in% (RF_final_reduced$coefnames %in% colnames(rare_data))) {
    name <- RF_final_reduced$coefnames[!(RF_final_reduced$coefnames
                                         %in% colnames(rare_data))]
    ##this is easier code to understand
    #name2=setdiff(RF_final_reduced$coefnames,colnames(rare_data))
    #identical(name,name2)#check that output identical (sense check for suggested alteration).
    # Warn about missing taxa
    #if(!is.null(name)){- then complete the warning messages,same code as below}

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

# Check if taxa used for RF predictions are present in individual samples, warn
  # about missing taxa making IQI less reliable
    #check this code#############===
  #need to define name here too ,in case not called above
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
  #move prediction to column 1
  rare_data=rare_data[, c("predicted_IQI", setdiff(names(rare_data), "predicted_IQI"))]
  rare_data2<<-rare_data
  return(rare_data)
}
