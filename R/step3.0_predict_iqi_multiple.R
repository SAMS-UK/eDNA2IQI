#' Generates the predicted IQI values based on the optimized random forest
#' models.
#'
#' @param folder A string. Location of raw fastq data files.
#' @param AnnotatedASVs A data frame. Contains the bacterial taxa read counts which
#' should be used to predict the IQI values.
#'
#' @importFrom utils glob2rx
#'
#' @return A data frame. As inputted, with added 'predicted_IQI' column
#' @export
#' @section Example usage: myPREDs <- step3.0_predict_iqi(folder = "C:/full/file/path/to/fastq/S16_reads", S16_reads = myTAXA)
#'
teststep3.0_predict_iqi_multiple <- function(folder, AnnotatedASVs) {
  
  messageColour("Function starting: 'predict_iqi_multiple'
  Predicting IQI using optimized random forest models  \n\n", "message")
  step3_start_time <- as.character(Sys.time())
  message("Step 3 Start Date/Time:",step3_start_time)
  
  
   ####Get model if missing####
  #Check for correct reference database, download if missing

  if (!file.exists(file.path(find.package("eDNA2IQI"), "extdata/", "RFModel_rarefy_5000_taxalevel_family_BMB_Run1_5_multiple.rda"))) {
    downloadExternal(auto_download = auto_download)

  }
   if (!file.exists(file.path(find.package("eDNA2IQI"), "extdata/", "RFModel_rarefy_5000_taxalevel_family_BMB_Run1_5_multiple.rda"))) {
      # If the file does not exist, stop execution with an error message
      stop()
    }
  
  ####Check Folder name and add / at the end if its not there
  if (substr(folder, nchar(folder), nchar(folder)) != "/") {
    folder <- paste0(folder, "/")
  }
  
  #rownames(AnnotatedASVs)=AnnotatedASVs$SampleID
  #AnnotatedASVs$SampleID=NULL
  
  #glob2rx is a utils::function
  #it is used to convert wildcard expressions (globs) to regular expressions.
  # Find trained random forests, extract rarefaction rate
  potenRfs <- dir(file.path(find.package("eDNA2IQI"), "extdata/"),
                  full.names = TRUE)
  
  # Load trained random forests
  RF=subset(potenRfs, grepl(glob2rx("*/*_BMB_Run*.rda"), potenRfs,ignore.case = TRUE))
  LOB=readRDS(RF);
  LOB=LOB[order(sapply(LOB,function(x) x[[2]]))]#order by RMSE
  RF=tolower(RF)
  rarefaction_rate <- as.numeric(stringr::str_extract(RF, "(?<=rarefy_)[0-9]+"))
  taxalevel <- as.character(stringr::str_extract(RF, "(?<=_taxalevel_)[^_]+"))
  

  # Check for samples  with fewer reads than rarefaction rate
  too_low <- which(rowSums(AnnotatedASVs) < rarefaction_rate)
  read_totals <- rowSums(AnnotatedASVs)
  readcounts_df <- data.frame("Denoised_Reads" = read_totals, row.names = rownames(AnnotatedASVs))
  
  AnnotatedASVs=AnnotatedASVs[rowSums(AnnotatedASVs)>rarefaction_rate,]
  if(nrow(AnnotatedASVs)==0){message("\nNo samples exceeded rarefaction limit, breaking");return(NULL)}
  
  #create frame as NULL, so it doesnt break when empty
  excluded_samples <- NULL
  
  if (!rlang::is_empty(too_low)) {
    
    # Store samples with low reads in a separate data frame
    excluded_samples <- data.frame(SampleID = names(too_low),
                                   MeanIQI = NA,
                                   stringsAsFactors = FALSE)
    row.names(excluded_samples) <- excluded_samples$SampleID
    excluded_samples$SampleID <- NULL
    
    messageColour("\nRead count under ", "warning")
    messageColour(rarefaction_rate, "warning")
    messageColour(" for the following samples. IQI was not calculated.\n", "warning")
    
    messageColour(paste(names(too_low), collapse = "\n"), "warningSamples")
    
    #remove samples under rarefaction
    AnnotatedASVs <- AnnotatedASVs[!rownames(AnnotatedASVs)%in%names(too_low),]
    
  }
  
  if(nrow(AnnotatedASVs)==0){
    message("\nNo samples exceeded rarefaction limit. Exiting function.")
    return(NULL)
  }
  
  set.seed(123)
  rare_data <- suppressWarnings(as.data.frame(vegan::rrarefy(AnnotatedASVs, rarefaction_rate)))
  for_barplot <- rare_data
  
  #write rarefied dataframe
  utils::write.csv(rare_data, file = file.path(folder,
                                               "/outputData/rarefied_taxa_allocated_reads_", taxalevel, ".csv", fsep = ""),row.names = TRUE)
  
  # Clean sample names
  #rownames(rare_data)=sub("_R1.*", "", rownames(rare_data),ignore.case = TRUE)
  
  
   ##make predictions
  

  #rounding function
  SIG=function(X){signif(X,4)}
  
  #LOOP THROUGH THE RFS IN THE LIST.
  CP2=vector()#to contain all taxa used in each model then combine, tabulate for output.
  collated_predictions=data.frame(matrix(nrow=nrow(rare_data),ncol=0))#to contain predicted IQI.
  
  #i=1
  for (i in 1:length(LOB)){
    RF_final_reduced=LOB[[i]][[1]]
    rare_data2=rare_data[,c(colnames(rare_data) %in%RF_final_reduced$xNames)]
    # Check if taxa used for RF predictions are present in any sample
    name2=setdiff(RF_final_reduced$xNames,colnames(rare_data2))
    # Add columns with zero values for missing taxa
    if(!is.null(name2)){rare_data2[name2]=0}
    
    #### Make predictions with the trained random forest####
    predicted_IQI=SIG(stats::predict(RF_final_reduced, rare_data2))
    NosFeatures=length(LOB[[i]][[1]]$xNames)
    collated_predictions=cbind(collated_predictions,predicted_IQI)
    CP=LOB[[i]][[1]]$xNames
    CP2=c(CP,CP2)
  }
  
  colnames(collated_predictions) <- paste0("RF_", 1:length(LOB))
  collated_predictions$MeanIQI  =SIG(rowMeans(collated_predictions))
  collated_predictions$SdevIQI  =SIG(apply(collated_predictions[1:length(LOB)], 1, sd, na.rm=TRUE))
  collated_predictions$MaxIQI   =SIG(apply(collated_predictions[1:length(LOB)], 1, max, na.rm=TRUE))
  collated_predictions$MinIQI   =SIG(apply(collated_predictions[1:length(LOB)], 1, min, na.rm=TRUE))
  collated_predictions$RangeIQI =SIG(collated_predictions$MaxIQI-collated_predictions$MinIQI)
  collated_predictions$SampleID=rownames(AnnotatedASVs)
  
  A=grep("SampleID", colnames(collated_predictions))#bring SampleID to 1st column.
  B=setdiff(1:ncol(collated_predictions),A)
  collated_predictions=collated_predictions[,c(A,B)]
  
  # Clean sample names
  collated_predictions$SampleID=sub("_R1.*", "", collated_predictions$SampleID,ignore.case = TRUE)
  
  #save predictions
  
  if (!is.null(excluded_samples)) {
    
    final_data <- bind_rows(collated_predictions, excluded_samples)
  } else {
    final_data <- collated_predictions
  }
  
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
  
  ####tidy and exit
  AllTaxa=as.data.frame(table(CP2))#tabulate key taxa as utilised in multiple RF models.
  #write file
  utils::write.csv(AllTaxa, file = file.path(folder,
                                                "/outputData/AllTaxa.csv", fsep = ""),row.names = TRUE)
												
  return(list(Predicted_IQI=collated_predictions))
}
