#' Wrapper function to perform all steps of Reading, denoising,
#' allocating taxa and predicting IQI
#'
#'
#' @param folder A string. Location of raw fastq data files.
#' @param parameteroptions A dataframe. Values for all function parameters.
#' @param auto_download Boolean to override console prompt to automatically
#'   download external dependencies. This can be useful for automating or
#'   testing this function.
#'
#' @export
#' @return dataframes containing denoised ASVs (myASVs),
#'                               allocated taxa (myTaxa) and
#'                               predicted IQI values (myPreds)
#'
#' @section Example usage:
#' edna2iqi(folder = "C:/full/file/path/to/fastq/", parameteroptions = parameteroptions)
#'

edna2iqi=function(folder, parameteroptions, auto_download = FALSE) {

####Check Folder name and add / at the end if its not there

if (substr(folder, nchar(folder), nchar(folder)) != "/") {
  folder <- paste0(folder, "/")
}

#instate parameteroptions
 if (missing(parameteroptions)) {
    utils::data("parameteroptions", envir = environment())
  }

  message("Wrapper function: edna2iqi")

  if (file.exists(file.path(paste(folder, "outputData/collated_asv_batches.rda", sep="")))) {
    #this will error if additional samples to the same batch are copied into the test folder, without resetting
    #could error check this, based on sample names comparison.
    #if samples are added, then the whole batch would need to be redone.

    messageColour("Data already denoised, delete 'intermediate' folder and files to repeat, \n\n", "warnMessage")

  } else

  {
  myASVs <- step1.0_readFastq(folder, parameteroptions, auto_download = auto_download)
  }


  if (file.exists(file.path(paste(folder, "outputData/taxaAllocatedReads_Family.rda", sep="")))) {
    #this will error if additional samples to the same batch are copied into the test folder, without resetting
    #could error check this, based on sample names comparison.
    #if samples are added, then the whole batch would need to be redone.

    messageColour("Data already taxa allocated, delete 'taxaAllocatedReads_*' files to repeat, \n\n", "warnMessage")
    load(file.path(paste(folder, "outputData/collated_asv_batches.rda", sep="")))
    myASVs <- Combined_ASV_batches

  } else

  {
  load(file.path(paste(folder, "outputData/collated_asv_batches.rda", sep="")))
  myASVs <- Combined_ASV_batches

  myTaxa <- step2.0_annotate_ASVs(folder,parameteroptions, myASVs)
  }

  load(file.path(paste(folder, "outputData/taxaAllocatedReads_Family.rda", sep="")))
  myTaxa <- S16_reads
  myPreds <- step3.0_predict_iqi(folder, myTaxa)

  message("edna2iqi finished")

  return(list(myASVs = myASVs, myTaxa = myTaxa, myPreds = myPreds))
}
