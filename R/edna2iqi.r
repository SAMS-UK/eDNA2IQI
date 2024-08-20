#' Wrapper function to perform all steps of Reading, denoising,
#' allocating taxa and predicting IQI
#'
#'
#' @param folder A string. Location of raw fastq data files.
#' @param parameteroptions A dataframe. Values for all function parameters.
#'
#' @export
#' @return dataframes containing denoised ASVs (myASVs),
#'                               allocated taxa (myTaxa) and
#'                               predicted IQI values (myPreds)
#'
#' @section Example usage:
#' edna2iqi(folder = "C:/full/file/path/to/fastq/", parameteroptions = parameteroptions)
#'

edna2iqi=function(folder, parameteroptions) {

#instate parameteroptions
 if (missing(parameteroptions)) {
    utils::data("parameteroptions", envir = environment())
  }

  message("Wrapper function: edna2iqi")

  myASVs <- step1.0_readFastq(folder, parameteroptions)

  myTaxa <- step2.0_annotate_ASVs(folder,parameteroptions, myASVs)

  myPreds <- step3.0_predict_iqi(myTaxa)

  message("edna2iqi finished")

  return(list(myASVs = myASVs, myTaxa = myTaxa, myPreds = myPreds))
}
