#' Title Step2.0_annotate_ASVs
#'
#' @param folder A string. Location of raw fastq files
#' @param parameteroptions A dataframe,
#' @param collated_asv_batches result from Step1..., dataframe
#'
#' @return A dataframe, written to file, combing ASVs across multiple batches
#' @export
#'
#' @section Example usage: myASVs <- step2.0_annotate_ASVs(folder,parameteroptions,collated_asv_batches)

step2.0_annotate_ASVs=function(folder,parameteroptions,collated_asv_batches)
{
  message("Function: Step2.0_annotate_ASV")
  message("rda used in taxa_allocation call")
  
  ####1.Get reference database if missing####
  #Check for correct reference database, download if missing
  
  if (!file.exists(file.path(find.package("eDNA2IQI"),#CORRECT VERSION, INSTALLS IF REQUIRED.
                             "/extdata/referenceDatabase_full.fa.gz"))) {
    url <- "https://thredds.sams.ac.uk/thredds/fileServer/Full_database/referenceDatabase_full.fa.gz"
    filePath <- file.path(find.package("eDNA2IQI"),
                          "/extdata/referenceDatabase_full.fa.gz", sep = "")
    options(timeout = 300)
    messageColour("Updating taxa reference database \n\n", "warnMessage")
    utils::download.file(url, filePath, quiet = TRUE)
  }

    taxalevel=extractParameters(parameteroptions,"global","taxalevel")
  
  message("Calling Step2.1_taxa_allocation")
  S16_reads=step2.1_taxa_allocation(folder, collated_asv_batches, taxalevel)#class(S16_reads)

  S16_readsB=S16_reads
  S16_readsB$SampleID=rownames(S16_readsB)
  A=grep("SampleID",colnames(S16_readsB))#move SampleID to first column
  S16_readsB=S16_readsB[,c(A,1:(A-1))]

  # Save dataframe
  message(paste("Writing annotated data here:",folder))
  utils::write.csv(S16_readsB, file = file.path(folder,
               "/outputData/taxaAllocatedReads_", taxalevel, ".csv", fsep = ""),row.names = FALSE)
  messageColour("Function finished: 'taxa_allocation' \n\n", "message")
  return(S16_reads)
}

