#' Title Step2.0_annotate_ASVs
#'
#' @param folder A string. Location of raw fastq files
#' @param parameteroptions A dataframe,
#' @param collated_asv_batches result from Step1..., dataframe
#' @param auto_download Boolean to override console prompt to automatically
#'   download external dependencies. This can be useful for automating or
#'   testing this function.
#'
#' @return A dataframe, written to file, combing ASVs across multiple batches
#' @export
#'
#' @section Example usage: myTAXA <- step2.0_annotate_ASVs(folder,parameteroptions,myASVs)

step2.0_annotate_ASVs=function(folder,parameteroptions,collated_asv_batches, auto_download = FALSE)
{
  message("Function: Step2.0_annotate_ASV")
  message("rda used in taxa_allocation call")

  #instate parameteroptions
  if (missing(parameteroptions)) {
    utils::data("parameteroptions", envir = environment())
  }

  #setFolder name and add / at the end if its not there----
  if (substr(folder, nchar(folder), nchar(folder)) != "/") {
    folder <- paste0(folder, "/")
  }

   #create output folder
   outputfolder=paste(folder,"outputData",sep="")
  if (!dir.exists(file.path(outputfolder))) {dir.create(file.path(outputfolder)) }

  #1.Get reference database if missing----
  #Check for correct reference database, download if missing

  if (!file.exists(file.path(find.package("eDNA2IQI"), "extdata", "referenceDatabase_full.fa.gz"))) {
    browser()
    downloadExternal(auto_download = auto_download)

  }
   if (!file.exists(file.path(find.package("eDNA2IQI"), "extdata", "referenceDatabase_full.fa.gz"))) {
      # If the file does not exist, stop execution with an error message
      stop()
    }
  #continue
  taxalevel=extractParameters(parameteroptions,"global","taxalevel")

  #2. step2.1 allocate taxa----
  message("Calling Step2.1_taxa_allocation")
  S16_reads=step2.1_taxa_allocation(folder, collated_asv_batches, taxalevel)

  S16_readsB=S16_reads
  S16_readsB$SampleID=rownames(S16_readsB)
  #move SampleID to first column
  A=grep("SampleID",colnames(S16_readsB))
  S16_readsB=S16_readsB[,c(A,1:(A-1))]

  #2.1 Save dataframe----
  message(paste("Writing annotated data here:",folder))
  utils::write.csv(S16_readsB, file = file.path(folder,
                                                "/outputData/taxaAllocatedReads_", taxalevel, ".csv", fsep = ""),row.names = FALSE)

  save(S16_reads, file = file.path(folder,"outputData/taxaAllocatedReads_", taxalevel, ".rda", fsep = ""))


  messageColour("Function finished: 'taxa_allocation' \n Taxa allocated reads written to file \n\n", "message")

  #3. drawBarplot of unrarefied reads----
  #Generate unrarefied read taxa plots
  drawBarplot(S16_readsB, folder)

  return(S16_reads)
}
