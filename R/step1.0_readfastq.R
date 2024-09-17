#' Read raw fastq.gz files, clean sequence reads, allocate taxa at level
#' specified. Save dataframe of S16 reads (as .rda), save dataframe of taxa
#' allocated reads (as .csv).
#'
#' @param folder A string. Location of raw fastq data files.
#' @param parameteroptions A dataframe. Values for all function parameters.
#'
#' @export
#' @return A data frame. Table of reads per taxa, sorted by sample.
#'
#' @section Example usage:
#' myASVs <- step1.0_readFastq(folder = "C:/full/file/path/to/fastq/", parameteroptions)
#'

step1.0_readFastq=function(folder, parameteroptions) {
  message("Function: step1.0_readFastq")

####get wd to return to at end
	return_to <- getwd()
####Check Folder name and add / at the end if its not there

if (substr(folder, nchar(folder), nchar(folder)) != "/") {
  folder <- paste0(folder, "/")
}

#### 0.5 Run step0.5capitaliseR1_R2 as default. ####

  step0.5_capitaliseR1_R2(folder)
  message("File name syntax checked, and corrected where necessary")


#### 1.Load parameteroptions dataframe if missing from function call ####

  if (missing(parameteroptions)) {
    utils::data("parameteroptions", envir = environment())#
  }

# Get taxalevel
  taxalevel=extractParameters(parameteroptions,"global","taxalevel")

#### 2.Check taxa level and generate folders ####

taxas <- c("Family", "Genus", "Species", "Order", "Class", "Phylum", "Kingdom")
  if (!taxalevel %in% taxas) {
    messageColour("Taxanomic level specified not valid \n
Choose from: \nKingdom, \nPhylum, \nClass, \nOrder, \nFamily, Genus, \nSpecies \n\nTaxa level is case-sensitive. \n\n", "warning")
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    stop()
  }

#### 3.Generate data folder to store results and intermediate folder to store ####
# intermediate step files
  outputfolder=paste(folder,"outputData",sep="")
  if (!dir.exists(file.path(outputfolder))) {dir.create(file.path(outputfolder)) }
  intermediatefolder=paste(outputfolder,"/intermediate",sep="")
  if (!dir.exists(intermediatefolder)) {dir.create(intermediatefolder)}
#### 4.Print parameters to file ####
  dataCleanfile <- file.path(folder, "outputData/dataCleaningDetails.txt")
  sink(dataCleanfile, append = TRUE)
  cat("Step 1 Start Date/Time:",as.character(Sys.time()))
  cat("\n\n")
  cat("eDNA2IQI Package Version:", as.character(utils::packageVersion("eDNA2IQI")))
  cat("\n\n")
  cat("Parameter options used: \n")
  suppressWarnings(utils::write.table(parameteroptions[c(1:3)], file = dataCleanfile,
                    row.names = FALSE, quote = FALSE, sep = "\t", append = TRUE,col.names = TRUE))
  cat("\n\n")
  sink()

#### 5.Check raw data quality - call Step_quality_check####
  doCheck=extractParameters(parameteroptions,"global","qualityCheck","Char2Vect")
  message("Quality checking raw data, Set qualityCheck [3,3] in parameteroptions to 'FALSE' to disable")
  if (doCheck) {step1.01_quality_check(folder)#outputs quality score graphics (no return)
  }
####


#### 6.extract data frame of unique instrument_run identifier, sort in batches ####
  #batches are separate members of list.
  sample_IDs=step1.1_extract_identifier(folder)

#### 7.MAIN LOOP denoise different batches separately, output CSV file, per batch ####
  #note that this passes batches of files, not individual files.
	set.seed(123)

    for (i in seq_along(sample_IDs))
    {
      #i=1#code testing
    (batch_name <- paste("ASV_reads_batch_", i, ".csv", sep = ""))
    (batch_name2 <-paste("ASV_reads_batch_", i, ".rda", sep = ""))
    mesBatch <- paste("Processing batch #", i, "\n", sep = "")
    messageColour(mesBatch, "message")
    sink(dataCleanfile, append = TRUE)
    cat("########################",batch_name,"#########################\n\n", sep = "")
    sink()

# Check if batch.csv already created, and skip if TRUE
    if (file.exists(file.path(paste(folder, "outputData/", batch_name, sep="")))) {
    #this will error if additional samples to the same batch are copied into the test folder, without resetting
    #could error check this, based on sample names comparison.
    #if samples are added, then the whole batch would need to be redone.

     messageColour("Data already denoised, delete 'intermediate' folder to repeat, \n\n", "warnMessage")

        } else ##################################### main function calls here:
          #calls to generate NS and trimmed and filtered.
        {
      #Load raw files, remove Ns, and save files into 'intermediate' folder
      #sample_IDs[[i]] - a list of sample names, per batch
      ####Step1.2 remove unspecified ####

      no_unspec <- step1.2_remove_unspecified(folder, sample_IDs[[i]])
      fwd_no_unspec <- no_unspec[[1]]; rev_no_unspec <- no_unspec[[2]]

      # Trim primer sequences from the files in intermediate folder

      ####Step1.3_trim_primers ####
      trimmed <- step1.3_trim_primers(folder, fwd_no_unspec, rev_no_unspec, parameteroptions)
      fwd_trimmed <- trimmed[[1]]; rev_trimmed <- trimmed[[2]]

      #### Step1.4 Filter, adjust with error model, merge, chimera check ####
      # loop cycles through sample_IDs list, batch-by-batch
      ASV_reads_batch=step1.4_filter_quality_chimera(folder, sample_IDs[[i]],
                                             fwd_trimmed, rev_trimmed,
                                             parameteroptions)

      #### Save dataframe of ASV reads ####
      if(inherits(ASV_reads_batch, "data.frame"))#error check, only writes if a dataframe returned
        {#batch_name is e.g. "ASV_reads_batch_1.csv";batch_name2 is the rda equivalent.
      utils::write.csv(ASV_reads_batch, file=file.path(folder, "outputData", batch_name),row.names = FALSE)
      save(ASV_reads_batch,      file=file.path(folder, "outputData", batch_name2))
        }
      }
      } #end of loop, writing a csv file of ASV reads to outputData folder, per batch of denoised samples, as CSV and rda files.
      ####


  CollateASVBatches=CollateASVBatches(folder, return_to, dataCleanfile)
  message("End of function: Step1.0_readFastq")
  return(CollateASVBatches)
}

CollateASVBatches=function(folder,return_to, dataCleanfile){
  message("Function: CollateASVBatches, for input to Step2.0")

  path=paste(folder,"outputData",sep="")
  asv_per_batch=list.files(path,pattern = "batch.*\\.rda$",full.names = TRUE)

  Combined_ASV_batches=data.frame()
  for (files_to_combine in asv_per_batch){
  A1=get(load(files_to_combine))
  Combined_ASV_batches=dplyr::bind_rows(A1,Combined_ASV_batches)
  }
  Combined_ASV_batches[is.na(Combined_ASV_batches)]=0
  #for user sense-check, sample names retained as column1.

  rownames(Combined_ASV_batches) <- Combined_ASV_batches$rowname
  Combined_ASV_batches$rowname <- NULL#keep all columns numerical.

  #Save final ASV reads dataframe, called 'compiled_asv_reads.rda'.
  save(Combined_ASV_batches, file = file.path(folder,"outputData","collated_asv_batches.rda"))

  #return to original wd
  setwd(return_to)


  sink(dataCleanfile, append = TRUE)
  cat("\n")
  cat("Step 1 End Date/Time:",as.character(Sys.time()))
  cat("\n\n")
  sink()


  message("Function: CollateASVBatches complete")
  return(Combined_ASV_batches)
}
