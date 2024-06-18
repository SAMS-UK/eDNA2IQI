#' Step1.4_filter_quality_chimera
#'
#' @param folder A string. Location of raw data files.
#' @param raw_files A string. List of files from same location.
#' @param fwd_trimmed A character vector. the file name of forward
#' reads of the input samples where primer sequences were removed.
#' @param rev_trimmed A character vector. Contains the file location of reverse
#' reads of the input samples where primer sequences were removed.
#' @param parameteroptions A dataframe.
#'
#' @export
#' @return A data frame.

#called by:
#raw_files are sample_IDs[[1]], these are batches of fastq on the same sequencing run, must be denoised collectively

Step1.4_filter_quality_chimera=function(folder, raw_files, fwd_trimmed, rev_trimmed,
                                   parameteroptions) {
 
  message("Function: Step1.4_filter_quality_chimera")
  message("generates ASVs from filtered files (in 'intermediate/filtered'")

  messageColour("Function starting: 'filter_quality_chimera' Removing low quality sequences \n\n", "message")
  #### 0.5 Write output and create files ####
  # Set file to store details of cleaning and filtering

  dataCleanfile <- file.path(folder, "outputData/dataCleaningDetails.txt")
  sink(dataCleanfile, append = TRUE)
  cat(as.character(Sys.time()))
  cat("\n\nFunction 1.4: filter_quality_chimera \n\n")
  sink()

  # Change folder of trimmed to filtered
  path_filtered <- file.path(folder, "outputData/intermediate", "filtered")

  

  if (!dir.exists(path_filtered)) {
    dir.create(path_filtered)
  }

  #create file_names in which to store the dada2 filter output.
  fwd_filtered <- sub("trimmed", "filtered", fwd_trimmed)#ok, replaces trimmed for filtered in string.
  rev_filtered <- sub("trimmed", "filtered", rev_trimmed)

  #### 1.Get dada2 filter parameters ####
  message("1.Get dada2 filter parameters")
  # Get parameter options


    maxn=ExtractParameters(parameteroptions,"filterAndTrim","maxN","Char2Vect")
    E1_E2=ExtractParameters(parameteroptions,"filterAndTrim","maxEE","Char2Vect")
    row_tL <- ExtractParameters(parameteroptions,"filterAndTrim","truncLen","Char2Vect")
    rp <- ExtractParameters(parameteroptions,"filterAndTrim","rm.phix","Char2Vect")
    tQ=ExtractParameters(parameteroptions,"filterAndTrim","truncQ","Char2Vect")
    minL=ExtractParameters(parameteroptions,"filterAndTrim","minLen","Char2Vect")
    maxL=ExtractParameters(parameteroptions,"filterAndTrim","maxLen","Char2Vect")
    c=ExtractParameters(parameteroptions,"filterAndTrim","compress","Char2Vect")
    idf=ExtractParameters(parameteroptions,"filterAndTrim","id.field","Char2Vect")
    minQ=ExtractParameters(parameteroptions,"filterAndTrim","minQ","Char2Vect")
    matchID=ExtractParameters(parameteroptions,"filterAndTrim","matchIDs","Char2Vect")
    multi=ExtractParameters(parameteroptions,"filterAndTrim","multithread","Char2Vect")
    n=ExtractParameters(parameteroptions,"filterAndTrim","n","Char2Vect")
    OMP=ExtractParameters(parameteroptions,"filterAndTrim","OMP","Char2Vect")
    orient=ExtractParameters(parameteroptions,"filterAndTrim","orient.fwd","Char2Vect")
    qtype=ExtractParameters(parameteroptions,"filterAndTrim","qualityType")
    rmlow=ExtractParameters(parameteroptions,"filterAndTrim","rm.lowcomplex","Char2Vect")
    tleft=ExtractParameters(parameteroptions,"filterAndTrim","trimLeft","Char2Vect")
    tright=ExtractParameters(parameteroptions,"filterAndTrim","trimRight","Char2Vect")
    v=ExtractParameters(parameteroptions,"filterAndTrim","verbose","Char2Vect")


  #### 2.Run dada2 filterAndTrim ####
  message("2.Run dada2 filterAndTrim")
  #fwd_trimmed file exists,
  #fwd_filtered path exists, but not files
  ####data2 creates fwd_filtered files.
  dada2::filterAndTrim(fwd_trimmed, fwd_filtered, rev_trimmed, rev_filtered,
                       maxN=maxn, maxEE=E1_E2, truncLen=row_tL,
                       rm.phix=rp, truncQ=tQ, compress=c, minLen=minL,
                       maxLen = maxL, id.field = idf, matchIDs = matchID,
                       minQ = minQ, multithread = multi, n = n, OMP = OMP,
                       orient.fwd = orient, qualityType = qtype,
                       rm.lowcomplex = rmlow, trimLeft = tleft,
                       trimRight = tright, verbose = v)

   # Calculate dada2 error rates
  #### 3.Learn-error parameters ####
  message("3.Get learn-error parameters")
  (nB=ExtractParameters(parameteroptions,"learnErrors","nbases","Char2Vect"))
  (rnd=ExtractParameters(parameteroptions,"learnErrors","randomise","Char2Vect"))
  (multi=ExtractParameters(parameteroptions,"learnErrors","multithread","Char2Vect"))
  (mc=ExtractParameters(parameteroptions,"learnErrors","MAX_CONSIST","Char2Vect"))
  (oc=ExtractParameters(parameteroptions,"learnErrors","OMEGA_C","Char2Vect"))
  (qtype=ExtractParameters(parameteroptions,"learnErrors","qualityType"))
  (v=ExtractParameters(parameteroptions,"learnErrors","verbose","Char2Vect"))

 
  #re-extract sample names in case samples have been removed (no reads passing filter)
  fwd_actual <- stringr::str_subset(dir(path_filtered), "_R1_001.fastq.gz")
  rev_actual <- stringr::str_subset(dir(path_filtered), "_R2_001.fastq.gz")
  fwd_expected <- stringr::str_subset(raw_files, "_R1_001.fastq.gz")
  rev_expected <- stringr::str_subset(raw_files, "_R2_001.fastq.gz")
  fwd_filtered <- intersect(fwd_actual, fwd_expected)
  fwd_filtered <- file.path(path_filtered, stringr::str_subset(fwd_filtered,
                                                               "_R1_001.fastq.gz"))
  rev_filtered <- intersect(rev_actual, rev_expected)
  rev_filtered <- file.path(path_filtered, stringr::str_subset(rev_filtered,
                                                               "_R2_001.fastq.gz"))


  #### 4.Run error model, object created ####
  message("4.Run error model, object created")
  message("forward error:")
  message("This is not memory intenstive")
  error_fwd <- dada2::learnErrors(fwd_filtered, multithread = multi,
                                  nbases = nB, randomize = rnd,
                                  verbose = v,
                                  MAX_CONSIST = mc, OMEGA_C = oc,
                                  qualityType = qtype)
  message("reverse error:")
  error_rev <- dada2::learnErrors(rev_filtered, multithread = multi,
                                  nbases = nB, randomize = rnd,
                                  verbose = v,
                                  MAX_CONSIST = mc, OMEGA_C = oc,
                                  qualityType = qtype)

  #### 5.Get dada2 inference parameters ####
  message("5.Get dada2 inference parameters")
  # Get parameter options

  (eEF=ExtractParameters(parameteroptions,"dada","errorEstimationFunction"))
  (multi=ExtractParameters(parameteroptions,"dada","multithread","Char2Vect"))
  (v=ExtractParameters(parameteroptions,"dada","verbose","Char2Vect"))
  (selfC=ExtractParameters(parameteroptions,"dada","selfConsist","Char2Vect"))
  (pool=ExtractParameters(parameteroptions,"dada","pool","Char2Vect"))


  #### 6.Run error model on the filtered reads, in the filtered directory ####
  message("6.Run error model on the filtered reads, in the filtered directory")
  message("forward error being processed - not memory intensive (~12GB")

  dada2_fwd <- dada2::dada(fwd_filtered, err = error_fwd, multithread = multi,
                           verbose = v, selfConsist = selfC, pool = pool)
  message("reverse error being processed")
  dada2_rev <- dada2::dada(rev_filtered, err = error_rev, multithread = multi,
                           verbose = v, selfConsist = selfC, pool = pool)

  #### 7.Write detail of cleaning #####
  message("7.Write detail of cleaning")
  dataCleanfile <- file.path(folder, "outputData/dataCleaningDetails.txt")
  fw_in_L <- c()
  fw_out_L <- c()
  sink(dataCleanfile, append = TRUE)
  cat(as.character(Sys.time()))
  

  #determine the difference in reads between trimmed and filtered.
  #and report
  for (x in seq_along(fwd_filtered)) {
    fw_in_L[x] <- length(ShortRead::id(ShortRead::readFastq(fwd_trimmed[x])))
    fw_out_L[x] <- length(ShortRead::id(ShortRead::readFastq(fwd_filtered[x])))

  cat("\nFiles ", sub(".*/trimmed/", "", sub(".fastq.gz", "", fwd_trimmed[x])),#clean-up file name
        " & ", sub(".*/trimmed/", "", sub(".fastq.gz", "", rev_trimmed[x])),
        " - ", fw_in_L[x], " reads quality filtered to ", fw_out_L[x],
        " reads (", round(100 * (fw_out_L[x] / fw_in_L[x]), 1), "%)",#output % reads passing filtration
        sep = "")
  }
  cat("\n\n")
  sink()

  # DATA FLAG: flag if data loss >70%
  names <- sub(".*/trimmed/", "", sub("_R1.*", "", fwd_trimmed))
  if (length(names[(fw_out_L / fw_in_L) < 0.3]) > 0) {#count how many samples lose more than 70%, if 1 or more then message
    messageColour("High loss of reads when removing low quality sequences, treat IQI predictions with caution.
                  \n\n", "warning")
    sink(dataCleanfile, append = TRUE)
    cat("High loss of reads when removing low quality sequences, data of low ")
    cat("quality.\nThe following samples lost >70% of reads.\nTreat IQI ")
    cat("predictions with caution.\n")
    for (y in names[(fw_out_L / fw_in_L) < 0.3]){
      cat(y, "\n")
    }
    cat("\n\n")
    sink()
  }

  #### 8.Merge forward and backward reads ####
  message("8.Merge forward and backward reads")
  # Get parameter options
  (MinO=ExtractParameters(parameteroptions,"mergePairs","minOverlap","Char2Vect"))
  (maxMis=ExtractParameters(parameteroptions,"mergePairs","maxMismatch","Char2Vect"))
  (propC=ExtractParameters(parameteroptions,"mergePairs","propagateCol","Char2Vect"))
  (rR=ExtractParameters(parameteroptions,"mergePairs","returnRejects","Char2Vect"))
  (v=ExtractParameters(parameteroptions,"mergePairs","verbose","Char2Vect"))
  (jC=ExtractParameters(parameteroptions,"mergePairs","justConcatenate","Char2Vect"))
  (tO=ExtractParameters(parameteroptions,"mergePairs","trimOverhang","Char2Vect"))

  #######merging reads
  message("merging reads")
  merged_reads <- dada2::mergePairs(dada2_fwd, fwd_filtered, dada2_rev,
                                    rev_filtered, verbose = v,
                                    minOverlap = MinO, maxMismatch = maxMis,
                                    propagateCol = propC, returnRejects = rR,
                                    justConcatenate = jC, trimOverhang = tO)

  # Construct ASV reads table and remove chimera reads
  
  message("making sequence table- very fast ")
  ASV_table_raw <- dada2::makeSequenceTable(merged_reads)
  #####write output to dataCleanfile ####
  fw_all <- c()
  fw_unique <- c()
  sink(dataCleanfile, append = TRUE)
  cat(as.character(Sys.time()))

  for (x in seq_along(fwd_filtered)) {
  
    fw_all[x]    <- sum(ASV_table_raw[x, ])
    fw_unique[x] <- sum(ASV_table_raw[x, ] > 0)

cat("\nFiles ", sub(".*/filtered/", "", sub(".fastq.gz", "", fwd_filtered[x])),
        " & ", sub(".*/filtered/", "", sub(".fastq.gz", "", rev_filtered[x])),
        " - ", fw_all[x], " reads, in ", fw_unique[x], " unique sequences",
        sep = "")
  }
  cat("\n\n")
  sink()

  # DATA FLAG: flag if data loss >70%
  if (length(names[(fw_all / fw_out_L) < 0.3]) > 0) {
    messageColour("High loss of reads when merging forward and backward reads,
                  data of low quality. Treat IQI predictions with caution.\n\n", "warning")
    sink(dataCleanfile, append = TRUE)
    cat(as.character(Sys.time()))
    cat("\nHigh loss of reads when merging forward and backward reads")
    cat("\nThe following samples lost >70% of reads.\nTreat IQI predictions with caution.\n")
    for (y in names[(fw_all / fw_out_L) < 0.3]){
      cat(y, "\n")
    }
    cat("\n\n")
    sink()
  }


  #### 9.Remove chimeras and print cleaning details ####
  message("9.Remove chimera read and print cleaning details")
  # Get parameter options
  (mth=ExtractParameters(parameteroptions,"removeBimeraDenovo","method"))
  (v=ExtractParameters(parameteroptions,"removeBimeraDenovo","verbose","Char2Vect"))

  #ASV_table_raw created by dada2 in above step ('makesequencetable').
  message("very CPU intensive >100%, note that multithread=TRUE, but fast step")
  ASV_table_no_chim <- dada2::removeBimeraDenovo(ASV_table_raw,
                                                 method = mth,
                                                 multithread = TRUE,
                                                 verbose = v)
  
  #####write output to dataCleanfile ####

  fw_all2 <- length(colnames(ASV_table_raw))
  fw_no_chim <- length(colnames(ASV_table_no_chim))

  sink(dataCleanfile, append = TRUE)
  cat(as.character(Sys.time()))
  cat("\n", fw_all2 - fw_no_chim, " sequences removed from ", fw_all2, " (",
  round(100 * (fw_all2 - fw_no_chim) / fw_all2, 1), "%), as chimeras, ")
  cat("across all samples in batch. \n\n\n", sep = "")
  sink()

  # DATA FLAG: flag if data loss >70%
  if (length(names[(rowSums(ASV_table_no_chim) / rowSums(ASV_table_raw)) < 0.3])
      > 0) {
   messageColour("High loss of reads when removing chimeras, data of low quality.
                  Treat IQI predictions with caution.\n\n", "warning")
    sink(dataCleanfile, append = TRUE)
    cat(as.character(Sys.time()))
    cat("\nHigh loss of reads when removing chimeras, data of low quality.")
    cat("\nThe following samples lost >70% of reads.\nTreat IQI ")
    cat("predictions with caution\nTW_needs debugging.\n")
    
    for (y in names[(rowSums(ASV_table_no_chim) / rowSums(ASV_table_raw)) < 0.3]){

      cat(y, "\n")
    }
    cat("\n\n")
    sink()
  }
  



  C=list.files(paste(folder,"outputData/intermediate/filtered/",sep=""))

 

  Retained_samples=rownames(ASV_table_no_chim)
  fwd_input2=file.path(folder,Retained_samples)
  fwd_input3=gsub("//", "/", fwd_input2)#shortcut, in case of folder/ or folder specified.

  fw_in_L <- c()#empty variable.  Vector of total counts for samples.
  for (x in seq_along(fwd_input3)) {
    fw_in_L[x] <- length(ShortRead::id(ShortRead::readFastq(fwd_input3[x])))#gives the total read count for all the samples
  }
  # DATA FLAG: flag if data loss >70%
  

  if (length(names[(rowSums(ASV_table_no_chim) / fw_in_L) < 0.3]) > 0) {
    messageColour("High loss of reads during denoising, data of low quality.
                   Treat IQI predictions with caution.", "warning")
    sink(dataCleanfile, append = TRUE)
    cat(as.character(Sys.time()))
    cat("\nTW_High loss of reads during denoising, data of low quality.\n")
    cat("The following samples lost >70% of reads.\nTreat IQI predictions with caution.\n")

   
    A=as.data.frame(ASV_table_no_chim)#ASV_table_no_chim is a matrix_array.
    (B=(rowSums(A)/fw_in_L)<0.3)
    C=names(B)
    for (y in C)
         {cat(y, "\n")}

    
    cat("\n\n")
    sink()
  }

  # DATA FLAG: flag if <10000 reads in merged samples
  if (length(names[rowSums(ASV_table_no_chim) < 10000]) > 0) {
    messageColour("Low merged read count for following samples.\n", "warning")
    messageColour("Confidence in prediction may be low.
               Consider investigation of raw read count and denoising statistics. \n",
              "warning")
    for (y in names[rowSums(ASV_table_no_chim) < 10000]){
      messageColour(y, "warningSamples")
      messageColour("\n", "warningSamples")
    }
    cat("\n")

    sink(dataCleanfile, append = TRUE)
    cat(as.character(Sys.time()))
    cat("\nLow merged read count for following samples.\n")
    cat("Confidence in prediction may be low.\n")
cat("Consider investigation of raw read count and denoising statistics.\n")
    for (y in names[rowSums(ASV_table_no_chim) < 10000]){
      cat(y, "\n")
    }
    cat("\n\n")
    sink()
  }

  rownames(ASV_table_no_chim) <- gsub("^.*/", "", fwd_filtered)


  ASV_table_no_chim<-tibble::rownames_to_column(as.data.frame(ASV_table_no_chim))

  # DATA FLAG: flag if <25 ASVs
  if (length(colnames(ASV_table_no_chim)) < 25) {
    messageColour("Low number of unique sequence variants. Treat IQI predictions with caution. \n\n", "warning")
    sink(dataCleanfile, append = TRUE)
    cat(as.character(Sys.time()))
    cat("\nLow number of unique sequence variants. \n")
    cat("Treat IQI predictions with caution.\n\n\n")
    sink()
  }
  #####
  messageColour("Function finished: 'filter_quality_chimera' \n\n", "message")
  return(ASV_table_no_chim)#class is tibble or data.frame
  #this is a 'batch' of files, from the same sequencing run.
}


