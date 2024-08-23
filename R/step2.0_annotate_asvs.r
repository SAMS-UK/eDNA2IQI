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

  #instate parameteroptions
 if (missing(parameteroptions)) {
    utils::data("parameteroptions", envir = environment())
  }

####Check Folder name and add / at the end if its not there

if (substr(folder, nchar(folder), nchar(folder)) != "/") {
  folder <- paste0(folder, "/")
}
  ####1.Get reference database if missing####
  #Check for correct reference database, download if missing

  if (!file.exists(file.path(find.package("eDNA2IQI"), "extdata", "referenceDatabase_full.fa.gz"))) {

    # Prompt the user for permission to download the file
    user_input <- readline(prompt = "The reference database file that needs to be stored locally is missing. Do you want to attempt to download it? (y/n): ")

    if (tolower(user_input) == "y") {

    tryCatch({
      # Attempt to download the file from SAMS Thredds server
      url <- "https://thredds.sams.ac.uk/thredds/fileServer/Full_database/referenceDatabase_full.fa.gz"
      filePath <- file.path(find.package("eDNA2IQI"), "extdata", "referenceDatabase_full.fa.gz")
      options(timeout = 300)
      messageColour("Updating taxa reference database from SAMS Thredds server \n\n", "warnMessage")
      utils::download.file(url, filePath, quiet = TRUE)
      message("SAMS Thredds server download successful")
    }, error = function(e) {
      # Handle error if SAMS Thredds download fails
      message("\nSAMS Thredds server download unsuccessful \n\n", conditionMessage(e))
      message("\n\n Attempting download from Zenodo, may be unstable in the future")

      tryCatch({
        # Attempt to download the file from Zenodo erver
        url <- "https://zenodo.org/records/4587955/files/silva_nr99_v138.1_train_set.fa.gz?download=1"
        filePath <- file.path(find.package("eDNA2IQI"), "extdata", "referenceDatabase_full.fa.gz")
        options(timeout = 300)
        messageColour("Updating taxa reference database from Zenodo server \n\n", "warnMessage")
        utils::download.file(url, filePath, quiet = TRUE)
        message("Zenodo server download successful")
      }, error = function(e) {
        # Handle error if Zenodo download fails
        message("\nZenodo server download unsuccessful \n\n", conditionMessage(e))
        messageColour("\n
                    Automatic downloads have failed \n
                    Try manual download from: \n
                    https://thredds.sams.ac.uk/thredds/fileServer/Full_database/referenceDatabase_full.fa.gz, \n
                    or \n
                    https://zenodo.org/records/4587955/files/silva_nr99_v138.1_train_set.fa.gz \n\n
                    and save file as `referenceDatabase_full.fa.gz` \n
                    in the /extdata/ directory in the eDNA2IQI R library files", "message")
        message("\n")
      })
    })
  } else {
    stop(messageColour("Reference database needs to be downloaded to continue with taxa allocation. \n
          Proceed with automatic download, or manual download from: \n
                    https://thredds.sams.ac.uk/thredds/fileServer/Full_database/referenceDatabase_full.fa.gz, \n
                    or \n
                    https://zenodo.org/records/4587955/files/silva_nr99_v138.1_train_set.fa.gz \n\n
                    and save file as `referenceDatabase_full.fa.gz` \n
                    in the /extdata/ directory in the eDNA2IQI R library files \n", "message"))
  }


  if (!file.exists(file.path(find.package("eDNA2IQI"), "extdata", "referenceDatabase_full.fa.gz"))) {
    # If the file does not exist, stop execution with an error message
    stop()
  }
}

  #continue
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

  save(S16_reads, file = file.path(folder,"outputData/taxaAllocatedReads_", taxalevel, ".rda", fsep = ""))


  messageColour("Function finished: 'taxa_allocation' \n Taxa allocated reads written to file \n\n", "message")

 ### Generate raw read taxa plots

	#shorten sample name
	S16_readsB <- S16_readsB %>%
	dplyr::mutate(SampleID = stringr::str_remove(SampleID, "_S.*"))

	#long format
	ST <- tidyr::pivot_longer(S16_readsB,cols=grep("Bact",colnames(S16_readsB)),
                 names_to = "Taxon",
                 values_to = "reads")

	#get total abundance
	total_abundance <- ST %>%
		dplyr::group_by(Taxon) %>%
		dplyr::summarise(Total_Abundance = sum(reads), .groups = 'drop') %>%
		dplyr::arrange(desc(Total_Abundance))

#get top 20
	top_20_taxa <- total_abundance %>%
		dplyr::top_n(20, Total_Abundance) %>%
		dplyr::pull(Taxon)

#collate non top 20 into others
	ST <- ST %>%
		dplyr::mutate(Taxon = ifelse(Taxon %in% top_20_taxa, Taxon, "Other_Bacteria"))

#collate others and order by sample
	reorder_data <- ST %>%
		dplyr::group_by(SampleID,Taxon) %>%
		dplyr::summarise(Total_Abundance = sum(reads), .groups = 'drop') %>%
		dplyr::arrange(SampleID,desc(Total_Abundance))

#add "Others" onto top_20_taxa vector
	top_20_taxa <- c(top_20_taxa, "Other_Bacteria")


#reorder
reordered_data <- reorder_data %>%
  dplyr::mutate(Taxon = factor(Taxon, levels = top_20_taxa)) %>%
  dplyr::arrange(SampleID, Taxon)


#make colour palette
palette1 <- RColorBrewer::brewer.pal(n = 12, "Set3")
palette2 <- RColorBrewer::brewer.pal(n = 9, "Paired")
custom_palette <- c(palette1, palette2)

#ggplot barplot
taxaplot <- ggplot2::ggplot(reordered_data, ggplot2::aes(x = SampleID, y = Total_Abundance, fill = Taxon)) +
 ggplot2::geom_bar(stat = "identity", width = 0.95) +
  ggplot2::scale_fill_manual(values = custom_palette) +
  ggplot2::scale_y_continuous(expand = c(0,0))+
  ggplot2::theme_classic() +
  ggplot2::theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.direction = "vertical",
    legend.title = ggplot2::element_text(size = 10),
    legend.text = ggplot2::element_text(size = 8),
    plot.margin = ggplot2::margin(1, 1, 2, 1, "cm"),
    axis.text.x = ggplot2::element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5)
  ) +
  ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1)) +
  ggplot2::labs(title = "Raw Read Counts of Top 20 Taxa and 'Other Bacteria'",
       x = "Sample",
       y = "Read Count",
       fill = "Taxa")

ggplot2::ggsave(filename = file.path(folder,"/outputData/raw_read_taxaplot.png"), plot = taxaplot, height = 8, units = "in")

messageColour("Raw Taxa Plot saved in Output folder \n\n", "message")

return(S16_reads)
}

