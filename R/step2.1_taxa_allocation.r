#' Convert sequence reads to bacteria taxa names using the SILVA v138.1 database
#' as a reference.
#'
#' @param folder A string. Location of raw data files.
#' @param collated_asv_batches A dataframe. Contains the S16 ASV reads, collated across batches.
#' @param taxalevel A string. Specifies the taxonomic level used for analysis.
#'
#' @export
#' @return A dataframe.
#'
step2.1_taxa_allocation <- function(folder, collated_asv_batches, taxalevel) {
  message("Function: Step2.1_taxa_allocation")
 
  

  # Use reference database to allocate taxa#
  reference <- file.path(find.package("eDNA2IQI"),
                         "extdata/referenceDatabase_full.fa.gz")

  taxa <- as.data.frame(dada2::assignTaxonomy(as.matrix(collated_asv_batches), reference,
                                             multithread = TRUE))

  # Remove categories that are not part of the taxon
  #tomw edit, 16 June 2024
  ReplaceInsertae=function(x){ifelse(grepl("incertae|unknown",x,ignore.case = TRUE),NA,x)}
  taxa[] <- lapply(taxa, function(x) ReplaceInsertae(x))
  taxa[] <- lapply(taxa, function(x) stringr::str_replace_all(x, ";", ""))
  taxa[] <- lapply(taxa, function(x) stringr::str_replace_all(x, " ", ""))
  taxa[] <- lapply(taxa, function(x) stringr::str_replace_all(x, "\\.", ""))
  taxa[] <- lapply(taxa, function(x) stringr::str_replace_all(x, "-", ""))
  #tomw edit end. 
                   
  # removes leading and tailing white_spaces
  taxa[] <- lapply(taxa, trimws, which = "both")
  # Set NA values to "Xunidentified"
  taxa <- as.data.frame(taxa)
  taxa[is.na(taxa)] <- "XUnidentified"

  # Add the taxa together
  taxa <- within(taxa, {
    Phylum <- paste(Kingdom, Phylum, sep = "_")
    Class <- paste(Phylum, Class, sep = "_")
    Order <- paste(Class, Order, sep = "_")
    Family <- paste(Order, Family, sep = "_")
    Genus <- paste(Family, Genus, sep = "_")
    Species <- paste(Family, Species, sep = "_")
  })

  # Merge ASV reads and bacteria, set taxa level specificity as column names,
  # remove other taxa
  taxa_combined <- rbind(collated_asv_batches, t(taxa))
  colnames(taxa_combined) <- taxa_combined[taxalevel, ]
  taxa_combined <- as.data.frame(taxa_combined[-which(rownames(taxa_combined)
                                                      %in% colnames(taxa)), ])
  taxa_combined <- data.frame(lapply(taxa_combined, as.numeric),
                              check.names = FALSE,
                              row.names = rownames(taxa_combined))

  # Remove columns with chloroplast or mitochondria in the name
  taxa_combined <- as.data.frame(t(rowsum(t(taxa_combined),
                                          group = colnames(taxa_combined),
                                          na.rm = TRUE)))
  # Sum the reads in columns with the same names

  taxa_combined <- taxa_combined %>%
    dplyr::select(-dplyr::contains("chloroplast")) %>%
    dplyr::select(-dplyr::contains("mitochondria"))

  S16_reads <- taxa_combined[, grepl("Bacteria", names(taxa_combined))]
  # Only keep Bacteria columns


  message("Function Step2.1_taxa_allocation finished")

  return(S16_reads)
}


                   
