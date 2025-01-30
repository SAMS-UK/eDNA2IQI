#' drawBarplot
#'
#' Takes S16_readsB from step2, and step3 and saves unrarefied and rarefied taxa barplot respectively
#'
#' @param S16_readsB annotated reads for plotting from
#' @param folder active folder
#' @return NULL
#' @importFrom stats median
#' @importFrom utils tail
#' @keywords internal

drawBarplot <- function(S16_readsB, folder) {

  #reduce names
  ShortenNames6 <- function(T) {
    # Remove everything from the first occurrence of "_XUnidentified" onwards
    cleaned_string <- sub("_XUnidentified.*", "", T)

    # Split the cleaned string by underscores into individual components
    parts <- unlist(strsplit(cleaned_string, "_"))

    # Ensure there are at least two parts to consider
    if (length(parts) >= 2) {
      # Get the last two parts before the first occurrence of "XUnidentified"
      last_two_parts <- tail(parts, n = 2)
      # Combine the last two parts with an underscore
      result <- paste(last_two_parts, collapse = "_")
    } else {
      # If there is only one part, just return that
      result <- parts[1]
    }

    return(result)
  }




#shorten sample name
	S16_readsB <- S16_readsB %>%
	dplyr::mutate(SampleID = stringr::str_remove(SampleID, "_S.*"))

	#long format
	ST <- tidyr::pivot_longer(S16_readsB,cols=grep("Bact",colnames(S16_readsB)),
                 names_to = "Taxon",
                 values_to = "reads")

	ST$Taxon <- sapply(ST$Taxon, ShortenNames6)

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
palette2 <- RColorBrewer::brewer.pal(n = 10, "Paired")
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



#if loop to save in correct place
if (median(rowSums(S16_readsB[,-1], na.rm = TRUE)) == 5000) {
ggplot2::ggsave(filename = file.path(folder,"/outputData/rarefied_read_taxaplot.png"), plot = taxaplot, height = 8, width = 15, units = "in")
} else{
ggplot2::ggsave(filename = file.path(folder,"/outputData/raw_read_taxaplot.png"), plot = taxaplot, height = 8, width = 15, units = "in")
}
messageColour(" Taxa Plot saved in Output folder \n\n", "message")

}
