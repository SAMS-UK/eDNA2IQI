#' Draw scatter plot
#'
#' Draw NMDS scatter plot of observed data compared to training data,
#' highlighting outliers.
#'
#' @param nmds_test output from internal nmds_bray() function.
#' @param folder active folder
#' @return NULL
#' @keywords internal

draw_nmds_plot <- function(nmds_test, folder) {
  # Plot applicability test output ---------------------------------------------
  # Create labels, titles, alpha etc variables to produce nice plot
  nmds_test$alpha <- 1
  # nmds_test$alpha[(nrow(AnnotatedASVs) + 1):nrow(nmds_test)] <- 0.3
  nmds_test$alpha[nmds_test$type == "Training"] <- 0.3
  nmds_test$label <- row.names(nmds_test)
  nmds_test$label[nmds_test$outlier == FALSE] <- ""
  #nmds_test$label[(nrow(AnnotatedASVs) + 1):nrow(nmds_test)] <- ""
  nmds_test$labela[nmds_test$type == "Training"] <- ""
#   nmds_test$outlier[(nrow(AnnotatedASVs) + 1):nrow(nmds_test)] <-
#     "Model training
# data"
  nmds_test$outlier[nmds_test$type == "Training"] <-
    "Model training
data"

  ggplot2::ggplot(nmds_test,
                  ggplot2::aes(
                    x = NMDS1,
                    y = NMDS2,
                    label = label,
                    colour = outlier
                  )) +
    ggplot2::geom_point(size = 2, alpha = nmds_test$alpha) +
    ggplot2::geom_text(vjust = -0.5,
                       hjust = 0.5,
                       size = 2) +
    ggplot2::labs(title = "NMDS of Taxa Data for each sample",
                  x = "NMDS1",
                  y = "NMDS2") +
    ggplot2::scale_color_manual(values = c("orange", "grey", "blue"))
  # Save the plot
  ggplot2::ggsave(file.path(folder, "/outputData/applicablity.png", fsep = ""))
  #  remove variables used for
  nmds_test <- dplyr::select(nmds_test, -label, -alpha)
  messageColour(" NMDS Plot saved in Output folder \n\n", "message")

}
