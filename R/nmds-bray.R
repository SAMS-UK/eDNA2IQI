#' Compare Non-parametric Multidimensional Scaling (NMDS) Centroid Distances
#'
#' The relationship between IQI and bacteria within the random forest model is
#' based on data collected to train the model. Sometimes if the data
#' observed is significantly different to the data in the training set, we
#' should have less confidence in the prediction.
#'
#' To test the applicability of the data observed for use with the model, we
#' compare the distance to NMDS coordinates of each sample in the observed data
#' to the 95th percentile distance to the centroid of samples in the training
#' set. Samples are flagged in the observed data which are more than 95th
#' percentile distance to centroid compared to the training set.
#'
#' Samples may not be applicable for the model due to environmental, sampling,
#' analysis or data processing anomalies. The cause should be investigated and
#' remedial steps taken where required.
#'
#' @param data Observed data
#' @param training_data Data from training set.
#' @param trymax Maximum number of random starts in search of stable solution
#' @param probs Set percentile distance
#' @param data_output_only Return observed data only, if FALSE, training data is
#'   also returned
#'
#' @return Dataframe of observed samples.
#' @keywords internal

nmds_bray <- function(
  data,
  training_data,
  trymax = 20,
  probs = 0.95,
  data_output_only = TRUE
) {
  message(
    "Internal Function Starting: 'nmds_bray()'
Testing the applicability of the data against the reference dataset used to train the model"
  )
  # Set seed for reproducibility
  set.seed(56)
  # Format 'new' data to test against the training data. Row names get removed
  # by next step so save to add back later.
  row_names <- row.names(data)
  data <- data[, names(data) %in% names(training_data)]
  # Check if taxa used for RF predictions are present in any sample
  missing_names <- setdiff(names(training_data), names(data))
  # Add columns with zero values for missing taxa
  if (!is.null(missing_names)) {
    data[missing_names] <- 0
  }
  data[is.na(data)] <- 0
  # Only include RF model taxa with reads
  data <- data[stats::complete.cases(data), ]
  row.names(data) <- row_names

  # Combine new data  points to training data.
  combined <- dplyr::bind_rows(data, training_data)

  # Run metaMDS
  nmds_result <- suppressMessages(vegan::metaMDS(
    combined,
    distance = "bray",
    k = 2,
    trymax = trymax,
    try = trymax,
    trace = FALSE
  ))

  site_scores <- as.data.frame(vegan::scores(nmds_result, display = "sites"))
  site_scores$Sample <- row.names(combined)

  # Get centroid for training data only and measure distance to that centroid
  training_scores <- site_scores[
    (nrow(site_scores) -
      nrow(training_data) +
      1):nrow(site_scores),
  ]

  centroid <- colMeans(training_scores[, c("NMDS1", "NMDS2")])
  site_scores$dist_to_centroid <- sqrt(
    (site_scores$NMDS1 - centroid[1])^2 +
      (site_scores$NMDS2 - centroid[2])^2
  )

  # Take the distance to centroid using the training data only quantile distance
  # to centroid
  quantile <- stats::quantile(
    site_scores$dist_to_centroid[
      (nrow(site_scores) -
        nrow(training_data) +
        1):nrow(site_scores)
    ],
    probs = probs
  )
  site_scores$quantile <- quantile
  site_scores$diff <- site_scores$quantile - site_scores$dist_to_centroid

  # If distance to centroid is more than the quantile mark as outlier
  site_scores <- site_scores %>%
    dplyr::mutate(outlier = ifelse(diff < 0, TRUE, FALSE))

  # Add Number of Taxa column (ntaxa)
  combined[combined > 0] <- 1
  combined$ntaxa <- rowSums(combined)
  site_scores$n_taxa <- combined$ntaxa
  # Add variable to show if training data or observed data
  site_scores$type <- ""
  site_scores$type[(nrow(data) + 1):nrow(site_scores)] <- "Training"
  site_scores$type[1:nrow(data)] <- "New data"

  # Return only values for new data, not training data
  if (data_output_only == TRUE) {
    site_scores <-
      site_scores[nrow(site_scores) - nrow(training_data):nrow(site_scores), ]
  } else {}
  return(site_scores)
}
