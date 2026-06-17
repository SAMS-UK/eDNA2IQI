#' Compare Non-parametric Multidimensional Scaling (NMDS) Centroid Distances
#'
#' The relationship between IQI and bacteria within the random forest model is
#' based on data collected to train the model. Sometimes if the data
#' observed is significantly different to the data in the training set, we
#' should have less confidence in the prediction.
#'
#' To test the applicability of the data observed for use with the model, we
#' compare the distance to the NMDS centroid of each sample in the observed data
#' to the 95th percentile distance to the centroid of samples in the training
#' set. Samples are flagged in the observed data which are more than 95th
#' percentile distance to centroid compared to the training set.
#'
#' @param data Observed data
#' @param training_data Data from training set.
#' @param metadata
#' @param trymax Maximum number of random starts in search of stable solution
#' @param probs Set percentile distance
#' @param data_output_only Return observed data only, if FALSE, training data is
#'   also returned
#'
#' @return Dataframe of observed samples.
#'
#' @examples
nmds_bray <- function(data,
                      training_data,
                      metadata,
                      trymax = 20,
                      probs = 0.95,
                      data_output_only = TRUE) {
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
  data <- data[complete.cases(data), ]
  row.names(data) <- row_names

  # Combine new data  points to training data.
  combined <- dplyr::bind_rows(data, training_data)

  # Run metaMDS
  nmds_result <- vegan::metaMDS(combined,
                         distance = "bray",
                         k = 2,
                         trymax = trymax,
                         try = trymax
  )

  site_scores <- as.data.frame(vegan::scores(nmds_result, display = "sites"))
  site_scores$Sample <- row.names(combined)

  # Get centroid for training data only and measure distance to that centroid
  training_scores <- site_scores[(nrow(site_scores) -
                                  nrow(training_data) + 1):nrow(site_scores), ]

  centroid <- colMeans(training_scores[, c("NMDS1", "NMDS2")])
  site_scores$dist_to_centroid <- sqrt((site_scores$NMDS1 - centroid[1])^2 +
                                         (site_scores$NMDS2 - centroid[2])^2)

  # Take the distance to centroid using the training data only quantile distance
  # to centroid
  quantile <- quantile(
    site_scores$dist_to_centroid[(nrow(site_scores) -
                                    nrow(training_data) + 1):nrow(site_scores)],
    probs = probs
  )
  site_scores$quantile <- quantile
  site_scores$diff <- site_scores$quantile - site_scores$dist_to_centroid

  # If distance to centroid is more than the quantile mark as outlier
  site_scores <- site_scores %>%
    dplyr::mutate(
      outlier = ifelse(diff < 0, TRUE, FALSE)
    )

  # Return only values for new data not training data
  if (data_output_only == TRUE) {
    site_scores <-
      site_scores[nrow(site_scores) - nrow(training_data):nrow(site_scores), ]
  } else {

  }
  return(site_scores)
}
