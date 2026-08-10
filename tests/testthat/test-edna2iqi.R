test_that("edna2iqi function works", {
  skip(
    "Needs threads creds to login (maybe add as secret to repo for testing purposes?"
  )

  path <- list.dirs(
    system.file("extdata/example_data/input", package = "eDNA2IQI"),
    full.names = TRUE
  )[1]

  test_data <- eDNA2IQI::edna2iqi(folder = path, auto_download = TRUE)
  preds <- test_data$myPreds
  testthat::expect_equal(round(preds$MeanIQI, 4)[2], round(c(0.3238), 4))
  # Test nmds test is removing outliers
  testthat::expect_equal(preds$MeanIQI[1], as.numeric(NA))
})
