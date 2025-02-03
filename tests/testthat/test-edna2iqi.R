test_that("edna2iqi function works", {
  
  path <- list.dirs(system.file("extdata/example_data/input", package = "eDNA2IQI"), full.names = TRUE)[1]
  test_data <- eDNA2IQI::edna2iqi(folder = path, auto_download = TRUE)
  preds <- test_data$myPreds
  
  testthat::expect_equal(round(preds$MeanIQI,4), round(c(0.3154),4))
  
})
