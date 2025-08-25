#' Example dataset: edna2iqi_paper
#'
#' A dataset bundled with the package for testing and demonstration.
#'
#' @format A list with two elements:
#' \describe{
#'   \item{data.frame}{A data frame with 762 observations and 867 variables. Each row corresponds to a sample, and each column corresponds to a taxonomic feature (e.g. ASVs, taxa).}
#'   \item{<other_element>}{A data frame of 762 observations and 13 variables. Each row corresponds to a sample, and each column corresponds to a metadata variable, including the actual IQI, farm name, etc.}
#' }
#'
#' @details
#' dataset used in internal testing and examples for the
#' `edna2iqiRFtrainer` package.
#' @source BactMetBar phase 1
"edna2iqi_paper"

#' Parameter options dataset
#'
#' A dataset listing all parameter options for denoising functions.
#'
#' @format ## `parameteroptions`
#' A data frame with 52 rows and 4 columns:
#' \describe{
#'   \item{Function}{Function name parameter belongs to}
#'   \item{Parameter}{Parameter name}
#'   \item{Value}{Parameter value}
#'   \item{Description}{Description of parameter, may include default value}
#' }
#'
#' @details
#' Dataset used in examples and testing to configure denoising steps in
#' the `edna2iqiRFtrainer` package.
"parameteroptions"
