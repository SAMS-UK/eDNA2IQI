#' Check quality of raw data, output graphs
#'
#' @param folder A string. Location of raw data files.
#' @return NULL output saved directly to folder. 
#'

step1.01_quality_check <- function(folder) {
  messageColour("Function starting: 'quality_check'
  Checking quality of fastq data \n\n", "message")

  # Get unexported function to plot graphs
  plotCycleQuality <- utils::getFromNamespace(".plotCycleQuality", "ShortRead")

  # Get file names to check
  fls <- dir(folder, pattern = "*fastq*", full.names = TRUE)

  # Generate quality checks report
  # Warnings supressed due to errors from node serialisation
  qaSummary <- suppressWarnings(ShortRead::qa(fls, type = "fastq"))
  perCycle <- suppressWarnings(qaSummary[["perCycle"]])

  # Average quality scores
  allSamples <- perCycle
  allSamples$quality$lane[grepl("_R1_", allSamples$quality$lane,ignore.case = TRUE)]="Forward reads"
  allSamples$quality$lane[grepl("_R2_", allSamples$quality$lane,ignore.case = TRUE)]="Reverse reads"

  # Plot scores
  path1 <- file.path(folder, "outputData", "qualityScores.png")
  per_plot<-plotCycleQuality(perCycle$quality)
  grDevices::png(path1, width = 1920, height = 1080); print(per_plot);

  grDevices::dev.off()

  path2 <- file.path(folder, "outputData", "qualityScores_averaged.png")
  avg_plot<-plotCycleQuality(allSamples$quality)
  grDevices::png(path2, width = 1920, height = 1080); print(avg_plot); grDevices::dev.off()

  messageColour("Details of checks saved in 'outputData' folder. \n\n", "message")
  messageColour("Function finished: 'quality_check' \n\n", "message")

}
