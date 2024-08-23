#' Check quality of raw data, output graphs
#'
#' @param folder A string. Location of raw data files.
#'
#'

#quality_check
step1.01_quality_check <- function(folder) {
  messageColour("Function starting: 'quality_check'
  Checking quality of fastq data \n\n", "message")

  # Get unexported function to plot graphs
  plotCycleQuality <- utils::getFromNamespace(".plotCycleQuality", "ShortRead")

  #class(plotCycleQuality)#function.

  # Get file names to check
  fls <- dir(folder, pattern = "*fastq*", full.names = TRUE)

  # Generate quality checks report
  qaSummary <- ShortRead::qa(fls, type = "fastq")
  perCycle <- qaSummary[["perCycle"]]

  # Average quality scores
  allSamples <- perCycle
  allSamples$quality$lane[grepl("_R1_", allSamples$quality$lane,ignore.case = TRUE)]="Forward reads"
  allSamples$quality$lane[grepl("_R2_", allSamples$quality$lane,ignore.case = TRUE)]="Reverse reads"


  # Plot scores
  #library(ragg)
  path1 <- file.path(folder, "outputData", "qualityScores.png")
  per_plot<-plotCycleQuality(perCycle$quality)
  #class(per_plot)
  grDevices::png(path1, width = 1920, height = 1080); print(per_plot);

  grDevices::dev.off()

  path2 <- file.path(folder, "outputData", "qualityScores_averaged.png")
  avg_plot<-plotCycleQuality(allSamples$quality)
  grDevices::png(path2, width = 1920, height = 1080); print(avg_plot); grDevices::dev.off()

  messageColour("Details of checks saved in 'outputData' folder. \n\n", "message")
  messageColour("Function finished: 'quality_check' \n\n", "message")

}
