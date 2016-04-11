#' @name nh_multiplotR
#' @title NewHybrids Multi-Plot
#'
#' @description \code{nh_plotR} plots the cumulative probabilities of assignment for each individual for multiple NewHybrids results
#' @param NHResults the directory in which the NewHybrids results (in individual folders as returned by parallelNH_XX) are located.
#' @export
#' @import ggplot2
#' @importFrom reshape2 melt

NHResults <- "~/Desktop/DFO Aquaculture Interaction/Nova Scotia hybrid Analysis/Nova Scotia Analysis and R integration testing/NH.ResultsNS/"

nh_multiplotR <- function(NHResults){


  tbPlot <- list.files(NHResults)

  if(length(grep(x = tbPlot, pattern = "Figures and Data")) > 0){tbPlot = tbPlot[-which(tbPlot == "Figures and Data")]}
  if(length(grep(x = tbPlot, pattern = "NewHybrids Plots")) > 0){tbPlot = tbPlot[-which(tbPlot == "NewHybrids Plots")]}

  where.to <- paste0(NHResults, "NewHybrids Plots/")
  dir.create(where.to)

  writeLines("Multi-PlotR Progress: \r")
  ##

  PlotProgress <- txtProgressBar(min = 0, max = length(tbPlot), style = 3)

  for(i in tbPlot){
    j = which(tbPlot %in% i) ## changes i to a number so can be used in the progress bar

    setTxtProgressBar(PlotProgress, j)

    LiamPloting <- list.files(paste0(NHResults, i))
    pzPlotingFind <- LiamPloting[grep("PofZ", LiamPloting)]

    png(filename = paste0(where.to, i, ".png"), width = 2400, height = 2400, res = 300)

    print(hybriddetective::nh_plotR(paste0(NHResults, i, "/", pzPlotingFind)))
      dev.off()

  }


}