#' @name preCheckR
#' @title Annotate Classes
#'
#' @description \code{preCheckR} detects cases where the NewHybrids MCMC chains fail to converge, and warns the user that the identified analyses should be repeated. When NewHybrids fails to converge, it is generally visible in cumulative probability plots (see nh_plotR) and identifiable in the PofZ file as an excess of F2 posterior probability. This function checks that the proportion of known Pure Population 1 or Pure Population 2 individuals assigned at least a
#' @param PreDir the directory in which the NewHybrids results (in individual folders as returned by parallelNH_XX) are located.
#' @param propCutOff The proportion of individuals in either Pure Population 1 OR Pure Population 2 allowed to exceed the F2 assignment threshold (PofZ). The default is 0.5.
#' @param PofZCutOff The threshold posterior probability of assignment (PofZ) to F2 above which Pure Population 1 or Pure Population 2 individuals are flagged to indicate possible non-convergence. The default is 0.1.
#' @export

preCheckR <- function(PreDir, propCutOff = 0.5, PofZCutOff=0.1){

  if(propCutOff >= 1 | propCutOff < 0){
    stop("propCutOff must be a proportion between 0 and 1. Please refer to help file")
  }

  if(PofZCutOff >= 1 | PofZCutOff < 0){
    stop("PofZCutOff must be a proportion between 0 and 1. Please refer to help file")
  }

  writeLines("PrecheckR Progress: \r")
  ##
  tbCheck <- list.files(PreDir)
  CheckProgress <- txtProgressBar(min = 0, max = length(tbCheck), style = 3)
  possibleProbs_amChecking <- NULL
  for(i in tbCheck){
    j = which(tbCheck %in% i) ## changes i to a number so can be used in the progress bar

    setTxtProgressBar(CheckProgress, j)

    LiamChecking <- list.files(paste0(PreDir, i))
    pzCheckingFind <- LiamChecking[grep("PofZ", LiamChecking)]
    amChecking <- read.table(paste0(PreDir, i, "/", pzCheckingFind), header = T)

    indCheckingFind <- LiamChecking[grep("individuals", LiamChecking)]
    indChecking <- read.table(paste0(PreDir, i, "", indCheckingFild))
    Output <- n_class(x = indChecking)


    Prop_Pure1_F2_PofZ = length(which(amChecking[c(1:Output[1,2]), 6] > PofZCutOff))/Output[1,2]
    Prop_Pure2_F2_PofZ = length(which(amChecking[c((Output[1,2]+1):Output[2,2]), 6] > PofZCutOff))/Output[2,2]

    if(Prop_Pure1_F2_PofZ > propCutOff | Prop_Pure2_F2_PofZ > propCutOff){
      possibleProbs_amChecking <- c(possibleProbs_amChecking, pzCheckingFind)
        }


  }

  writeLines("\r")

  if(length(possibleProbs_amChecking) < 1){
    print("Looks good bud, giv'er")
  }
  if(length(possibleProbs_amChecking) >= 1){
    print(paste("Bud, you might want to have a closer look at", possibleProbs_amChecking, sep = " "))
      }



}