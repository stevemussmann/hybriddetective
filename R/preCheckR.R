#' @name preCheckR
#' @title Annotate Classes
#'
#' @description \code{preCheckR} detects cases where the NewHybrids MCMC chains fail to converge, and warns the user that the identified analyses should be repeated.
#' @param PreDir the directory in which the NewHybrids results (in individual folders as returned by parallelNH_XX) are located
#' @export

preCheckR <- function(PreDir){


  ##
  tbCheck <- list.files(PreDir)
CheckProgress <- txtProgressBar(min = 0, max = length(tbCheck), style = 3)
  possibleProbs_amChecking <- NULL
  for(i in tbCheck){
    j = which(tbCheck %in% i)

    setTxtProgressBar(CheckProgress, j)

    LiamChecking <- list.files(paste0(PreDir, i))
    pzCheckingFind <- LiamChecking[grep("PofZ", LiamChecking)]
    amChecking <- read.table(paste0(PreDir, i, "/", pzCheckingFind), header = T)

    #possibleProbs_amChecking <- c(possibleProbs_amChecking, pzCheckingFind)

    if(sum(amChecking[,6])>sum(sum(amChecking[,3]),sum(amChecking[,4]))){
      possibleProbs_amChecking <- c(possibleProbs_amChecking, pzCheckingFind)
      }

   # if(sum(amChecking$X0.250.0.250.0.250.0.250)>sum(sum(amChecking$X1.000.0.000.0.000.0.000),sum(amChecking$X0.000.0.000.0.000.1.000))){
    #  possibleProps_amChecking <- c(possibleProps_amChecking, pzCheckingFind)
     # }

  }

print("
  ")

  if(length(possibleProbs_amChecking) < 1){
    print("Looks good bud, giv'er")
  }
  if(length(possibleProbs_amChecking) >= 1){
    print(paste("Bud, you might want to have a closer look at", possibleProbs_amChecking, sep = " "))
  }



}