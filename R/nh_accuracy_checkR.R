#' @name nh_accuracy_checkR
#' @title NewHybrids accuracy checker
#'
#' @description \code{nh_accuracy_checkR} calcultes the accuracy with which NewHybrids assigns (simulated) individuals to their known genotype frequency category. That is, what proportion of all individuals that are known to belong to a given category were classified as belonging to that category by NH. The proportion correct will be returned as a dataframe with the name "NH.accuracy", and can also be printed to the screen if desired
#' @param NHResultsDir A file path to the NewHybrids (PofZ file) result to be checked.
#' @param print.results a logical query for whether the function output should be printed to the screen in addition to exported as an object. Default is TRUE
#' @param all.hyb a logical query if the proportion of all indivuals known to be hybrids were assigned to a hybrid category regardless if the category was the correct one. Default is FALSE
#' @importFrom stringr str_detect str_split
#' @export



nh_accuracy_checkR <- function(NHResultsDir, print.results = TRUE, all.hyb = FALSE){



  NHResultsDir_Split <- unlist(str_split(string = NHResultsDir, pattern = "/"))

  if(length(which(NHResultsDir_Split == "")) > 0){
    NHResultsDir_Split <- NHResultsDir_Split[-which(NHResultsDir_Split == "")]
  }

  # NHResultsDir_Split <- NHResultsDir_Split[-which(str_detect(string = NHResultsDir_Split, pattern = "PofZ"))]

  NHResultsDir_Dir <- paste(NHResultsDir_Split, collapse = "/")
  NHResultsDir_Dir <- paste0(NHResultsDir_Dir, "/")

  NHResults_FileList <- list.files(NHResultsDir_Dir)

  NHResults_FileList_indivs <-NHResults_FileList[which(str_detect(string = NHResults_FileList, pattern = "_individuals"))]
  NHResults_FileList_results <- NHResults_FileList[which(str_detect(string = NHResults_FileList, pattern = "PofZ"))]

  NHResults <- read.table(file = paste0(NHResultsDir, NHResults_FileList_results), header = TRUE)[,-2]

  colnames(NHResults) <- c("ind", "Pure1", "Pure2", "F1", "F2", "BC1", "BC2")

  NHindivs <- paste0(NHResultsDir, NHResults_FileList_indivs)

  Output <- n_class(x = NHindivs)

  Pure1 <- Output[1,2]
  Pure2 <- Output[2,2]
  F1 <- Output[3,2]
  F2 <- Output[4,2]
  BC1 <- Output[5,2]
  BC2 <- Output[6,2]

  Pure1.inds <- 1:Pure1
  Pure2.inds <- (Pure1 + 1):(Pure1 + Pure2)
  F1.inds <- (Pure1 + Pure2 + 1):(Pure1 + Pure2 + F1)
  F2.inds <- (Pure1 + Pure2 + F1 + 1):(Pure1 + Pure2 + F1 + F2)
  BC1.inds <- (Pure1 + Pure2 + F1 + F2 + 1):(Pure1 + Pure2 + F1 + F2 + BC1)
  BC2.inds <- (Pure1 + Pure2 + F1 + F2 + BC1 + 1):length(NHResults)



NHResults$pHyb <- rowSums(NHResults[4:7])




if(sum(NHResults$Pure1[Pure1.inds]) > sum(NHResults$Pure2[Pure1.inds])){

   ## calculate the number of indiviudals in the known populations whose probability of being in the correct population exceeds
        ## given by NH exceeds the given level of stringency / the total number of individuals in that population

  ## proportion correct at p = 0.5
  prop.corr.pure1.p.5 <- length(which(NHResults$Pure1[Pure1.inds] > 0.5))/Pure1
  prop.corr.pure2.p.5 <- length(which(NHResults$Pure2[Pure2.inds] > 0.5))/Pure2
  prop.corr.F1.p.5 <- length(which(NHResults$F1[F1.inds] > 0.5))/F1
  prop.corr.F2.p.5 <- length(which(NHResults$F2[F2.inds] > 0.5))/F2
  prop.corr.pure1BC.p.5 <- length(which(NHResults$BC1[BC1.inds] > 0.5))/BC1
  prop.corr.pure2BC.p.5 <- length(which(NHResults$BC2[BC2.inds] > 0.5))/BC2

  ## proportion correct at p=0.75

  prop.corr.pure1.p.75 <- length(which(NHResults$Pure1[Pure1.inds] > 0.75))/Pure1
  prop.corr.pure2.p.75 <- length(which(NHResults$Pure2[Pure2.inds] > 0.75))/Pure2
  prop.corr.F1.p.75 <- length(which(NHResults$F1[F1.inds] > 0.75))/F1
  prop.corr.F2.p.75 <- length(which(NHResults$F2[F2.inds] > 0.75))/F2
  prop.corr.pure1BC.p.75 <- length(which(NHResults$BC1[BC1.inds] > 0.75))/BC1
  prop.corr.pure2BC.p.75 <- length(which(NHResults$BC2[BC2.inds] > 0.75))/BC2

  ## proportion correct at p=0.9

  prop.corr.pure1.p.9 <- length(which(NHResults$Pure1[Pure1.inds] > 0.9))/Pure1
  prop.corr.pure2.p.9 <- length(which(NHResults$Pure2[Pure2.inds] > 0.9))/Pure2
  prop.corr.F1.p.9 <- length(which(NHResults$F1[F1.inds] > 0.9))/F1
  prop.corr.F2.p.9 <- length(which(NHResults$F2[F2.inds] > 0.9))/F2
  prop.corr.pure1BC.p.9 <- length(which(NHResults$BC1[BC1.inds] > 0.9))/BC1
  prop.corr.pure2BC.p.9 <- length(which(NHResults$BC2[BC1.inds] > 0.9))/BC2


  ## should the proportion of all indiviuals known to be hyybrids which have been identified as hybrids be calculated at each level of stringency
  if(all.hyb == TRUE){
    hyb.det.ALL.p.5 <- length(which(NHResults$pHyb[min(F1):nrow(NHResults)] > 0.5))/(F1 + F2 + BC1 + BC2)
    hyb.det.ALL.p.75 <- length(which(NHResults$pHyb[min(F1):nrow(NHResults)] > 0.75))/(F1 + F2 + BC1 + BC2)
    hyb.det.ALL.p.9 <- length(which(NHResults$pHyb[min(F1):nrow(NHResults)] > 0.9))/(F1 + F2 + BC1 + BC2)
  }

}

### if the First population given to NH has been denoted Pop2 by NH
if(sum(NHResults$Pure1[Pure1.inds]) > sum(NHResults$Pure2[Pure1.inds])){

  ## proportion correct at p=0.5

  prop.corr.pure1.p.5 <- length(which(NHResults$Pure1[Pure2.inds] > 0.5))/Pure2
  prop.corr.pure2.p.5 <- length(which(NHResults$Pure2[Pure1.inds] > 0.5))/Pure1
  prop.corr.F1.p.5 <- length(which(NHResults$F1[F1.inds] > 0.5))/F1
  prop.corr.F2.p.5 <- length(which(NHResults$F2[F2.inds] > 0.5))/F2
  prop.corr.pure1BC.p.5 <- length(which(NHResults$BC1[BC2.inds] > 0.5))/BC2
  prop.corr.pure2BC.p.5 <- length(which(NHResults$BC2[BC1.inds] > 0.5))/BC1

  ## proportion correct at p=0.75

  prop.corr.pure1.p.5 <- length(which(NHResults$Pure1[Pure2.inds] > 0.75))/Pure2
  prop.corr.pure2.p.5 <- length(which(NHResults$Pure2[Pure1.inds] > 0.75))/Pure1
  prop.corr.F1.p.5 <- length(which(NHResults$F1[F1.inds] > 0.75))/F1
  prop.corr.F2.p.5 <- length(which(NHResults$F2[F2.inds] > 0.75))/F2
  prop.corr.pure1BC.p.5 <- length(which(NHResults$BC1[BC2.inds] > 0.75))/BC2
  prop.corr.pure2BC.p.5 <- length(which(NHResults$BC2[BC1.inds] > 0.75))/BC1

   ## proportion correct at p=0.9

  prop.corr.pure1.p.5 <- length(which(NHResults$Pure1[Pure2.inds] > 0.9))/Pure2
  prop.corr.pure2.p.5 <- length(which(NHResults$Pure2[Pure1.inds] > 0.9))/Pure1
  prop.corr.F1.p.5 <- length(which(NHResults$F1[F1.inds] > 0.9))/F1
  prop.corr.F2.p.5 <- length(which(NHResults$F2[F2.inds] > 0.9))/F2
  prop.corr.pure1BC.p.5 <- length(which(NHResults$BC1[BC2.inds] > 0.9))/BC2
  prop.corr.pure2BC.p.5 <- length(which(NHResults$BC2[BC1.inds] > 0.9))/BC1

  if(all.hyb == TRUE){
    hyb.det.ALL.p.5 <- length(which(NHResults$pHyb[min(F1):nrow(NHResults)] > 0.5))/(F1 + F2 + BC1 + BC2)
    hyb.det.ALL.p.75 <- length(which(NHResults$pHyb[min(F1):nrow(NHResults)] > 0.75))/(F1 + F2 + BC1 + BC2)
    hyb.det.ALL.p.9 <- length(which(NHResults$pHyb[min(F1):nrow(NHResults)] > 0.9))/(F1 + F2 + BC1 + BC2)
  }

}

if(print.results == TRUE){
print(paste("Pure Population 1 Prop Correct 0.5", prop.corr.pure1.p.5))
print(paste("Pure Population 2 Prop Correct 0.5", prop.corr.pure2.p.5))
print(paste("F1 Prop Correct 0.5", prop.corr.F1.p.5))
print(paste("F2 Prop Correct 0.5", prop.corr.F2.p.5))
print(paste("Pure Population 1 BC Prop Correct 0.5", prop.corr.pure1BC.p.5))
print(paste("Pure Population 2 BC Prop Correct 0.5", prop.corr.pure2BC.p.5))
if(all.hyb == TRUE){
  print(paste("Hybrid Detected Regardless of Class 0.5", hyb.det.ALL.p.5))
}
print("")
print(paste("Pure Population 1 Prop Correct 0.75", prop.corr.pure1.p.75))
print(paste("Pure Population 2 Prop Correct 0.75", prop.corr.pure2.p.75))
print(paste("F1 Prop Correct 0.75", prop.corr.F1.p.75))
print(paste("F2 Prop Correct 0.75", prop.corr.F2.p.75))
print(paste("Pure Population 1 BC Prop Correct 0.75", prop.corr.pure1BC.p.75))
print(paste("Pure Populatin 2 BC Prop Correct 0.75", prop.corr.pure2BC.p.75))
if(all.hyb == TRUE){
  print(paste("Hybrid Detected Regardless of Class 0.75", hyb.det.ALL.p.75))
}
print("")
print(paste("Pure Population 1 Prop Correct 0.9", prop.corr.pure1.p.9))
print(paste("Pure Population 2 Prop Correct 0.9", prop.corr.pure2.p.9))
print(paste("F1 Prop Correct 0.9", prop.corr.F1.p.9))
print(paste("F2 Prop Correct 0.9", prop.corr.F2.p.9))
print(paste("Pure Population 1 BC Prop Correct 0.9", prop.corr.pure1BC.p.9))
print(paste("Pure Population 2 BC Prop Correct 0.9", prop.corr.pure2BC.p.9))
if(all.hyb == TRUE){
  print(paste("Hybrid Detected Regardless of Class 0.9", hyb.det.ALL.p.9))
}

} ## end of whether or not to print the results if test


## if the proportion of all hybs assigned correctly has been calculated, output

if(all.hyb == TRUE){
p.5 <- rbind(prop.corr.pure1.p.5, prop.corr.pure2.p.5, prop.corr.F1.p.5, prop.corr.F2.p.5, prop.corr.pure1BC.p.5,
  prop.corr.pure2BC.p.5,  hyb.det.ALL.p.5)
p.75 <- rbind(prop.corr.pure1.p.75, prop.corr.pure2.p.75, prop.corr.F1.p.75, prop.corr.F2.p.75, prop.corr.pure1BC.p.75,
  prop.corr.pure2BC.p.75,  hyb.det.ALL.p.75)
p.9 <- rbind(prop.corr.pure1.p.9, prop.corr.pure2.p.9, prop.corr.F1.p.9, prop.corr.F2.p.9, prop.corr.pure1BC.p.9,
  prop.corr.pure2BC.p.9,  hyb.det.ALL.p.9)

corr.names <- c("Pure1","Pure2", "F1", "F2", "Pure1BC", "Pure2BC", "All.Hyb")
prop.corr.output <- cbind(corr.names, p.5, p.75, p.9)
colnames(prop.corr.output)[2:4] <- c("p0.5", "p0.75", "p0.9")
}

## if the proportion of all individuals known to be hybs have not been calculated -- output
if(all.hyb == FALSE){
p.5 <- rbind(prop.corr.pure1.p.5, prop.corr.pure2.p.5, prop.corr.F1.p.5, prop.corr.F2.p.5, prop.corr.pure1BC.p.5,
  prop.corr.pure2BC.p.5)
p.75 <- rbind(prop.corr.pure1.p.75, prop.corr.pure2.p.75, prop.corr.F1.p.75, prop.corr.F2.p.75, prop.corr.pure1BC.p.75,
  prop.corr.pure2BC.p.75)
p.9 <- rbind(prop.corr.pure1.p.9, prop.corr.pure2.p.9, prop.corr.F1.p.9, prop.corr.F2.p.9, prop.corr.pure1BC.p.9,
  prop.corr.pure2BC.p.9)

corr.names <- c("Pure1","Pure2", "F1", "F2", "Pure1BC", "Pure2BC")
prop.corr.output <- cbind(corr.names, p.5, p.75, p.9)
colnames(prop.corr.output)[2:4] <- c("p0.5", "p0.75", "p0.9")
}


name.assign <- "NH.accuracy"
print(name.assign)
assign(x = name.assign, value = prop.corr.output, envir = globalenv())



}