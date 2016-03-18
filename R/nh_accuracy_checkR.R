#' @name nh_accuracy_checkR
#' @title NewHybrids accuracy checker
#'
#' @description \code{nh_accuracy_checkR} calcultes the accuracy with which NewHybrids assigns (simulated) individuals to their known genotype frequency category. That is, what proportion of all individuals that are known to belong to a given category were classified as belonging to that category by NH. The proportion correct will be returned as a dataframe with the name "NH.accuracy", and can also be printed to the screen if desired
#' @param NHResults a dataframe that is the PofZ object
#' @param print.results a logical query for whether the results should be printed to the screen in addition to exported as an object. Default is TRUE
#' @param all.hyb a logical query if the proportion of all indivuals known to be hybrids were assigned to a hybrid category regardless if the category was the correct one. Default is FALSE
#' @export


NHresults <- read.table("~/Desktop/DFO Aquaculture Interaction/Nova Scotia hybrid Analysis/Nova Scotia Simulated Proportional Sampling/NStop240ForSimulation_S1R1_NH.txt_Results/NStop240ForSimulation_S1R1_NH.txt_PofZ.txt", header = TRUE)[,-2]

if(50% of P1 > 0.1 == TRUE | 50% of P2 > 0.1 == TRUE){FLAG}

nh_accuracy_checkR <- function(NHResults, print.results = TRUE, all.hyb = FALSE){


  NHResults <- read.table(NHResults, header = TRUE)[,-2]

  colnames(NHResults) <- c()


num.sim <- nrow(NHResults)/(length(NHResults)-1)
NHResults$pHyb <- rowSums(NHResults[4:7])


if(sum(NHResults$Pure1[1:200]) > sum(NHResults$Pure2[1:200])){

   ## calculate the number of indiviudals in the known populations whose probability of being in the correct population exceeds
        ## given by NH exceeds the given level of stringency / the total number of individuals in that population

  ## proportion correct at p = 0.5
  prop.corr.farm.p.5 <- length(which(NHResults$Pure1[1:200] > 0.5))/num.sim
  prop.corr.wild.p.5 <- length(which(NHResults$Pure2[201:400] > 0.5))/num.sim
  prop.corr.F1.p.5 <- length(which(NHResults$F1[401:600] > 0.5))/num.sim
  prop.corr.F2.p.5 <- length(which(NHResults$F2[601:800] > 0.5))/num.sim
  prop.corr.farmBC.p.5 <- length(which(NHResults$BC1[801:1000] > 0.5))/num.sim
  prop.corr.wildBC.p.5 <- length(which(NHResults$BC2[1001:1200] > 0.5))/num.sim

  ## proportion correct at p=0.75

  prop.corr.farm.p.75 <- length(which(NHResults$Pure1[1:200] > 0.75))/num.sim
  prop.corr.wild.p.75 <- length(which(NHResults$Pure2[201:400] > 0.75))/num.sim
  prop.corr.F1.p.75 <- length(which(NHResults$F1[401:600] > 0.75))/num.sim
  prop.corr.F2.p.75 <- length(which(NHResults$F2[601:800] > 0.75))/num.sim
  prop.corr.farmBC.p.75 <- length(which(NHResults$BC1[801:1000] > 0.75))/num.sim
  prop.corr.wildBC.p.75 <- length(which(NHResults$BC2[1001:1200] > 0.75))/num.sim

  ## proportion correct at p=0.9

  prop.corr.farm.p.9 <- length(which(NHResults$Pure1[1:200] > 0.9))/num.sim
  prop.corr.wild.p.9 <- length(which(NHResults$Pure2[201:400] > 0.9))/num.sim
  prop.corr.F1.p.9 <- length(which(NHResults$F1[401:600] > 0.9))/num.sim
  prop.corr.F2.p.9 <- length(which(NHResults$F2[601:800] > 0.9))/num.sim
  prop.corr.farmBC.p.9 <- length(which(NHResults$BC1[801:1000] > 0.9))/num.sim
  prop.corr.wildBC.p.9 <- length(which(NHResults$BC2[1001:1200] > 0.9))/num.sim


  ## should the proportion of all indiviuals known to be hyybrids which have been identified as hybrids be calculated at each level of stringency
  if(all.hyb == TRUE){
    hyb.det.ALL.p.5 <- length(which(NHResults$pHyb[401:1200] > 0.5))/(num.sim*4)
    hyb.det.ALL.p.75 <- length(which(NHResults$pHyb[401:1200] > 0.75))/(num.sim*4)
    hyb.det.ALL.p.9 <- length(which(NHResults$pHyb[401:1200] > 0.9))/(num.sim*4)
  }

}

### if the First population given to NH has been denoted Pop2 by NH
if(sum(NHResults$Pure1[1:200]) < sum(NHResults$Pure2[1:200])){

  ## proportion correct at p=0.5

   prop.corr.farm.p.5 <- length(which(NHResults$Pure1[201:400] > 0.5))/num.sim
  prop.corr.wild.p.5 <- length(which(NHResults$Pure2[1:200] > 0.5))/num.sim
  prop.corr.F1.p.5 <- length(which(NHResults$F1[401:600] > 0.5))/num.sim
  prop.corr.F2.p.5 <- length(which(NHResults$F2[601:800] > 0.5))/num.sim
  prop.corr.farmBC.p.5 <- length(which(NHResults$BC1[1001:1200] > 0.5))/num.sim
  prop.corr.wildBC.p.5 <- length(which(NHResults$BC2[801:1000] > 0.5))/num.sim

  ## proportion correct at p=0.75

  prop.corr.farm.p.75 <- length(which(NHResults$Pure1[201:400] > 0.75))/num.sim
  prop.corr.wild.p.75 <- length(which(NHResults$Pure2[1:200] > 0.75))/num.sim
  prop.corr.F1.p.75 <- length(which(NHResults$F1[401:600] > 0.75))/num.sim
  prop.corr.F2.p.75 <- length(which(NHResults$F2[601:800] > 0.75))/num.sim
  prop.corr.farmBC.p.75 <- length(which(NHResults$BC1[1001:1200] > 0.75))/num.sim
  prop.corr.wildBC.p.75 <- length(which(NHResults$BC2[801:1000] > 0.75))/num.sim

   ## proportion correct at p=0.9

  prop.corr.farm.p.9 <- length(which(NHResults$Pure1[201:400] > 0.9))/num.sim
  prop.corr.wild.p.9 <- length(which(NHResults$Pure2[1:200] > 0.9))/num.sim
  prop.corr.F1.p.9 <- length(which(NHResults$F1[401:600] > 0.9))/num.sim
  prop.corr.F2.p.9 <- length(which(NHResults$F2[601:800] > 0.9))/num.sim
  prop.corr.farmBC.p.9 <- length(which(NHResults$BC1[1001:1200] > 0.9))/num.sim
  prop.corr.wildBC.p.9 <- length(which(NHResults$BC2[801:1000] > 0.9))/num.sim

  if(all.hyb == TRUE){
    hyb.det.ALL.p.5 <- length(which(NHResults$pHyb[401:1200] > 0.5))/(num.sim*4)
    hyb.det.ALL.p.75 <- length(which(NHResults$pHyb[401:1200] > 0.75))/(num.sim*4)
    hyb.det.ALL.p.9 <- length(which(NHResults$pHyb[401:1200] > 0.9))/(num.sim*4)
  }

}

if(print.results == TRUE){
print(paste("Farm Prop Correct 0.5", prop.corr.farm.p.5))
print(paste("Wild Prop Correct 0.5", prop.corr.wild.p.5))
print(paste("F1 Prop Correct 0.5", prop.corr.F1.p.5))
print(paste("F2 Prop Correct 0.5", prop.corr.F2.p.5))
print(paste("Farm BC Prop Correct 0.5", prop.corr.farmBC.p.5))
print(paste("Wild BC Prop Correct 0.5", prop.corr.wildBC.p.5))
if(all.hyb == TRUE){
  print(paste("Hybrid Detected Regardless of Class 0.5", hyb.det.ALL.p.5))
}
print("")
print(paste("Farm Prop Correct 0.75", prop.corr.farm.p.75))
print(paste("Wild Prop Correct 0.75", prop.corr.wild.p.75))
print(paste("F1 Prop Correct 0.75", prop.corr.F1.p.75))
print(paste("F2 Prop Correct 0.75", prop.corr.F2.p.75))
print(paste("Farm BC Prop Correct 0.75", prop.corr.farmBC.p.75))
print(paste("Wild BC Prop Correct 0.75", prop.corr.wildBC.p.75))
if(all.hyb == TRUE){
  print(paste("Hybrid Detected Regardless of Class 0.75", hyb.det.ALL.p.75))
}
print("")
print(paste("Farm Prop Correct 0.9", prop.corr.farm.p.9))
print(paste("Wild Prop Correct 0.9", prop.corr.wild.p.9))
print(paste("F1 Prop Correct 0.9", prop.corr.F1.p.9))
print(paste("F2 Prop Correct 0.9", prop.corr.F2.p.9))
print(paste("Farm BC Prop Correct 0.9", prop.corr.farmBC.p.9))
print(paste("Wild BC Prop Correct 0.9", prop.corr.wildBC.p.9))
if(all.hyb == TRUE){
  print(paste("Hybrid Detected Regardless of Class 0.9", hyb.det.ALL.p.9))
}

} ## end of whether or not to print the results if test


## if the proportion of all hybs assigned correctly has been calculated, output

if(all.hyb == TRUE){
p.5 <- rbind(prop.corr.farm.p.5, prop.corr.wild.p.5, prop.corr.F1.p.5, prop.corr.F2.p.5, prop.corr.farmBC.p.5,
  prop.corr.wildBC.p.5,  hyb.det.ALL.p.5)
p.75 <- rbind(prop.corr.farm.p.75, prop.corr.wild.p.75, prop.corr.F1.p.75, prop.corr.F2.p.75, prop.corr.farmBC.p.75,
  prop.corr.wildBC.p.75,  hyb.det.ALL.p.75)
p.9 <- rbind(prop.corr.farm.p.9, prop.corr.wild.p.9, prop.corr.F1.p.9, prop.corr.F2.p.9, prop.corr.farmBC.p.9,
  prop.corr.wildBC.p.9,  hyb.det.ALL.p.9)

corr.names <- c("Farm","Wild", "F1", "F2", "FarmBC", "WildBC", "All.Hyb")
prop.corr.output <- cbind(corr.names, p.5, p.75, p.9)
colnames(prop.corr.output)[2:4] <- c("p0.5", "p0.75", "p0.9")
}

## if the proportion of all individuals known to be hybs have not been calculated -- output
if(all.hyb == FALSE){
p.5 <- rbind(prop.corr.farm.p.5, prop.corr.wild.p.5, prop.corr.F1.p.5, prop.corr.F2.p.5, prop.corr.farmBC.p.5,
  prop.corr.wildBC.p.5)
p.75 <- rbind(prop.corr.farm.p.75, prop.corr.wild.p.75, prop.corr.F1.p.75, prop.corr.F2.p.75, prop.corr.farmBC.p.75,
  prop.corr.wildBC.p.75)
p.9 <- rbind(prop.corr.farm.p.9, prop.corr.wild.p.9, prop.corr.F1.p.9, prop.corr.F2.p.9, prop.corr.farmBC.p.9,
  prop.corr.wildBC.p.9)

corr.names <- c("Farm","Wild", "F1", "F2", "FarmBC", "WildBC")
prop.corr.output <- cbind(corr.names, p.5, p.75, p.9)
colnames(prop.corr.output)[2:4] <- c("p0.5", "p0.75", "p0.9")
}


name.assign <- "NH.accuracy"
print(name.assign)
assign(x = name.assign, value = prop.corr.output, envir = globalenv())



}