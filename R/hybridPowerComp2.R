#' @name hybridPowerComp2
#' @title Assignment power comparison among different SNP subsets using NewHybrids simulated datasets2.
#' @description Evaluates the accuracy with which NewHybrids assigns individuals of known hybrid class to the correct hybrid class in simulated datasets at varying levels of stringency (PofZ). The code will write graphical and numerical results to the directory provided by the user.
#' @param dir File path to the directory in which the NewHybrids results (in individual folders as returned by parallelNH_XX) are located.
#' @param filetag An optional character vector to be applied as the name of the outputs.
#' @param Thresholds A vector of thresholds which will be added to the plots showing the assignment success for different levels of probability of a given class estimated by NewHybrids. Default is (NULL) so if nothing is specified it will not add this to the output plots (success ~ threshold by class).
#' @param samplesize The number of individuals per NewHybrids class. By (default: NULL) this data will be extracted from the "*individuals.txt" output from parallelnewhybrids if present in the same folder as the PofZ file. This can also explicitly defined as a vector (6 values corresponding to # in P1,P2,F1,F2,BC1,BC2) or a path to a *_Individuals.txt.
#' @param CT The threshold posterior probability of assignment (PofZ) to F2 above which Pure Population 1 or Pure Population 2 individuals are flagged to indicate possible non-convergence. The default is 0.1.
#' @param CTI The proportion of individuals in either Pure Population 1 OR Pure Population 2 allowed to exceed the F2 assignment threshold (PofZCutOff). The default is 0.5.
#' @rdname hybridpowercomp2
#' @import ggplot2
#' @import magrittr
#' @importFrom dplyr filter summarise ungroup group_by
#' @importFrom grid arrow unit
#' @importFrom stringr str_extract
#' @importFrom reshape2 melt
#' @importFrom  scales alpha
#' @export

#
 #  library(ggplot2)
 #  library(magrittr)
 #  library(dplyr)
 #  library(stringr)
 #  library(reshape2)
 #  library(grid)
 #  library(scales)
 #  library(hybriddetective)
 #
 #
 #
 # dir = "~/Desktop/DFO Aquaculture Interaction/Nova Scotia hybrid Analysis/Nova Scotia Analysis and R integration testing/NSTop48-1000-WithZed//"
 #  filetag=""
 #  Thresholds=c(0.5,0.6,0.7,0.8,0.9)
 #  addThresh=FALSE
 #  samplesize=NULL
 #  CT=0.1
 #  CTI=0.5



hybridPowerComp2 <-function(dir,filetag="",Thresholds=c(0.5,0.6,0.7,0.8,0.9),addThresh=FALSE,samplesize=NULL,CT=0.1,CTI=0.5){




  #set directory for which holds the New Hybrids output folders
  filedir <- dir
  lfiles <- setdiff(list.files(dir),c("Figures and Data", "NewHybrids Plots")) #ignores Figures folder in case this is run more than once and in case plots made
    if(length(which(list.files(dir)=="Figures and Data"))==0){dir.create(paste0(dir,"Figures and Data"))} # if there isn't a 'Figures and Data' folder for output create one
    if(length(which(list.files(paste0(dir,"Figures and Data"))=="pdf"))==0){dir.create(paste0(dir,"Figures and Data/pdf"))} #create a folder for pdfs
    if(length(which(list.files(paste0(dir,"Figures and Data"))=="jpg"))==0){dir.create(paste0(dir,"Figures and Data/jpg"))} #create a folder for jpgs
    if(length(which(list.files(paste0(dir,"Figures and Data"))=="data"))==0){dir.create(paste0(dir,"Figures and Data/data"))} #create a folder for data

  #Convergence checker
  arethereproblems = "no"

    # Collate the output from New Hybrids together ('p of z' files)
        output <- NULL
        for (i in lfiles)
        {
              tempfiles <- list.files(paste0(filedir,i))
              pzfile <- tempfiles[grep("PofZ",tempfiles)]
              tempfile <- read.table(paste0(filedir,i,"/",pzfile),head=T)

              LociandAlleles <- tempfiles[grep("LociAndAlleles", tempfiles)]
              LandAfile <- readChar(paste0(filedir, i, "/", LociandAlleles), file.info(paste0(filedir, i, "/", LociandAlleles))$size)
              numLociExt <- str_extract(string = LandAfile, pattern = paste0("from ", "[:digit:]{1,5}", " loci"))
              numLociWorking <- gsub(x = numLociExt, pattern = "from ", replacement = "")
              numLociWorking <- as.numeric(gsub(x = numLociWorking, pattern = " loci", replacement = ""))

              #identify the simulation and repeat info
              S_ident <- gsub("_","",str_extract(pzfile,paste0("_S","[:digit:]{1}","R","[:digit:]{1}","_")))
              tempfile$sim <- substring(S_ident,1,2)
              tempfile$rep <- substring(S_ident,3,4)
              tempfile$nLoci <- numLociWorking

              tempfile <- tempfile[,-grep("IndivName",colnames(tempfile))] #delete IndivName

              #rename the columns
              colnames(tempfile) <- c("Indv","Pure1","Pure2","F1","F2","BC1","BC2","sim","rep","nLoci")
              tempfile=tempfile[,c("Indv","sim","rep","nLoci","Pure1","Pure2","F1","F2","BC1","BC2")]# reorder

                #Get the samplesize for a given class
                IndividualsPath <- tempfiles[grep("individuals.txt", tempfiles)]
                  if(length(samplesize)==1 & is.numeric(samplesize)){samplesize <- rep(samplesize,6)}
                  if(length(samplesize)==1 & !is.numeric(samplesize)){samplesize <- as.vector(n_class(samplesize)[,2])}
                  if(is.null(samplesize)){samplesize <- as.vector(n_class(paste0(filedir, i, "/", IndividualsPath)))[,2]}

                #common order
                if(sum(tempfile[1:samplesize[1],"Pure1"],na.rm=T)<sum(tempfile[1:samplesize[1],"Pure2"],na.rm=T)){
                  pure1 <- tempfile$Pure2;pure2 <- tempfile$Pure1
                  bc1 <- tempfile$BC2;bc2 <- tempfile$BC1

                  tempfile$Pure1 <- pure1;tempfile$Pure2 <- pure2
                  tempfile$BC1 <- bc1;tempfile$BC2 <- bc2
                }

              #Filter for convervence issues. Based on the 'convergence filter (CT) and % of indviduals permited to fail (CTI)
              #Here we look at the "pure 1 and 2 populations for
                if(length(which(tempfile[1:samplesize[1],"F2"]>CT))/length(1:samplesize[1])>CTI &
                        length(which(tempfile[(samplesize[1]+1):samplesize[2],"F2"]>CT))/length((samplesize[1]+1):samplesize[2])>CTI){

                      tempfile[,5:length(tempfile)]=NA #replace data with NAs
                      print(paste("Possible non-convergence detected in", pzfile))
                      arethereproblems = "Yes"}

              output <- rbind(output,tempfile)

          }#end of for loop

        #If convergence issues were flagged then the process stops here.
        if(arethereproblems == "Yes")
          {
            stop("Please remove, or re-run those results for which non-convergence was detected", call. = F)
          }


  temp = output
  tempinds <- n_class(paste0(filedir, lfiles[1], "/", IndividualsPath))
  out.inds.class = NULL

  ## get the correct numbers of individuals for known categories
  for(j in 1:length(samplesize)){
    make.class <- rep(as.character(tempinds[j,1]), times = tempinds[j,2])
    out.inds.class <- c(out.inds.class, make.class)
    } # END J Loop

  temp$known <- out.inds.class
  colnames(temp)[c(5,6)] <- c("P1", "P2") ## rename these two columns to match column names between all data frames made
  temp$max.class <- NA ## add a column that will be the class to which each indiviudal is assigned

  ## calculate what the most probable (highest PofZ) genotype frequency category is for each indivivudal
  for(k in 1:nrow(temp)){
    temp$max.class[k]=names(which.max(temp[k, 5:10]))
  }



  ########################
  ## Calculate ACCURACY ##
  ########################

    writeLines("Calculating Accuracy
        ")

  ## ACCURACY =  number assigned correctly / total number assigned  -> for each category
  ### Check for every value of critical PofZ from 0.5 to 0.99

  ### set up progress bar
  stat1_check_progress <- txtProgressBar(min = 0.5, max = 0.99, style = 3)

  temp2 <- temp
  final.stats.1 <- NULL
  for(l in 50:99/100){
    # i = 0.75
    setTxtProgressBar(stat1_check_progress, l) ## call progress bar
    temp2$isgood <- "YES" ## set all isgood to "YES" - default is correct <- will switch to "NO" if PofZ value < critical value
    temp2$domatch <- "NO" ## set all domatch to "NO" -  default is incorrect <- will swith to "YES" if they match


    ### Check if the assignment matches the known class, and if the PofZ value for the assignment is greater than threshold
    for(m in 1:nrow(temp2)){

      if(as.character(temp2[m, "known"]) == as.character(temp2[m, "max.class"])){ ## evaluate match
        temp2[m, "domatch"] = "YES"
          }

      if(temp2[m, 5:10][which.max(temp2[m, 5:10])] < l){ ## evaluate PofZ threshold
        temp2[m, "isgood"] = "NO"
          }

        } ### END m loop - checks if call and known class match, and if the assignment PofZ is greater than the threshold

      ## create logical to retain only instances where assignment matches known category, and the PofZ is greater than the critical
      only_good_results <- temp2$isgood == "YES" & temp2$domatch == "YES"

      temp3 <- temp2[only_good_results, ]
      temp3$known <- factor(x = temp3$known, levels = c("P1", "P2", "F1", "F2", "BC1", "BC2")) ## change factor levels <- not sure that this is needed, because table() seems to screw it up anyways

      ## create logical that looks if the assignment is retained (above PofZ), no matter if it is correct (matches known) or not
      only_matched_results <- temp2$isgood   == "YES"
      temp4 <- temp2[only_matched_results, ]

      ## NEED to look at the number of individuals assigned to each category by: simulation, replicate, and number of loci
      ### WANT to average WITHIN simulation
      ## table() will create counts of individuals
      ## CREATE tables that are simulations by replicates, within each genotype category, within each number of loci
      ## table[simulation, replicate, category, loci]

      bigtable <- table(temp3$sim, temp3$rep, temp3$max.class, temp3$nLoci) ### CORRECT TABLE - assigned correctly at critical PofZ (numerator)
      bigtable2 <- table(temp4$sim, temp4$rep, temp4$max.class, temp4$nLoci) ### Assigned TABLE - assigned to group at critical PofZ (denominator)

      loci.groups <- unlist(dimnames(bigtable)[length(dimnames(bigtable))]) ## what are the numbers of loci - allows to be variable, and will be populated from the data provided
      known.groups <- unlist(dimnames(bigtable)[(length(dimnames(bigtable))-1)]) ## what are the genotype frequency classes evaluated - allows to be variable, and will be populated from data provided
      num.sims <- length(unique(temp3$sim)) ### how man simulations were conducted - allows to be variable, and will be populated from the data provided


      ### check the accuracy for each level of numbrr of loci - checks by loci number, then classes
      # o = 1
      by.loci.stats = NULL
      for(o in 1:length(loci.groups)){

        temp.by.loci <- bigtable[, , , loci.groups[o]] ## subset the CORRECT table for the kth number of loci - numerator
        temp.by.loci2 <- bigtable2[, , , loci.groups[o]] ## subset the Assigned table for the kth number of loci - denominator

        ### within the oth loci group, calculate values for each genotype frequency group
        out.stats = NULL

          ## start p loop
          for(p in 1:length(known.groups)){

            temp.by.group <- temp.by.loci[, , known.groups[p]] ## subset table by the hth genotype frequency category
            temp.by.group2 <- temp.by.loci2[, , known.groups[p]] ### subset table by the hth genotype frequency category
            ## calculate the means for each simulation - know that the data will now be a rectangular table with simulations as rows, and replicates as columns. <- so, calculate means of rows using apply, then divide the number correct by the number assigned in total to get the accuracy
            means <- apply(X = temp.by.group, MARGIN = 1, FUN = mean)/apply(X = temp.by.group2, MARGIN = 1, FUN = mean)

            ### create a dataframe made up of the genotype frequency class, the simulations, and the means for those simulations for that genotype frequency class
            hold.stats <- data.frame(rep(known.groups[p], times = length(row.names(temp.by.group))), row.names(temp.by.group), means)
            colnames(hold.stats)[1:2] <- c("known", "sim")
            out.stats <- rbind(out.stats, hold.stats) ## bind all the results for the genotype frequency classes together
              } # END p loop

        ## add the number of loci that the results just calculated represent
        out.stats$nloci <- loci.groups[o]
        ## bind the results just calculated for the kth number of loci, to all resutls calcualted so far - remember though, this is just within one value of critical PofZ, will be NULLed above to replicate calculation for each PofZ
        by.loci.stats = rbind(by.loci.stats, out.stats)
        ### assign the results for the ith critical PofZ to the k groups of loci to a new variable because by.loci.stats will be NULLed each k loop
        final.stats.hold <- by.loci.stats

        }# END o LOOP

      ## assign the lth PofZ critical value to the results output from teh o loop,
      final.stats.hold$PofZ <- l
      ## bind all i critical PofZ results togethwer
      final.stats.1 <- rbind(final.stats.1, final.stats.hold)

      } ## END ACCURACY LOOP


      ## calculate and plot two different things
      ## first the accuracy at critical PofZ > 0.5, 0.75 and 0.9
      ### then accuracy at all critical PofZ 0.5 < PofZ < 0.99

        ## rename while creating fucntion to not have to re-run loop over and over if make mistake
        final.stats.2 <- final.stats.1
        ### make number of loci a factor with levels increasing from smallest to largest value
        final.stats.2$nloci <- factor(x = final.stats.2$nloci, levels = ordered(unique(as.numeric(final.stats.2$nloci))))
        ## change the level values of the genotype frequency categories
        final.stats.2$known <- factor(x = final.stats.2$known, levels = c("P1", "P2", "F1", "F2", "BC1", "BC2"))

        ## create a logical value that will choose only PofZ = 0.5, 0.75, and 0.9
        pofzeds <- final.stats.2$PofZ == 0.5 | final.stats.2$PofZ == 0.75 | final.stats.2$PofZ == 0.9
        ## get these values
        final.stats.pofzeds <- final.stats.2[pofzeds, ]

        ## plot boxplot of accuracy at PofZ =0.5, 0.75 and 0.9
       accuracy_boxplot <- ggplot(final.stats.pofzeds, aes(x = known, y = means, fill = known)) + geom_boxplot() + facet_grid(PofZ~nloci) + labs(x = "Genotype Frequency Class", y = "Proportion of Assignments Correct") +
           scale_fill_brewer(palette = "Dark2") +
           theme(panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"), legend.position = "none", strip.background = element_rect(, colour = "black", fill = "white"), strip.text.x = element_text(colour = "black"), strip.text.y = element_text(colour = "black"))

           #     #Save plot
    if(filetag!=""){ggsave(paste0(dir,"Figures and Data/pdf/", filetag, "_AccuracyBoxPlot.pdf"), accuracy_boxplot, height = 10, width = 10)}else
      {ggsave(paste0(dir, "Figures and Data/pdf/AccuracyBoxPlot.pdf"), accuracy_boxplot, height = 10, width = 10)}

    if(filetag!=""){ggsave(paste0(dir,"Figures and Data/jpg/", filetag, "_AccuracyBoxPlot.jpg"), accuracy_boxplot, height = 10, width = 10)}else
      {ggsave(paste0(dir,"Figures and Data/jpg/AccuracyBoxPlot.jpg"),accuracy_boxplot,height = 10,width = 10)}

    if(filetag!=""){write.csv(final.stats.pofzeds, paste0(dir,"Figures and Data/data/", filetag,"_AccuracyBoxPlotData.csv"), row.names = FALSE, quote = FALSE)}else
      {write.csv(final.stats.pofzeds, paste0(dir,"Figures and Data/data/AccuracyBoxPlotData.csv"), row.names = FALSE, quote = FALSE)}



        ## now - plot the accuracy at all PofZ between 0.5 and 0.99
        testsum <- summarise(group_by(final.stats.2, PofZ, nloci, known), mean(means))
        accuracy_sd_hold <- summarise(group_by(final.stats.2, PofZ, nloci, known), sd(means))
        testsum[5] <- accuracy_sd_hold[4]
        # testsum[4] <- as.numeric(as.character(testsum[4]))
        # testsum[5] <- as.numeric(as.character(testsum[5]))
        colnames(testsum)[4]= "means"
        colnames(testsum)[5]= "sd"
        # testsum$sd <- droplevels(testsum$sd)
        testsum$sdPos <- testsum$means+testsum$sd
        testsum$sdNeg <- testsum$means-testsum$sd
        testsum <- data.frame(testsum)


        ## line plot - accuracy with SD
       accuracy_lineplotSD <- ggplot(testsum) +
           geom_line(aes(x = PofZ, y = means, colour = known), lwd = 1.25) + geom_line(aes(y = sdNeg, x = PofZ, colour = known), linetype = 2) + geom_line(aes(y = sdPos, x = PofZ, colour = known), linetype = 2) +
          facet_grid(.~nloci) + theme(panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),
          legend.position="bottom", strip.background = element_rect(colour = "black", fill = "white")) +
          scale_color_brewer(palette = "Dark2")+
          labs(x = "Critical PofZ Threshold", y = expression("Proportion of Assignments Correct "%+-%"sd"), col="Genotype Frequency Class") + ylim(0, 1)

               #     #Save plot
    if(filetag!=""){ggsave(paste0(dir,"Figures and Data/pdf/", filetag, "_AccuracyLinePlotSD.pdf"), accuracy_lineplotSD, height = 10, width = 10)}else
      {ggsave(paste0(dir, "Figures and Data/pdf/AccuracyLinePlotSD.pdf"), accuracy_lineplotSD, height = 10, width = 10)}

    if(filetag!=""){ggsave(paste0(dir,"Figures and Data/jpg/", filetag, "_AccuracyLinePlotSD.jpg"), accuracy_lineplotSD, height = 10, width = 10)}else
      {ggsave(paste0(dir,"Figures and Data/jpg/AccuracyLinePlotSD.jpg"),accuracy_lineplotSD, height = 10, width = 10)}

    # if(filetag!=""){write.csv(final.stats.pofzeds, paste0(dir,"Figures and Data/data/", filetag,"_AccuracyBoxPlotData.csv"), row.names = FALSE, quote = FALSE)}else
    #   {write.csv(final.stats.pofzeds, paste0(dir,"Figures and Data/data/AccuracyBoxPlotData.csv"), row.names = FALSE, quote = FALSE)}

        ## line plot - accuracy no SD
       accuracy_lineplot <-  ggplot(testsum) +
           geom_line(aes(x = PofZ, y = means, colour = known), lwd = 1.25) + scale_color_brewer(palette = "Dark2")+
          facet_grid(.~nloci) + theme(panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),
          legend.position="bottom", strip.background = element_rect(colour = "black", fill = "white")) +
          labs(x = "Critical PofZ Threshold", y = "Proportion of Assignments Correct ", col="Genotype Frequency Class") + ylim(0, 1)


       if(filetag!=""){ggsave(paste0(dir,"Figures and Data/pdf/", filetag, "_AccuracyLinePlot.pdf"), accuracy_lineplot, height = 10, width = 10)}else
      {ggsave(paste0(dir, "Figures and Data/pdf/AccuracyLinePlot.pdf"), accuracy_lineplot, height = 10, width = 10)}

    if(filetag!=""){ggsave(paste0(dir,"Figures and Data/jpg/", filetag, "_AccuracyLinePlot.jpg"), accuracy_lineplot, height = 10, width = 10)}else
      {ggsave(paste0(dir,"Figures and Data/jpg/AccuracyLinePlotSD.jpg"),accuracy_lineplotSD, height = 10, width = 10)}

    if(filetag!=""){write.csv(testsum, paste0(dir,"Figures and Data/data/", filetag,"_AccuracyLinePlotData.csv"), row.names = FALSE, quote = FALSE)}else
      {write.csv(testsum, paste0(dir,"Figures and Data/data/AccuracyLinePlotData.csv"), row.names = FALSE, quote = FALSE)}




       accuracy_lineplot_ClassFacet_SD <- ggplot(testsum) +
           geom_line(aes(x = PofZ, y = means, colour = nloci), lwd = 1.25) + geom_line(aes(y = sdNeg, x = PofZ, colour = nloci), linetype = 2) + geom_line(aes(y = sdPos, x = PofZ, colour = nloci), linetype = 2) +
          facet_grid(.~known) + theme(panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),
          legend.position="bottom", strip.background = element_rect(colour = "black", fill = "white")) +
          scale_color_brewer(palette = "Dark2")+
          labs(x = "Critical PofZ Threshold", y = expression("Proportion of Assignments Correct "%+-%"sd"), col="Panel Size (Loci)") + ylim(0, 1)

       if(filetag!=""){ggsave(paste0(dir,"Figures and Data/pdf/", filetag, "_AccuracyLinePlot_ClassFacetSD.pdf"), accuracy_lineplot_ClassFacet_SD, height = 10, width = 10)}else
      {ggsave(paste0(dir, "Figures and Data/pdf/AccuracyLinePlot_ClassFacetSD.pdf"), accuracy_lineplot_ClassFacet_SD, height = 10, width = 10)}

    if(filetag!=""){ggsave(paste0(dir,"Figures and Data/jpg/", filetag, "_AccuracyLinePlot_ClassFacetSD.jpg"), accuracy_lineplot_ClassFacet_SD, height = 10, width = 10)}else
      {ggsave(paste0(dir,"Figures and Data/jpg/AccuracyLinePlot_ClassFacetSD.jpg"),accuracy_lineplot_ClassFacet_SD, height = 10, width = 10)}

    if(filetag!=""){write.csv(testsum, paste0(dir,"Figures and Data/data/", filetag,"_AccuracyLinePlot_ClassFacetSDData.csv"), row.names = FALSE, quote = FALSE)}else
      {write.csv(testsum, paste0(dir,"Figures and Data/data/AccuracyLinePlot_ClassFacetSDData.csv"), row.names = FALSE, quote = FALSE)}


      accuracy_lineplot_ClassFacet <-  ggplot(testsum) +
           geom_line(aes(x = PofZ, y = means, colour = nloci), lwd = 1.25) +
          facet_grid(.~known) + theme(panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),
          legend.position="bottom", strip.background = element_rect(colour = "black", fill = "white")) +
          scale_color_brewer(palette = "Dark2")+
          labs(x = "Critical PofZ Threshold", y = expression("Proportion of Assignments Correct "%+-%"sd"), col="Panel Size (Loci)") + ylim(0, 1)

       if(filetag!=""){ggsave(paste0(dir,"Figures and Data/pdf/", filetag, "_AccuracyLinePlot_ClassFacet.pdf"), accuracy_lineplot_ClassFacet, height = 10, width = 10)}else
      {ggsave(paste0(dir, "Figures and Data/pdf/AccuracyLinePlot_ClassFacet.pdf"), accuracy_lineplot_ClassFacet, height = 10, width = 10)}

    if(filetag!=""){ggsave(paste0(dir,"Figures and Data/jpg/", filetag, "_AccuracyLinePlot_ClassFacet.jpg"), accuracy_lineplot_ClassFacet, height = 10, width = 10)}else
      {ggsave(paste0(dir,"Figures and Data/jpg/AccuracyLinePlot_ClassFacet.jpg"),accuracy_lineplot_ClassFacet, height = 10, width = 10)}

    if(filetag!=""){write.csv(testsum, paste0(dir,"Figures and Data/data/", filetag,"_AccuracyLinePlot_ClassFacetData.csv"), row.names = FALSE, quote = FALSE)}else
      {write.csv(testsum, paste0(dir,"Figures and Data/data/AccuracyLinePlot_ClassFacetData.csv"), row.names = FALSE, quote = FALSE)}


        ##################
        ## TYPE I ERROR ##
        ##################


  ### Calculate Type I Error - The number of individuals wrongly called hybrid over the actual number of purebreds in the sample

  writeLines("
           Calculating Type I Error
           ")


         bad.res.out <- NULL
         good.res.out <- NULL
    ### TYPE I ERROR
       TypeIout <- NULL
            stat1_check_progress <- txtProgressBar(min = 0.5, max = 0.99, style = 3)
        tempType1 <- temp
        final.stats.1 <- NULL
        for(r in 50:99/100){
            # i = 0.75

            setTxtProgressBar(stat1_check_progress, r)
            tempType1$isgood <- "YES" ## set all isgood to "YES" - default is correct <- will switch to "NO" if PofZ value < critical value
            tempType1$domatch <- "NO" ## set all domatch to "NO" -  default is incorrect <- will swith to "YES" if they match
                  ##

              for(s in 1:nrow(tempType1)){

                # # print(j)
                ## I think order here is important - evaluate match, then if it is greater than critical PofZ - but I can't remember why
                  ##
                if(as.character(tempType1[s, "known"]) == as.character(tempType1[s, "max.class"])){ ## evaluate match
                  tempType1[s, "domatch"] = "YES"
                }

                if(tempType1[s, 5:10][which.max(tempType1[s, 5:10])] < r){ ## evaluate PofZ threshold
                  tempType1[s, "isgood"] = "NO"
                }

              } ## END s Loop

          out.stats.Type1 = NULL
          for(u in 1:length(loci.groups)){

            test3 <- tempType1
            test3 <- test3[test3$nLoci == loci.groups[u], ]
            test3 <- test3[(test3$known == "P1" | test3$known == "P2"), ]

            goodpos <- test3$max.class %in% c("P1","P2")  & test3$isgood == "YES"
            # badpos <- test3$max.class !="P1" & test3$max.class != "P2" & test3$isgood == "YES"
            badpos <- !test3$max.class %in% c("P1","P2")  & test3$isgood == "YES"
            ################ RYAN - Can you put piping in here?????
            good.res <- test3[goodpos, ]
            bad.res <- test3[badpos, ]

            # bad.res.out <- rbind(bad.res.out, bad.res)
            # good.res.out <- rbind(good.res.out, good.res)

            numsim <- length(unique(test3$sim))

            temp.hold.matrix <- matrix(nrow = 2, ncol = numsim)
            colnames(temp.hold.matrix) <- unique(test3$sim)

            good.table <- table(good.res$sim)
            good.temp <- t(data.frame(good.table))
            colnames(good.temp) <- good.temp[1,]
            # good.temp <- good.temp[-1, ]

            bad.table <- table(bad.res$sim)
            bad.temp <- t(data.frame(bad.table))
            colnames(bad.temp) <- bad.temp[1,]
            # bad.temp <- bad.temp[-1, ]

            temp.hold.matrix <- matrix(nrow = 2, ncol = numsim)
            temp.hold.matrix[1, ] <- 0
            # temp.hold.matrix[2, ] <- good.table
            temp.hold.matrix[2, ] <- 0

            match.cols.bad <- which(colnames(good.temp)%in%colnames(bad.temp))
            match.cols.good <- which(colnames(good.temp)%in%unique(test3$sim))

            # temp.hold.matrix[1, match.cols.bad] <- bad.temp[2, ]
            # temp.hold.matrix[2, match.cols.good] <- good.temp[2, ]

            bad.table.dims <- dim(bad.temp)[2]

              if(bad.table.dims > 0){

                for(z in 1:length(match.cols.bad)){
                  temp.hold.matrix[1, match.cols.bad[z]] = bad.temp[2, z]
                  } ## END FOR
                    } ## END IF

            good.table.dims <- dim(good.temp)[2]

              if(good.table.dims > 0){

                for(g in 1:max(good.table.dims)){
                  temp.hold.matrix[2, match.cols.good[g]] = good.temp[2, g]
                  } ## END FOR
                    } ## END IF

            colnames(temp.hold.matrix) <- colnames(good.temp)

            typeI <- as.numeric(temp.hold.matrix[1, ])/as.numeric(temp.hold.matrix[2, ])
            typeI <- cbind(loci.groups[u], typeI)
            out.stats.Type1 = rbind(out.stats.Type1, typeI)

            } ## END u loop


        out.stats.Type1 <- data.frame(out.stats.Type1)
        colnames(out.stats.Type1) <- c("Loci", "Prop")

        # TypeIout.hold <- out.stats.Type1
        TypeIout.hold <- data.frame(Loci = out.stats.Type1$Loci, Prop = out.stats.Type1$Prop, PofZ = r)
        TypeIout <- rbind(TypeIout, TypeIout.hold)

        } ### END TYPE I ERROR LOOP


        TypeIout$Loci <- factor(x = TypeIout$Loci, levels = ordered(unique(as.numeric(as.character(TypeIout$Loci)))))

        TypeIout.PofZeds <- TypeIout$PofZ == 0.5 | TypeIout$PofZ == 0.75 | TypeIout$PofZ == 0.9
        TypeIout.PofZeds <- TypeIout[TypeIout.PofZeds, ]
        TypeIout.PofZeds$Prop <- as.numeric(as.character(TypeIout.PofZeds$Prop))
        TypeIout.PofZeds$PofZ <- as.factor(TypeIout.PofZeds$PofZ)

        TypeI_BoxPlot <- ggplot(TypeIout.PofZeds, aes(x = PofZ, y = Prop)) + geom_boxplot(fill = "grey75") + facet_grid(.~Loci) +
          theme(panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),
          legend.position = "none", strip.background = element_rect(colour = "black", fill = "white")) +
          labs(x = "Critical PofZ Threshold", y = "Type I Error Proportion")

    if(filetag!=""){ggsave(paste0(dir,"Figures and Data/pdf/", filetag, "_TypeIBoxPlot.pdf"), TypeI_BoxPlot, height = 10, width = 10)}else
      {ggsave(paste0(dir, "Figures and Data/pdf/TypeIBoxPlot.pdf"), TypeI_BoxPlot, height = 10, width = 10)}

    if(filetag!=""){ggsave(paste0(dir,"Figures and Data/jpg/", filetag, "_TypeIBoxPlot.jpg"), TypeI_BoxPlot, height = 10, width = 10)}else
      {ggsave(paste0(dir,"Figures and Data/jpg/TypeIBoxPlot.jpg"), TypeI_BoxPlot, height = 10, width = 10)}

    if(filetag!=""){write.csv(TypeIout.PofZeds, paste0(dir,"Figures and Data/data/", filetag,"_TypeIBoxPlotData.csv"), row.names = FALSE, quote = FALSE)}else
      {write.csv(TypeIout.PofZeds, paste0(dir,"Figures and Data/data/TypeIBoxPlotData.csv"), row.names = FALSE, quote = FALSE)}


        ## now - plot the TYPE I ERROR at all PofZ between 0.5 and 0.99
        TypeIout$Prop <- as.numeric(as.character(TypeIout$Prop))
        testsum_TypeI <- summarise(group_by(TypeIout, PofZ, Loci), mean(Prop))
        hold.sd <- summarise(group_by(TypeIout, PofZ, Loci), sd(Prop))
        colnames(testsum_TypeI)[3] <- "means"
        testsum_TypeI[,4] <- hold.sd[,3]
        colnames(testsum_TypeI)[4] <- "sd"
        testsum_TypeI$sd <- as.numeric(as.character(testsum_TypeI$sd))

        max.y = max(testsum_TypeI$means + testsum_TypeI$sd)

         ## line plot - Type I Error
         typeI_lineplot <- ggplot(testsum_TypeI) + geom_line(aes(x = PofZ, y = means, colour = Loci), size = 2) + geom_line( aes(y = (means+sd), x = PofZ, colour = Loci), linetype = 2) +
           geom_line( aes(y = (means-sd), x = PofZ, colour = Loci), linetype = 2) + facet_grid(.~Loci) +  ylim(ymin = 0, ymax = max.y) +
           labs(x = "Critical PofZ Threshold", y = expression("Type I Error Proportion "%+-%"sd")) +
           theme(panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),
             legend.position = "none", strip.background = element_rect(colour = "black", fill = "white"), text = element_text(colour = "black"))

    if(filetag!=""){ggsave(paste0(dir,"Figures and Data/pdf/", filetag, "_TypeILinePlot.pdf"), typeI_lineplot, height = 10, width = 10)}else
      {ggsave(paste0(dir, "Figures and Data/pdf/TypeILinePlot.pdf"), typeI_lineplot, height = 10, width = 10)}

    if(filetag!=""){ggsave(paste0(dir,"Figures and Data/jpg/", filetag, "_TypeILinePlot.jpg"), typeI_lineplot, height = 10, width = 10)}else
      {ggsave(paste0(dir,"Figures and Data/jpg/TypeILinePlot.jpg"), typeI_lineplot, height = 10, width = 10)}

    if(filetag!=""){write.csv(testsum_TypeI, paste0(dir,"Figures and Data/data/", filetag,"_TypeILinePlotData.csv"), row.names = FALSE, quote = FALSE)}else
      {write.csv(testsum_TypeI, paste0(dir,"Figures and Data/data/TypeILinePlotData.csv"), row.names = FALSE, quote = FALSE)}

         ############################
         ## MEAN AND SD PofZ SCORE ##
         ############################

       ## average and SD the  replicate runs of each simulation in New Hybrids. Filter is just a holder for the dplyr:: call
      sim_data <- as.data.frame(dplyr::filter(output)%>%group_by(nLoci,sim,Indv)%>%summarise(Pure1_sd=sd(Pure1),Pure1=mean(Pure1),
                                                                       Pure2_sd=sd(Pure2),Pure2=mean(Pure2),
                                                                       F1_sd=sd(F1),F1=mean(F1),
                                                                       F2_sd=sd(F2),F2=mean(F2),
                                                                       BC1_sd=sd(BC1),BC1=mean(BC1),
                                                                       BC2_sd=sd(BC2),BC2=mean(BC2))%>%ungroup())

    #pull out just the means of the replicates
      sim_means <- sim_data[,-grep("_sd",colnames(sim_data))]

    ## assign the classes to the data
      classvec <- rep(c("Pure1","Pure2","F1","F2","BC1","BC2"),times=samplesize)
      classvec <- rep(classvec,times=nrow(sim_means)/length(classvec))
      sim_means$class <- classvec

    #Compare the simulations using boxplots
      boxdata <- NULL
      for (i in unique(sim_means$nLoci))
        {
        for (j in unique(sim_means$class))
          {
          temp <- filter(sim_means,class==j,nLoci==i)
          tout <- data.frame(sim=temp$sim,Indv=temp$Indv,class=j,nLoci=i,value=temp[,j])
          boxdata <- rbind(boxdata,tout)
          } ## END J Loop
        } ## END I LOOP

    boxdata$nLoci=factor(boxdata$nLoci)

    # Create pot
    meanPofZ_AllClasses_BoxPlot <-
      ggplot(boxdata,aes(x=nLoci,y=value,fill=sim))+geom_boxplot(alpha = 0.8, outlier.size = 0) + facet_wrap(~class,nrow=3,scales="free_y")+
      labs(y="PofZ Score",x="Panel Size (Loci)") + scale_fill_manual(values=c("grey75","grey75","grey75"))+
      theme(panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),
  legend.position = "none", strip.background = element_rect(colour = "black", fill = "white"), text = element_text(colour = "black"))

    #save plot
    if(filetag!=""){ggsave(paste0(dir,"Figures and Data/pdf/",filetag,"_meanPofZ_AllClasses_BoxPlot.pdf"),  meanPofZ_AllClasses_BoxPlot, height = 8, width = 10)}else
    {ggsave(paste0(dir,"Figures and Data/pdf/meanPofZ_AllClasses_BoxPlot.pdf"),  meanPofZ_AllClasses_BoxPlot, height = 8, width = 10)}

    if(filetag!=""){ggsave(paste0(dir,"Figures and Data/jpg/",filetag,"_meanPofZ_AllClasses_BoxPlot.jpg"),  meanPofZ_AllClasses_BoxPlot, height = 8, width = 10)}else
    {ggsave(paste0(dir,"Figures and Data/jpg/meanPofZ_AllClasses_BoxPlot.jpg"),  meanPofZ_AllClasses_BoxPlot, height = 8, width = 10)}

    if(filetag!=""){write.csv(boxdata, paste0(dir,"Figures and Data/data/",filetag,"_meanPofZ_AllClasses_BoxPlotData.csv"), row.names = FALSE)}else
    {write.csv(boxdata, paste0(dir,"Figures and Data/data/meanPofZ_AllClasses_BoxPlotData.csv"), row.names = FALSE)}

    #Combined loci
    sim_means2 <- sim_means
    sim_means2$hybrid <- rowSums(sim_means2[,c("F1","F2","BC1","BC2")])
    sim_means2[which(sim_means2$class =="Pure1"),"hybrid"]=sim_means2[which(sim_means2$class =="Pure1"),"Pure1"] #add values of the Pure
    sim_means2[which(sim_means2$class =="Pure2"),"hybrid"]=sim_means2[which(sim_means2$class =="Pure2"),"Pure2"] #add values of the Pure

    sim_means2$hclass <- "Hybrid"
    sim_means2[which(sim_means$class=="Pure1"),"hclass"] <- "Pure1"
    sim_means2[which(sim_means$class=="Pure2"),"hclass"] <- "Pure2"

    sim_means2$hclass <- factor(sim_means2$hclass,levels=c("Pure1","Pure2","Hybrid"))

    meanPofZ_PureHyb_BoxPlot <-
      ggplot(sim_means2,aes(x=factor(nLoci),y=hybrid,fill=sim))+geom_boxplot(alpha=0.8,outlier.size = 0)+theme_bw()+facet_wrap(~hclass,nrow=3,scales="free_y")+
      labs(y="PofZ Score",x="Panel Size (Loci)") + scale_fill_manual(values=c("grey75","grey75","grey75"))+
      theme(panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),
  legend.position = "none", strip.background = element_rect(colour = "black", fill = "white"), text = element_text(colour = "black"))

    #Save plot
    if(filetag!=""){ggsave(paste0(dir,"Figures and Data/pdf/",filetag,"_meanPofZ_PureHyb_BoxPlot.pdf"), meanPofZ_PureHyb_BoxPlot, height = 8, width = 8)}else
    {ggsave(paste0(dir,"Figures and Data/pdf/meanPofZ_PureHyb_BoxPlot.pdf"), meanPofZ_PureHyb_BoxPlot, height = 8, width = 8)}

    if(filetag!=""){ggsave(paste0(dir,"Figures and Data/jpg/",filetag,"_meanPofZ_PureHyb_BoxPlot.jpg"), meanPofZ_PureHyb_BoxPlot, height = 8, width = 8)}else
    {ggsave(paste0(dir,"Figures and Data/jpg/meanPofZ_PureHyb_BoxPlot.jpg"),meanPofZ_PureHyb_BoxPlot,height = 8,width = 8)}

    if(filetag!=""){write.csv(sim_means2, paste0(dir,"Figures and Data/data/",filetag,"_meanPofZ_PureHyb_BoxPlotData.csv"), row.names = FALSE)}else
    {write.csv(sim_means2, paste0(dir,"Figures and Data/data/meanPofZ_PureHyb_BoxPlotData.csv"), row.names = FALSE)}


    ###########################
    ### CALCULATE EFFICIENCY ##
    ###########################

    ## Look at assignment success as a function of threshold probability
    num.sim <- length(which(sim_means$sim=="S1"))/6/length(unique(sim_means$nLoci))

    ## find he dim of each class in a given sim (length = sum n_class$n
    classvec2 <- rep(c("Pure1","Pure2","F1","F2","BC1","BC2"),times=samplesize)

    ProbOutput <- NULL
    for (s in unique(sim_means$nLoci)){
        lsub <- filter(sim_means,nLoci == s)
        for(i in unique(sim_means$sim)){
          tempsub <- filter(lsub,sim==i)
            for(q in 50:99/100){ # probability of 50 - 99%

              p1.p <- length(which(tempsub[which(classvec2=="Pure1"),"Pure1"] > q))/samplesize[1]
              p2.p <- length(which(tempsub[which(classvec2=="Pure2"),"Pure2"] > q))/samplesize[2]
              F1.p <- length(which(tempsub[which(classvec2=="F1"),"F1"] > q))/samplesize[3]
              F2.p <- length(which(tempsub[which(classvec2=="F2"),"F2"] > q))/samplesize[4]
              BC1.p <- length(which(tempsub[which(classvec2=="BC1"),"BC1"] > q))/samplesize[5]
              BC2.p <- length(which(tempsub[which(classvec2=="BC2"),"BC2"] > q))/samplesize[6]
              tempout <- data.frame(nLoci=s,sim=i,level=q,prob=c(p1.p, p2.p,F1.p,F2.p,BC1.p,BC2.p),
                                    class=c("Pure1","Pure2","F1","F2","BC1","BC2"))
              ProbOutput <- rbind(ProbOutput,tempout)

            } # end q loop
          } # end i loop
        } # end s loop

    #combined hybrid probabilities _______ EFFICIENCY?????
    ProbOutput2 <- NULL
    for (s in unique(sim_means$nLoci)){
      lsub <- filter(sim_means,nLoci == s)
      for(i in unique(sim_means$sim)){
        tempsub <- filter(lsub,sim==i)
        tempsub$phyb <- rowSums(tempsub[,c("F1","F2","BC1","BC2")])

        for(q in 50:99/100){ # probability of 50 - 99%
          p1.p <- length(which(tempsub[which(classvec2=="Pure1"),"Pure1"] > q))/samplesize[1]
          p2.p <- length(which(tempsub[which(classvec2=="Pure2"),"Pure2"] > q))/samplesize[2]
          Hybrid <- length(which(tempsub[which(classvec2%in%c("F1","F2","BC1","BC2")),"phyb"]>q))/sum(samplesize[3:6])
          tempout <- data.frame(nLoci=s,sim=i,level=q,prob=c(p1.p, p2.p,Hybrid),
                                class=c("Pure1","Pure2","Hybrid"))
          ProbOutput2 <- rbind(ProbOutput2,tempout)
        } # end q loop
      } # end i loop
    } # end s loop

    # get the mean and standard error for the estimates of assignment succes based on NH probabilty among simulations
      FinalData <- data.frame(ProbOutput%>%group_by(nLoci,level,class)%>%summarise(mprob = mean(prob,na.rm=T),
                                                                            sdprob = sd(prob,na.rm=T))%>%ungroup())
      FinalData$class <- factor(FinalData$class, levels=c("Pure1","Pure2","F1","F2","BC1","BC2")) # NH class

    # set plotting levels
      FinalData$group <- "Pure"
      FinalData[which(FinalData$class %in% c("BC1","BC2")),"group"] <- "Back-cross"
      FinalData[which(FinalData$class %in% c("F1","F2")),"group"] <- "Generational hybrids"

      FinalData$group <-  factor(FinalData$group,levels=c("Pure","Generational hybrids","Back-cross"))
      FinalData$class <- factor(FinalData$class,levels=c("Pure1","Pure2","F1","F2","BC1","BC2"))

    #plot the class groupings
      Efficiency_AllClass_LinePlot <-
        ggplot(FinalData, aes(x = level, y = mprob, col = class)) + geom_line(lwd = 1.25) + facet_grid( group ~ nLoci, scales = "free_y") +
        theme(panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),legend.position="bottom",strip.background = element_rect(fill="white",colour = "black"), text = element_text(colour = "black"))+
        scale_color_brewer(palette = "Dark2")+
        labs(x="Critical PofZ Threshold",y="Efficiency",col="Genotype Frequency Class")



      #Save plot
      if(filetag!=""){ggsave(paste0(dir,"Figures and Data/pdf/",filetag,"_Efficiency_AllClass_LinePlot.pdf"), Efficiency_AllClass_LinePlot, height = 10, width = 8)} else
      {ggsave(paste0(dir,"Figures and Data/pdf/Efficiency_AllClass_LinePlot.pdf"), Efficiency_AllClass_LinePlot, height = 10, width = 8)}

      if(filetag!=""){ggsave(paste0(dir,"Figures and Data/jpg/",filetag,"_Efficiency_AllClass_LinePlot.jpg"), Efficiency_AllClass_LinePlot, height = 10, width = 8)} else
      {ggsave(paste0(dir,"Figures and Data/jpg/Efficiency_AllClass_LinePlot.jpg"), Efficiency_AllClass_LinePlot, height = 10, width = 8)}

      if(filetag!=""){write.csv(FinalData, paste0(dir,"Figures and Data/data/", filetag ,"_Efficiency_AllClass_LinePlotData.csv"), row.names = FALSE)}else
    {write.csv(FinalData, paste0(dir,"Figures and Data/data/Efficiency_AllClass_LinePlotData.csv"), row.names = FALSE)}

      #ComboHybrids ------------
      FinalData2 <- data.frame(ProbOutput2%>%group_by(nLoci,level,class)%>%summarise(mprob = mean(prob,na.rm=T),
                                                                               sdprob = sd(prob,na.rm=T))%>%ungroup())
      FinalData2$class <- factor(FinalData2$class, levels=c("Pure1","Pure2","Hybrid")) # set plotting levels

      Efficiency_PureHyb_LinePlot <-
        ggplot(FinalData2)+
        geom_line(aes(x=level,y=mprob,col=class),lwd=1.25)+
        geom_line(aes(x=level,y=mprob+sdprob,col=class),lty=2)+
        geom_line(aes(x=level,y=mprob-sdprob,col=class),lty=2)+
        facet_grid(~nLoci)+
        theme(panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),legend.position="bottom",strip.background = element_rect(fill="white",colour = "black"), text = element_text(colour = "black"))+
        scale_color_brewer(palette = "Dark2")+
        labs(x="Critical PofZ Threshold",y=expression("Efficiency "%+-%"sd"),col="Genotype Frequency Class") + ylim(0, 1)

      if(filetag!=""){ggsave(paste0(dir,"Figures and Data/pdf/",filetag,"_Efficiency_PureHyb_LinePlot.pdf"), Efficiency_PureHyb_LinePlot, height = 8, width = 10)} else
      {ggsave(paste0(dir,"Figures and Data/pdf/Efficiency_PureHyb_LinePlot.pdf"), Efficiency_PureHyb_LinePlot, height = 8, width = 10)}

      if(filetag!=""){ggsave(paste0(dir,"Figures and Data/jpg/",filetag,"_Efficiency_PureHyb_LinePlot.jpg"), Efficiency_PureHyb_LinePlot, height = 8, width = 10)} else
      {ggsave(paste0(dir,"Figures and Data/jpg/Efficiency_PureHyb_LinePlot.jpg"), Efficiency_PureHyb_LinePlot,height = 8,width = 10)}

      if(filetag!=""){write.csv(FinalData2, paste0(dir,"Figures and Data/data/", filetag ,"_Efficiency_PureHyb_LinePlotData.csv"), row.names = FALSE)}else
      {write.csv(FinalData2, paste0(dir,"Figures and Data/data/Efficiency_PureHyb_LinePlotData.csv"), row.names = FALSE)}

    #plot if no threshold specified
      if(addThresh){
        Efficiency_AllClass_LinePlot_ClassFacet_thresh <-
          ggplot(data=FinalData)+geom_line(aes(x=level,y=mprob,col=factor(nLoci)),lwd=1.25)+
          geom_line(aes(x=level,y=mprob+sdprob,col=factor(nLoci)),lty=2)+
          geom_line(aes(x=level,y=mprob-sdprob,col=factor(nLoci)),lty=2)+
          facet_wrap(~class,nrow=3,scales="free_y")+theme(strip.background = element_rect(fill="white",colour = "black"))+
          scale_color_brewer(palette = "Dark2")+
          theme(panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),legend.position="bottom",strip.background = element_rect(fill="white",colour = "black"), text = element_text(colour = "black"))+
          labs(x="Critical PofZ Threshold",y=expression("Efficiency "%+-%"sd"),col="Panel Size (Loci)")+geom_vline(xintercept = Thresholds, lty=2) +
          guides(colour = guide_legend(override.aes = list(shape = 15))) +
          ylim(0, 1)

        if(filetag!=""){ggsave(paste0(dir,"Figures and Data/pdf/",filetag,"_Efficiency_AllClass_LinePlot_ClassFacet_thresh.pdf"), Efficiency_AllClass_LinePlot_ClassFacet_thresh, height = 8, width = 10)} else
        {ggsave(paste0(dir,"Figures and Data/pdf/Efficiency_AllClass_LinePlot_ClassFacet_thresh.pdf"), Efficiency_AllClass_LinePlot_ClassFacet_thresh, height = 8, width = 10)}

        if(filetag!=""){ggsave(paste0(dir,"Figures and Data/jpg/",filetag,"_Efficiency_AllClass_LinePlot_ClassFacet_thresh.jpg"), Efficiency_AllClass_LinePlot_ClassFacet_thresh, height = 8, width = 10)} else
        {ggsave(paste0(dir,"Figures and Data/jpg/Efficiency_AllClass_LinePlot_ClassFacet_thresh.jpg"), Efficiency_AllClass_LinePlot_ClassFacet_thresh, height = 8, width = 10)}

        if(filetag!=""){write.csv(FinalData, paste0(dir,"Figures and Data/data/", filetag ,"_Efficiency_AllClass_LinePlot_ClassFacet_threshData.csv"), row.names = FALSE)}else
        {write.csv(FinalData2, paste0(dir,"Figures and Data/data/Efficiency_AllClass_LinePlot_ClassFacet_threshData.csv"), row.names = FALSE)}

        }

      if(!addThresh){
        Efficiency_AllClass_LinePlot_ClassFacet_NOthresh <-
          ggplot(data=FinalData)+geom_line(aes(x=level,y=mprob,col=factor(nLoci)),lwd=1.25)+
          geom_line(aes(x=level,y=mprob+sdprob,col=factor(nLoci)),lty=2)+
          geom_line(aes(x=level,y=mprob-sdprob,col=factor(nLoci)),lty=2)+
           facet_wrap(~class,nrow=3,scales="free_y")+theme(strip.background = element_rect(fill="white",colour = "black"))+
        scale_color_brewer(palette = "Dark2")+
        theme(panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),legend.position="bottom",strip.background = element_rect(fill="white",colour = "black"), text = element_text(colour = "black"))+
        labs(x="Critical PofZ Threshold",y=expression("Efficiency "%+-%"sd"),col="Panel Size (Loci)")+
        ylim(0, 1)

        if(filetag!=""){ggsave(paste0(dir,"Figures and Data/pdf/",filetag,"_Efficiency_AllClass_LinePlot_ClassFacet_NOthresh.pdf"), Efficiency_AllClass_LinePlot_ClassFacet_NOthresh, height = 8, width = 10)} else
        {ggsave(paste0(dir,"Figures and Data/pdf/Efficiency_AllClass_LinePlot_ClassFacet_NOthresh.pdf"), Efficiency_AllClass_LinePlot_ClassFacet_NOthresh, height = 8, width = 10)}

        if(filetag!=""){ggsave(paste0(dir,"Figures and Data/jpg/",filetag,"_Efficiency_AllClass_LinePlot_ClassFacet_NOthresh.jpg"), Efficiency_AllClass_LinePlot_ClassFacet_NOthresh, height = 8, width = 10)} else
        {ggsave(paste0(dir,"Figures and Data/jpg/Efficiency_AllClass_LinePlot_ClassFacet_NOthresh.jpg"), Efficiency_AllClass_LinePlot_ClassFacet_NOthresh, height = 8, width = 10)}

        if(filetag!=""){write.csv(FinalData, paste0(dir,"Figures and Data/data/", filetag ,"_Efficiency_AllClass_LinePlot_ClassFacet_NOthreshData.csv"), row.names = FALSE)}else
        {write.csv(FinalData, paste0(dir,"Figures and Data/data/Efficiency_AllClass_LinePlot_ClassFacet_NOthreshData.csv"), row.names = FALSE)}


        }


      ## combined hybrids

      if(addThresh){
        Efficiency_PureHyb_LinePlot_ClassFacet_thresh <-
          ggplot(data=FinalData2)+geom_line(aes(x=level,y=mprob,col=factor(nLoci)),lwd=1.25)+
          geom_line(aes(x=level,y=mprob+sdprob,col=factor(nLoci)),lty=2)+
          geom_line(aes(x=level,y=mprob-sdprob,col=factor(nLoci)),lty=2)+
          theme_bw()+facet_grid(~class)+theme(strip.background = element_rect(fill="white",colour = "black"))+
          scale_color_brewer(palette = "Dark2")+
          theme(panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),legend.position="bottom",strip.background = element_rect(fill="white",colour = "black"), text = element_text(colour = "black")) +
          labs(x="Critical PofZ Threshold",y=expression("Efficiency "%+-%"sd"),col="Panel Size (Loci)")+
          ylim(0, 1) + geom_vline(xintercept = Thresholds, lty=2)

        if(filetag!=""){ggsave(paste0(dir,"Figures and Data/pdf/",filetag,"_Efficiency_PureHyb_LinePlot_ClassFacet_thresh.pdf"), Efficiency_PureHyb_LinePlot_ClassFacet_thresh, height = 8, width = 10)} else
        {ggsave(paste0(dir,"Figures and Data/pdf/Efficiency_PureHyb_LinePlot_ClassFacet_thresh.pdf"), Efficiency_PureHyb_LinePlot_ClassFacet_thresh, height = 8, width = 10)}

        if(filetag!=""){ggsave(paste0(dir,"Figures and Data/jpg/",filetag,"_Efficiency_PureHyb_LinePlot_ClassFacet_thresh.jpg"), Efficiency_PureHyb_LinePlot_ClassFacet_thresh, height = 8, width = 10)} else
        {ggsave(paste0(dir,"Figures and Data/jpg/Efficiency_PureHyb_LinePlot_ClassFacet_thresh.jpg"), Efficiency_PureHyb_LinePlot_ClassFacet_thresh, height = 8, width = 10)}

        if(filetag!=""){write.csv(FinalData2, paste0(dir,"Figures and Data/data/", filetag ,"_Efficiency_PureHyb_LinePlot_ClassFacet_threshData.csv"), row.names = FALSE)}else
        {write.csv(FinalData2, paste0(dir,"Figures and Data/data/Efficiency_PureHyb_LinePlot_ClassFacet_threshData.csv"), row.names = FALSE)}


        }

      if(!addThresh){
        Efficiency_PureHyb_LinePlot_ClassFacet_NOthresh <-
          ggplot(data=FinalData2)+geom_line(aes(x=level,y=mprob,col=factor(nLoci)),lwd=1.25)+
          geom_line(aes(x=level,y=mprob+sdprob,col=factor(nLoci)),lty=2)+
          geom_line(aes(x=level,y=mprob-sdprob,col=factor(nLoci)),lty=2)+
          theme_bw()+facet_grid(~class)+theme(strip.background = element_rect(fill="white",colour = "black"))+
          scale_color_brewer(palette = "Dark2")+
             theme(panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),legend.position="bottom",strip.background = element_rect(fill="white",colour = "black"), text = element_text(colour = "black")) +
          labs(x="Critical PofZ Threshold",y=expression("Efficiency "%+-%"sd"),col="Panel Size (Loci)")+
          ylim(0, 1)

        if(filetag!=""){ggsave(paste0(dir,"Figures and Data/pdf/",filetag,"_Efficiency_PureHyb_LinePlot_ClassFacet_NOthresh.pdf"), Efficiency_PureHyb_LinePlot_ClassFacet_NOthresh, height = 8, width = 10)} else
        {ggsave(paste0(dir,"Figures and Data/pdf/Efficiency_PureHyb_LinePlot_ClassFacet_NOthresh.pdf"), Efficiency_PureHyb_LinePlot_ClassFacet_NOthresh, height = 8, width = 10)}

        if(filetag!=""){ggsave(paste0(dir,"Figures and Data/jpg/",filetag,"_Efficiency_PureHyb_LinePlot_ClassFacet_NOthresh.jpg"), Efficiency_PureHyb_LinePlot_ClassFacet_NOthresh, height = 8, width = 10)} else
        {ggsave(paste0(dir,"Figures and Data/jpg/Efficiency_PureHyb_LinePlot_ClassFacet_NOthresh.jpg"), Efficiency_PureHyb_LinePlot_ClassFacet_NOthresh, height = 8, width = 10)}

        if(filetag!=""){write.csv(FinalData2, paste0(dir,"Figures and Data/data/", filetag ,"_Efficiency_PureHyb_LinePlot_ClassFacet_NOthreshData.csv"), row.names = FALSE)}else
        {write.csv(FinalData2, paste0(dir,"Figures and Data/data/Efficiency_PureHyb_LinePlot_ClassFacet_NOthreshData.csv"), row.names = FALSE)}
      }



      ## mean plot ----------

      #facet labels
      FinalData$threshold <- paste0(FinalData$level*100,"%")

       Efficiency_ByThreshold_LinePlot_AllClass <-
        ggplot(filter(FinalData,level %in% Thresholds),aes(x=factor(nLoci),y=mprob,col=class,group=class))+
        geom_point(size=2.5)+geom_path(lwd=0.9)+
        geom_errorbar(aes(ymin=mprob-sdprob,ymax=mprob+sdprob),width=0.1)+
        facet_grid(~level)+theme_bw()+
        labs(x="Panel Size (Loci)",y=expression("Efficiency "%+-%"sd"),col="Genotype Frequency Class",group="")+scale_color_brewer(palette = "Dark2")+
        theme(panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),legend.position="bottom",strip.background = element_rect(fill="white",colour = "black"), text = element_text(colour = "black"))

        #Save plot
        if(filetag!=""){ggsave(paste0(dir,"Figures and Data/pdf/",filetag,"_Efficiency_ByThreshold_LinePlot_AllClass.pdf"), Efficiency_ByThreshold_LinePlot_AllClass, height = 8, width = 10)} else
        {ggsave(paste0(dir,"Figures and Data/pdf/Efficiency_ByThreshold_LinePlot_AllClass.pdf"), Efficiency_ByThreshold_LinePlot_AllClass, height = 8, width = 10)}

        if(filetag!=""){ggsave(paste0(dir,"Figures and Data/jpg/",filetag,"_Efficiency_ByThreshold_LinePlot_AllClass.jpg"), Efficiency_ByThreshold_LinePlot_AllClass, height = 8, width = 10)} else
        {ggsave(paste0(dir,"Figures and Data/jpg/Efficiency_ByThreshold_LinePlot_AllClass.jpg"), Efficiency_ByThreshold_LinePlot_AllClass, height = 8, width = 10)}

       if(filetag!=""){write.csv(FinalData, paste0(dir,"Figures and Data/data/", filetag ,"_Efficiency_PureHyb_LinePlot_ClassFacet_NOthreshData.csv"), row.names = FALSE)}else
        {write.csv(FinalData, paste0(dir,"Figures and Data/data/Efficiency_PureHyb_LinePlot_ClassFacet_NOthreshData.csv"), row.names = FALSE)}

        #Combined Hybrids
        FinalData2$threshold <- paste0(FinalData2$level*100,"%")

        Efficiency_ByThreshold_LinePlot_PureHyb <-
          ggplot(filter(FinalData2,level %in% Thresholds),aes(x=factor(nLoci),y=mprob,col=class,group=class))+
          geom_point(size=2.5)+geom_path(lwd=0.9)+
          geom_errorbar(aes(ymin=mprob-sdprob,ymax=mprob+sdprob),width=0.1)+
          facet_grid(~level)+theme_bw()+
          labs(x="Panel Size (Loci)",y=expression("Efficiency "%+-%"sd"),col="Genotype Frequency Class",group="")+scale_color_brewer(palette = "Dark2")+
        theme(panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),legend.position="bottom",strip.background = element_rect(fill="white",colour = "black"), text = element_text(colour = "black"))

        #Save plot
        if(filetag!=""){ggsave(paste0(dir,"Figures and Data/pdf/",filetag,"_Efficiency_ByThreshold_LinePlot_PureHyb.pdf"), Efficiency_ByThreshold_LinePlot_PureHyb, height = 8, width = 10)} else
        {ggsave(paste0(dir,"Figures and Data/pdf/Efficiency_ByThreshold_LinePlot_PureHyb.pdf"), Efficiency_ByThreshold_LinePlot_PureHyb, height = 8, width = 10)}

        if(filetag!=""){ggsave(paste0(dir,"Figures and Data/jpg/",filetag,"_Efficiency_ByThreshold_LinePlot_PureHyb.jpg"), Efficiency_ByThreshold_LinePlot_PureHyb, height = 8, width = 10)} else
        {ggsave(paste0(dir,"Figures and Data/jpg/Efficiency_ByThreshold_LinePlot_PureHyb.jpg"), Efficiency_ByThreshold_LinePlot_PureHyb, height = 8, width = 10)}

        if(filetag!=""){write.csv(FinalData2, paste0(dir,"Figures and Data/data/", filetag ,"_Efficiency_ByThreshold_LinePlot_PureHybData.csv"), row.names = FALSE)}else
        {write.csv(FinalData2, paste0(dir,"Figures and Data/data/Efficiency_ByThreshold_LinePlot_PureHybData.csv"), row.names = FALSE)}


        ###################
        ## TYPE II ERROR ##
        ###################


    writeLines("
           Calculating Type II Error
           ")


            stat11_check_progress <- txtProgressBar(min = 0, max = length(unique(sim_means$nLoci)), style = 3)
            prog_convert <- unique(sim_means$nLoci)


    ## Misclassification 'type II' error ------------
      classnames <- c("Pure1","Pure2","F1","F2","BC1","BC2")
        missout <- NULL
        for (s in unique(sim_means$nLoci)){
          lsub <- filter(sim_means,nLoci == s)
          setTxtProgressBar(stat11_check_progress, which(prog_convert == s))
          for(i in unique(sim_means$sim)){
            tempsub <- filter(lsub,sim==i)
            for(q in 50:99/100){ # probability of 50 - 99%
              tempq <- tempsub
              tempq$missclass <- classnames[apply(tempq[,classnames],1,which.max)] # what is the class of the highest NH probability

              tempq$missval <- 999 #place holder
              for(w in 1:nrow(tempq))
              {
                tempq[w,"missval"] <- tempq[w,classnames[which.max(tempq[w,classnames])]]
              }

              temp1 <- tempq[which(tempq$class!=tempq$missclass & tempq$missval>=q),] #dataset with missclassifications


              dummydf=data.frame(Var1 = c("Pure1","Pure2","F1","F2","BC1","BC2"),dummy=NA) # dataframe for dummy values

              for (z in classnames){
                temp2 <- filter(temp1,class == z)
                if(nrow(temp2)>0){
                  temp3 <- as.data.frame(table(temp2$missclass)/samplesize[which(classnames==z)]) # percentage of samples miss classed to a given class of a given type of class (i)

                  temp4 <-  merge(dummydf,temp3,by="Var1",all.x = TRUE)

                  tempout <- data.frame(nLoci=s,sim=i,level=q,class=z,
                                        mclass_P1=temp4[which(temp4$Var1 == "Pure1"),"Freq"],
                                        mclass_P2=temp4[which(temp4$Var1 == "Pure2"),"Freq"],
                                        mclass_F1=temp4[which(temp4$Var1 == "F1"),"Freq"],
                                        mclass_F2=temp4[which(temp4$Var1 == "F2"),"Freq"],
                                        mclass_BC1=temp4[which(temp4$Var1 == "BC1"),"Freq"],
                                        mclass_BC2=temp4[which(temp4$Var1 == "BC2"),"Freq"])
                } else
                {tempout <- data.frame(nLoci=s,sim=i,level=q,class=z,
                                       mclass_P1=NA,
                                       mclass_P2=NA,
                                       mclass_F1=NA,
                                       mclass_F2=NA,
                                       mclass_BC1=NA,
                                       mclass_BC2=NA)}
                missout <- rbind(missout,tempout)

              } # end of z loop
            } # end q loop
          } # end i loop
        } # end s loop

    #calcluate the means among simulations
    miss_mean <- data.frame(missout%>%group_by(nLoci,level,class)%>%summarise(mprobP1 = mean(mclass_P1,na.rm=T),
                                                                   sdprobP1 = sd(mclass_P1,na.rm=T),
                                                                   mprobP2 = mean(mclass_P2,na.rm=T),
                                                                   sdprobP2 = sd(mclass_P2,na.rm=T),
                                                                   mprobF1 = mean(mclass_F1,na.rm=T),
                                                                   sdprobF1 = sd(mclass_F1,na.rm=T),
                                                                   mprobF2 = mean(mclass_F2,na.rm=T),
                                                                   sdprobF2 = sd(mclass_F2,na.rm=T),
                                                                   mprobBC1 = mean(mclass_BC1,na.rm=T),
                                                                   sdprobBC1 = sd(mclass_BC1,na.rm=T),
                                                                   mprobBC2 = mean(mclass_BC2,na.rm=T),
                                                                   sdprobBC2 = sd(mclass_BC2,na.rm=T))%>%ungroup())

    miss_mean[is.na(miss_mean)]=NA #replace NaN's with NAs

    #merge with the other data
    FinalData3 <- merge(miss_mean,FinalData,by=c("nLoci","level","class"))

    PlotData <- melt(FinalData3[c("nLoci","level","class","mprobP1","mprobP2","mprobF1","mprobF2","mprobBC1","mprobBC2")],id.vars=c("nLoci","level","class"))
    PlotDatasd <- melt(FinalData3[c("nLoci","level","class","sdprobP1","sdprobP2","sdprobF1","sdprobF2","sdprobBC1","sdprobBC2")],id.vars=c("nLoci","level","class"))
    PlotData$sd <- PlotDatasd$value
    PlotData$variable <- as.character(PlotData$variable)
    PlotData[which(PlotData$variable == "mprobP1"),"variable"]="Pure1"
    PlotData[which(PlotData$variable == "mprobP2"),"variable"]="Pure2"
    PlotData[which(PlotData$variable == "mprobF1"),"variable"]="F1"
    PlotData[which(PlotData$variable == "mprobF2"),"variable"]="F2"
    PlotData[which(PlotData$variable == "mprobBC1"),"variable"]="BC1"
    PlotData[which(PlotData$variable == "mprobBC2"),"variable"]="BC2"
    PlotData$variable=factor(PlotData$variable,levels=c("Pure1","Pure2","F1","F2","BC1","BC2"))

    #Create the plots
    for (i in unique(PlotData$class)){
        temp.plot <- ggplot(filter(PlotData,class==i))+
        geom_line(aes(x=level,y=value,col=variable),lwd=1.25)+
        geom_line(aes(x=level,y=value+sd,col=variable),lty=2)+
        geom_line(aes(x=level,y=value-sd,col=variable),lty=2)+
        facet_grid(~nLoci,scales="free_y")+
        theme_bw()+scale_color_brewer(palette = "Dark2")+
        theme(legend.position="bottom",strip.background = element_rect(fill="white",colour = "black"))+
        labs(x="Probability threshold",y=paste0("Proportion ",i," misassigned \u00B1 sd"),col="Classification")
temp.plot
        # if(filetag!=""){ggsave(paste0(dir,"Figures and Data/pdf/",filetag,"_",i,"_MissAssignment~z-nloci.pdf"),temp.plot,height = 6,width = 8)} else
        # {ggsave(paste0(dir,paste0("Figures and Data/pdf/",i,"_MissAssignment~z-nloci.pdf"),temp.plot,height = 6,width = 8))}
        #
        # if(filetag!=""){ggsave(paste0(dir,"Figures and Data/jpg/",filetag,"_",i,"_MissAssignment~z-nloci.jpg"),temp.plot,height = 6,width = 8)} else
        # {ggsave(paste0(dir,paste0("Figures and Data/jpg/",i,"_MissAssignment~z-nloci.jpg"),temp.plot,height = 6,width = 8))}

    }



    ## Create summary booklet

    if(filetag!=""){
      pdf(file = paste0(dir,"Figures and Data/pdf/",filetag,"_OutputBooklet.pdf"))
      print(p1);print(h1)
      print(p3);print(h3)
      print(p4);print(h4)
      print(p5);print(h5)
      for (i in unique(PlotData$class)){
        temp.plot <- ggplot(filter(PlotData,class==i))+
          geom_line(aes(x=level,y=value,col=variable),lwd=1.25)+
          geom_line(aes(x=level,y=value+sd,col=variable),lty=2)+
          geom_line(aes(x=level,y=value-sd,col=variable),lty=2)+
          facet_grid(~nLoci,scales="free_y")+
          theme_bw()+scale_color_brewer(palette = "Dark2")+
          theme(legend.position="bottom",strip.background = element_rect(fill="white",colour = "black"))+
          labs(x="Probability threshold",y=paste0("Proportion ",i," misassigned \u00B1 sd"),col="Classification")
          print(temp.plot)
          }
      dev.off()
    } else {
      pdf(file = paste0(dir,"Figures and Data/pdf/",filetag,"_OutputBooklet.pdf"))
      print(p1);print(h1)
      print(p3);print(h3)
      print(p4);print(h4)
      print(p5);print(h5)
      for (i in unique(PlotData$class)){
        temp.plot <- ggplot(filter(PlotData,class==i))+
          geom_line(aes(x=level,y=value,col=variable),lwd=1.25)+
          geom_line(aes(x=level,y=value+sd,col=variable),lty=2)+
          geom_line(aes(x=level,y=value-sd,col=variable),lty=2)+
          facet_grid(~nLoci,scales="free_y")+
          theme_bw()+scale_color_brewer(palette = "Dark2")+
          theme(legend.position="bottom",strip.background = element_rect(fill="white",colour = "black"))+
          labs(x="Probability threshold",y=paste0("Proportion ",i," misassigned \u00B1 sd"),col="Classification")
        print(temp.plot)
      }
      dev.off()
    }

    ## clean workspace
    rm(list=setdiff(ls(), c("p1","p3","p4","p5","h1","h3","h4","h5",
                            "PlotData","boxdata","FinalData","FinalData2","FinalData3","sim_means2","Thresholds","filetag","dir")))

    #Save workspace image
    if(filetag!="")
    {save(list = ls(envir = environment(), all.names = TRUE),
         file = paste0(dir,"Figures and Data/data/",filetag,"_WorkSpace.RData"),
         envir = environment())}else
         {save(list = ls(envir = environment(), all.names = TRUE),
                file = paste0(dir,"Figures and Data/data/WorkSpace.RData"),
               envir = environment())}


} #end function
