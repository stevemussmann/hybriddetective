#' @name hybridPowerComp3
#' @title Assignment power comparison among different SNP subsets using NewHybrids simulated datasets3.
#' @description Evaluates the accuracy with which NewHybrids assigns individuals of known hybrid class to the correct hybrid class in simulated datasets at varying levels of stringency (PofZ). The code will write graphical and numerical results to the directory provided by the user.
#' @param dir File path to the directory in which the NewHybrids results (in individual folders as returned by parallelNH_XX) are located.
#' @param filetag An optional character vector to be applied as the name of the outputs.
#' @param Thresholds A vector of thresholds which will be added to the plots showing the assignment success for different levels of probability of a given class estimated by NewHybrids. Default is (NULL) so if nothing is specified it will not add this to the output plots (success ~ threshold by class).
#' @param samplesize The number of individuals per NewHybrids class. By (default: NULL) this data will be extracted from the "*individuals.txt" output from parallelnewhybrids if present in the same folder as the PofZ file. This can also explicitly defined as a vector (6 values corresponding to # in P1,P2,F1,F2,BC1,BC2) or a path to a *_Individuals.txt.
#' @param CT The threshold posterior probability of assignment (PofZ) to F2 above which Pure Population 1 or Pure Population 2 individuals are flagged to indicate possible non-convergence. The default is 0.1.
#' @param CTI The proportion of individuals in either Pure Population 1 OR Pure Population 2 allowed to exceed the F2 assignment threshold (PofZCutOff). The default is 0.5.
#' @rdname hybridpowercomp3
#' @import ggplot2
#' @import magrittr
#' @importFrom dplyr filter summarise ungroup group_by do
#' @importFrom grid arrow unit
#' @importFrom stringr str_extract
#' @importFrom reshape2 melt
#' @importFrom  scales alpha
#' @export
#'
hybridPowerComp3 <-function(dir,filetag="",Thresholds=c(0.5,0.6,0.7,0.8,0.9),addThresh=FALSE,samplesize=NULL,CT=0.1,CTI=0.5){



 # library(ggplot2)
 # library(magrittr)
 # library(dplyr)
 # library(stringr)
 # library(reshape2)
 # library(grid)
 # library(scales)
 # library(hybriddetective)
 #
 #
 #
 # dir = "~/Desktop/DFO Aquaculture Interaction/South West Rivers Analysis/South Coast - Proportional Sampling Analysis/NL South West Fixed Linkages NH/NH.Results/"
 #  filetag=""
 #  Thresholds=c(0.5,0.6,0.7,0.8,0.9)
 #  addThresh=FALSE
 #  samplesize=NULL
 #  CT=0.1
 #  CTI=0.5

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

          ### ANALYSIS REPs - if used hybriddetective to simulate the data, unique replicates will have different numbers of loci, a unique simulation number, and a replicate nubmer
              tempfiles <- list.files(paste0(filedir,i)) ## Get names of all files in the directory
              pzfile <- tempfiles[grep("PofZ",tempfiles)] ## Which of the files is the PofZ file?
              tempfile <- read.table(paste0(filedir,i,"/",pzfile),head=T) ### Read the PofZ file in

              ## Each analysis must have an accompanying LociAndAlleles file <- use this to figure out how many loci there are
              LociandAlleles <- tempfiles[grep("LociAndAlleles", tempfiles)]
              LandAfile <- readChar(paste0(filedir, i, "/", LociandAlleles), file.info(paste0(filedir, i, "/", LociandAlleles))$size)
              numLociExt <- stringr::str_extract(string = LandAfile, pattern = paste0("from ", "[:digit:]{1,5}", " loci")) ### find the string with the number of loci, extract
              numLociWorking <- gsub(x = numLociExt, pattern = "from ", replacement = "") ## remove "from"
              numLociWorking <- as.numeric(gsub(x = numLociWorking, pattern = " loci", replacement = "")) ### This is how many loci there are

              #identify the simulation and repeat info
              S_ident <- gsub("_","",str_extract(pzfile,paste0("_S","[:digit:]{1}","R","[:digit:]{1}","_"))) ### if used hybriddetective, will have S#_R# - exctract
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
    tempinds <- hybriddetective::n_class(paste0(filedir, lfiles[1], "/", IndividualsPath))
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
        } ## END K Loop


          ########################
          ## Calculate ACCURACY ##
          ########################

          #### CHECK HERE - ACCURACY LOOP - Can speed up???

          writeLines("Calculating Accuracy
          ")

          ## ACCURACY =  number assigned correctly / total number assigned  -> for each category

          temp2 <- temp
          tempAccuracy <- temp
          tempAccuracy$isgood=TRUE ### Set all isgood to TRUE - isgood is a check if the assigned PofZ > the critical PofZ - will change to FALSE in function
          tempAccuracy$domatch=FALSE ### Set all domatch to FALSE - domatch is check if the assigned class == known class. Will change to TRUE in function


          #Create long form data for dplyr loop
          PofZVector_Accuracy <- rep(50:99/100,each=nrow(tempAccuracy)) #vector of PofZs

          tempAccuracyLong <- do.call("rbind", replicate(length(50:99), tempAccuracy, simplify = FALSE))
          tempAccuracyLong$group <- as.character(PofZVector_Accuracy)
          tempAccuracyLong$pofz <- PofZVector_Accuracy

          AccuracyData <- tempAccuracyLong%>%
                          dplyr::group_by(pofz,nLoci,max.class)%>%
                          dplyr::do(accuracyfunction(.))%>%
                          dplyr::ungroup()%>%data.frame()

          #Set up plotting data
          AccuracyData$nLoci <- factor(x = AccuracyData$nLoci, levels = ordered(unique(as.numeric(AccuracyData$nLoci))))
          AccuracyData$max.class <- factor(x = AccuracyData$max.class, levels = c("P1", "P2", "F1", "F2", "BC1", "BC2"))


          accuracy_boxplot <-
            ggplot(filter(AccuracyData,pofz %in% c(0.5,0.75,0.9)), aes(x = max.class, y = means, fill = max.class)) +
            geom_boxplot() + facet_grid(pofz~nLoci) +
            labs(x = "Genotype Frequency Class", y = "Proportion of Assignments Correct") +
            scale_fill_brewer(palette = "Dark2") +
            theme(panel.background = element_rect(fill = "white", colour = "black"),
              plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),
              legend.position = "none", strip.background = element_rect(, colour = "black", fill = "white"),
              strip.text.x = element_text(colour = "black"), strip.text.y = element_text(colour = "black"))

          #Save plot and data
          if(filetag!=""){ggsave(paste0(dir,"Figures and Data/pdf/", filetag, "_AccuracyBoxPlot.pdf"), accuracy_boxplot, height = 10, width = 10)}else
            {ggsave(paste0(dir, "Figures and Data/pdf/AccuracyBoxPlot.pdf"), accuracy_boxplot, height = 10, width = 10)}

          if(filetag!=""){ggsave(paste0(dir,"Figures and Data/jpg/", filetag, "_AccuracyBoxPlot.jpg"), accuracy_boxplot, height = 10, width = 10)}else
            {ggsave(paste0(dir,"Figures and Data/jpg/AccuracyBoxPlot.jpg"),accuracy_boxplot,height = 10,width = 10)}

          if(filetag!=""){write.csv(AccuracyData, paste0(dir,"Figures and Data/data/", filetag,"_AccuracyBoxPlotData.csv"), row.names = FALSE, quote = FALSE)}else
            {write.csv(AccuracyData, paste0(dir,"Figures and Data/data/AccuracyBoxPlotData.csv"), row.names = FALSE, quote = FALSE)}


          #Summary among simulations --
          SummaryAccuracy <- AccuracyData%>%
                            dplyr::group_by(pofz,nLoci,max.class)%>%
                            dplyr::summarise(mean=mean(means,na.rm=T),
                            sd=sd(means,na.rm=T))%>%
                            dplyr::ungroup()%>%data.frame()

          get.y.min.AccuracyLine <- min((SummaryAccuracy$mean - SummaryAccuracy$sd), na.rm = TRUE)

          ## line plot - accuracy with SD
          accuracy_lineplotSD <-
            ggplot(SummaryAccuracy) +
            geom_line(aes(x = pofz, y = mean, colour = max.class), lwd = 1.25) +
            geom_line(aes(y = mean+sd, x = pofz, colour = max.class), linetype = 2) +
            geom_line(aes(y = mean-sd, x = pofz, colour = max.class), linetype = 2) +
            facet_wrap(~nLoci, ncol = 3) +
            theme(panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),
              legend.position="bottom", strip.background = element_rect(colour = "black", fill = "white")) +
            scale_color_brewer(palette = "Dark2") +
            labs(x = "Critical PofZ Threshold", y = expression("Proportion of Assignments Correct "%+-%"sd"), col="Genotype Frequency Class") +
            ylim(get.y.min.AccuracyLine, 1)

            ##Save plot
            if(filetag!=""){ggsave(paste0(dir,"Figures and Data/pdf/", filetag, "_AccuracyLinePlotSD.pdf"), accuracy_lineplotSD, height = 10, width = 10)}else
              {ggsave(paste0(dir, "Figures and Data/pdf/AccuracyLinePlotSD.pdf"), accuracy_lineplotSD, height = 10, width = 10)}

            if(filetag!=""){ggsave(paste0(dir,"Figures and Data/jpg/", filetag, "_AccuracyLinePlotSD.jpg"), accuracy_lineplotSD, height = 10, width = 10)}else
              {ggsave(paste0(dir,"Figures and Data/jpg/AccuracyLinePlotSD.jpg"),accuracy_lineplotSD, height = 10, width = 10)}

            if(filetag!=""){write.csv(SummaryAccuracy, paste0(dir,"Figures and Data/data/", filetag,"_SummaryAccuracy.csv"), row.names = FALSE, quote = FALSE)}else
              {write.csv(SummaryAccuracy, paste0(dir,"Figures and Data/data/SummaryAccuracy.csv"), row.names = FALSE, quote = FALSE)}

          ## line plot - accuracy no SD
          accuracy_lineplot <-
            ggplot(SummaryAccuracy) +
            geom_line(aes(x = pofz, y = mean, colour = max.class), lwd = 1.25) +
            scale_color_brewer(palette = "Dark2")+
            facet_wrap(~nLoci, ncol = 3) +
            theme(panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),
              legend.position="bottom", strip.background = element_rect(colour = "black", fill = "white")) +
            labs(x = "Critical PofZ Threshold", y = "Proportion of Assignments Correct ", col="Genotype Frequency Class") +
            ylim(get.y.min.AccuracyLine, 1)


          if(filetag!=""){ggsave(paste0(dir,"Figures and Data/pdf/", filetag, "_AccuracyLinePlot.pdf"), accuracy_lineplot, height = 10, width = 10)}else
            {ggsave(paste0(dir, "Figures and Data/pdf/AccuracyLinePlot.pdf"), accuracy_lineplot, height = 10, width = 10)}

          if(filetag!=""){ggsave(paste0(dir,"Figures and Data/jpg/", filetag, "_AccuracyLinePlot.jpg"), accuracy_lineplot, height = 10, width = 10)}else
            {ggsave(paste0(dir,"Figures and Data/jpg/AccuracyLinePlotSD.jpg"),accuracy_lineplotSD, height = 10, width = 10)}

          # if(filetag!=""){write.csv(testsum, paste0(dir,"Figures and Data/data/", filetag,"_AccuracyLinePlotData.csv"), row.names = FALSE, quote = FALSE)}else
          #   {write.csv(testsum, paste0(dir,"Figures and Data/data/AccuracyLinePlotData.csv"), row.names = FALSE, quote = FALSE)}


          accuracy_lineplot_ClassFacet_SD <-
            ggplot(SummaryAccuracy) +
            geom_line(aes(x = pofz, y = mean, colour = nLoci), lwd = 1.25) +
            geom_line(aes(y = mean-sd, x = pofz, colour = nLoci), linetype = 2) +
            geom_line(aes(y = mean+sd, x = pofz, colour = nLoci), linetype = 2) +
            facet_wrap(~max.class, nrow = 3) +
            theme(panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),
              legend.position="bottom", strip.background = element_rect(colour = "black", fill = "white")) +
            scale_color_brewer(palette = "Dark2")+
            labs(x = "Critical PofZ Threshold", y = expression("Proportion of Assignments Correct "%+-%"sd"), col="Panel Size (Loci)") +
            ylim(get.y.min.AccuracyLine, 1)

            if(filetag!=""){ggsave(paste0(dir,"Figures and Data/pdf/", filetag, "_AccuracyLinePlot_ClassFacetSD.pdf"), accuracy_lineplot_ClassFacet_SD, height = 10, width = 10)}else
              {ggsave(paste0(dir, "Figures and Data/pdf/AccuracyLinePlot_ClassFacetSD.pdf"), accuracy_lineplot_ClassFacet_SD, height = 10, width = 10)}

            if(filetag!=""){ggsave(paste0(dir,"Figures and Data/jpg/", filetag, "_AccuracyLinePlot_ClassFacetSD.jpg"), accuracy_lineplot_ClassFacet_SD, height = 10, width = 10)}else
              {ggsave(paste0(dir,"Figures and Data/jpg/AccuracyLinePlot_ClassFacetSD.jpg"),accuracy_lineplot_ClassFacet_SD, height = 10, width = 10)}

            # if(filetag!=""){write.csv(testsum, paste0(dir,"Figures and Data/data/", filetag,"_AccuracyLinePlot_ClassFacetSDData.csv"), row.names = FALSE, quote = FALSE)}else
            #   {write.csv(testsum, paste0(dir,"Figures and Data/data/AccuracyLinePlot_ClassFacetSDData.csv"), row.names = FALSE, quote = FALSE)}


          accuracy_lineplot_ClassFacet <-
            ggplot(SummaryAccuracy) +
            geom_line(aes(x = pofz, y = mean, colour = nLoci), lwd = 1.25) +
            facet_wrap(~max.class, nrow = 3) +
            theme(panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),
              legend.position="bottom", strip.background = element_rect(colour = "black", fill = "white")) +
            scale_color_brewer(palette = "Dark2")+
            labs(x = "Critical PofZ Threshold", y = expression("Proportion of Assignments Correct "%+-%"sd"), col="Panel Size (Loci)") +
            ylim(get.y.min.AccuracyLine, 1)

          if(filetag!=""){ggsave(paste0(dir,"Figures and Data/pdf/", filetag, "_AccuracyLinePlot_ClassFacet.pdf"), accuracy_lineplot_ClassFacet, height = 10, width = 10)}else
            {ggsave(paste0(dir, "Figures and Data/pdf/AccuracyLinePlot_ClassFacet.pdf"), accuracy_lineplot_ClassFacet, height = 10, width = 10)}

          if(filetag!=""){ggsave(paste0(dir,"Figures and Data/jpg/", filetag, "_AccuracyLinePlot_ClassFacet.jpg"), accuracy_lineplot_ClassFacet, height = 10, width = 10)}else
            {ggsave(paste0(dir,"Figures and Data/jpg/AccuracyLinePlot_ClassFacet.jpg"),accuracy_lineplot_ClassFacet, height = 10, width = 10)}

          # if(filetag!=""){write.csv(testsum, paste0(dir,"Figures and Data/data/", filetag,"_AccuracyLinePlot_ClassFacetData.csv"), row.names = FALSE, quote = FALSE)}else
          #   {write.csv(testsum, paste0(dir,"Figures and Data/data/AccuracyLinePlot_ClassFacetData.csv"), row.names = FALSE, quote = FALSE)}

          get.y.min.AccuracyThreshold <- min(dplyr::filter(SummaryAccuracy,pofz %in% Thresholds)$mean-dplyr::filter(SummaryAccuracy,pofz %in% Thresholds)$sd)

          Accuracy_ByThreshold_LinePlot_AllClass <-
            ggplot(dplyr::filter(SummaryAccuracy,pofz %in% Thresholds), aes(x=factor(nLoci),y=mean,col=max.class,group=max.class))+
            geom_point(size=2.5)+geom_path(lwd=0.9) +
            geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.1) +
            facet_grid(~pofz)+
            labs(x="Panel Size (Loci)",y=expression("Proportion of Assignments Correct "%+-%"sd"),col="Genotype Frequency Class",group="") +
            scale_color_brewer(palette = "Dark2")+
            theme(panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"),
              panel.grid.major = element_line(colour = "grey90"),legend.position="bottom",strip.background = element_rect(fill="white",colour = "black"), text = element_text(colour = "black")) +
            ylim(get.y.min.AccuracyThreshold, 1)


          if(filetag!=""){ggsave(paste0(dir,"Figures and Data/pdf/", filetag, "_Accuracy_ByThreshold_LinePlot_AllClass.pdf"), Accuracy_ByThreshold_LinePlot_AllClass, height = 10, width = 10)}else
            {ggsave(paste0(dir, "Figures and Data/pdf/Accuracy_ByThreshold_LinePlot_AllClass.pdf"), Accuracy_ByThreshold_LinePlot_AllClass, height = 10, width = 10)}

          if(filetag!=""){ggsave(paste0(dir,"Figures and Data/jpg/", filetag, "_Accuracy_ByThreshold_LinePlot_AllClass.jpg"), Accuracy_ByThreshold_LinePlot_AllClass, height = 10, width = 10)}else
            {ggsave(paste0(dir,"Figures and Data/jpg/Accuracy_ByThreshold_LinePlot_AllClass.jpg"),Accuracy_ByThreshold_LinePlot_AllClass, height = 10, width = 10)}

          # if(filetag!=""){write.csv(testsum, paste0(dir,"Figures and Data/data/", filetag,"_Accuracy_ByThreshold_LinePlot_AllClass.csv"), row.names = FALSE, quote = FALSE)}else
          #   {write.csv(testsum, paste0(dir,"Figures and Data/data/Accuracy_ByThreshold_LinePlot_AllClass.csv"), row.names = FALSE, quote = FALSE)}


          #ComboHybrids ------------
          ComboHybridAccuracy <- AccuracyData
          ComboHybridAccuracy$class <- as.character(AccuracyData$max.class)
          ComboHybridAccuracy$class[ComboHybridAccuracy$class %in% c("F1", "F2", "BC1", "BC2")] = "Hybrid"
          ComboHybridAccuracy$class[ComboHybridAccuracy$class == "P1"] = "Pure1"
          ComboHybridAccuracy$class[ComboHybridAccuracy$class == "P2"] = "Pure2"

          ComboHybridAccuracy <- ComboHybridAccuracy%>%
                            dplyr::group_by(nLoci,pofz,class)%>%
                            dplyr::summarise(mprob = mean(means,na.rm=T),
                                             sdprob = sd(means,na.rm=T))%>%
                            dplyr::ungroup()%>%data.frame()

          ComboHybridAccuracy$class <- factor(ComboHybridAccuracy$class, levels=c("Pure1","Pure2","Hybrid")) # set plotting levels

          Accuracy_ByThreshold_LinePlot_PureHyb <-
            ggplot(dplyr::filter(ComboHybridAccuracy, pofz %in% Thresholds), aes(x = factor(nLoci), y = mprob, col = class, group = class)) +
            geom_point(size=2.5) +
            geom_path(lwd=0.9) +
            geom_errorbar(aes(ymin = (mprob - sdprob),ymax = (mprob + sdprob)), width = 0.1) +
            facet_grid(~pofz) +
            labs(x = "Panel Size (Loci)", y = expression("Proportion of Assignments Correct "%+-%"sd"),
              col = "Genotype Frequency Class", group = "") +
            scale_color_brewer(palette = "Dark2") +
            theme(panel.background = element_rect(fill = "white", colour = "black"),
              plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),
              legend.position = "bottom", strip.background = element_rect(fill = "white", colour = "black"),
              text = element_text(colour = "black"))


          if(filetag!=""){ggsave(paste0(dir,"Figures and Data/pdf/", filetag, "_Accuracy_ByThreshold_LinePlot_PureHyb.pdf"), Accuracy_ByThreshold_LinePlot_PureHyb, height = 10, width = 10)}else
            {ggsave(paste0(dir, "Figures and Data/pdf/Accuracy_ByThreshold_LinePlot_PureHyb.pdf"), Accuracy_ByThreshold_LinePlot_PureHyb, height = 10, width = 10)}

          if(filetag!=""){ggsave(paste0(dir,"Figures and Data/jpg/", filetag, "_Accuracy_ByThreshold_LinePlot_PureHyb.jpg"), Accuracy_ByThreshold_LinePlot_PureHyb, height = 10, width = 10)}else
            {ggsave(paste0(dir,"Figures and Data/jpg/Accuracy_ByThreshold_LinePlot_PureHyb.jpg"),Accuracy_ByThreshold_LinePlot_PureHyb, height = 10, width = 10)}

          if(filetag!=""){write.csv(final.stats_pipe, paste0(dir,"Figures and Data/data/", filetag,"_Accuracy_ByThreshold_LinePlot_PureHyb.csv"), row.names = FALSE, quote = FALSE)}else
            {write.csv(dplyr::filter(ComboHybridAccuracy, pofz %in% Thresholds), paste0(dir,"Figures and Data/data/Accuracy_ByThreshold_LinePlot_PureHyb.csv"), row.names = FALSE, quote = FALSE)}


              ##################
              ## TYPE I ERROR ##
              ##################


              ### Calculate Type I Error - The number of individuals wrongly called hybrid over the actual number of purebreds in the sample

              writeLines("
                Calculating Type I Error
              ")

              ### TYPE I ERROR

              tempType1 <- temp
              tempType1$isgood <- TRUE ## set all isgood to "YES" - default is correct <- will switch to "NO" if PofZ value < critical value
              tempType1$domatch <- FALSE ## set all domatch to "NO" -  default is incorrect <- will swith to "YES" if they match

              #Create long form data for dplyr loop
              PofZVector_Type1 <- rep(50:99/100,each=nrow(tempType1)) #vector of PofZs

              tempType1Long <- do.call("rbind", replicate(length(50:99), tempType1, simplify = FALSE))
              tempType1Long$group <- as.character(PofZVector_Type1)
              tempType1Long$pofz <- PofZVector_Type1

              #dplyr loop
              TypeIout <- tempType1Long%>%
                    dplyr::group_by(group,nLoci)%>%dplyr::do(type1function(.))%>%
                    dplyr::ungroup()%>%data.frame()

              colnames(TypeIout) <- c("PofZ","Loci","Prop")

              TypeIout$Loci <- factor(x = TypeIout$Loci, levels = ordered(unique(as.numeric(as.character(TypeIout$Loci)))))

              TypeIout.PofZeds <- TypeIout
              TypeIout.PofZeds$Prop <- as.numeric(as.character(TypeIout.PofZeds$Prop))

              TypeI_BoxPlot <-
                ggplot(filter(TypeIout.PofZeds,PofZ %in% c("0.5","0.75","0.9")), aes(x = PofZ, y = Prop)) +
                geom_boxplot(fill = "grey75") +
                facet_grid(.~Loci) +
                theme(panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"),
                  panel.grid.major = element_line(colour = "grey90"),
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

              testsum_TypeI <- TypeIout%>%dplyr::group_by(PofZ,Loci)%>%
                dplyr::summarise(means=mean(Prop),sd=sd(Prop))%>%
                dplyr::ungroup()%>%data.frame()

              #set plot limits
              max.y = max(testsum_TypeI$means + testsum_TypeI$sd)

              ## line plot - Type I Error
              typeI_lineplot <-
              ggplot(testsum_TypeI) +
                geom_line(aes(x = as.numeric(PofZ), y = means, group = factor(Loci)), size = 2) +
                geom_line( aes(y = (means+sd), x = as.numeric(PofZ), group = factor(Loci)), lty = 2) +
                geom_line( aes(y = (means-sd), x = as.numeric(PofZ), group = factor(Loci)), lty = 2) +
                facet_wrap(~Loci, ncol = 3) +
                ylim(ymin = 0, ymax = max.y) +
                labs(x = "Critical PofZ Threshold", y = expression("Type I Error Proportion "%+-%"sd")) +
                theme(panel.background = element_rect(fill = "white", colour = "black"),
                  plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),
                  legend.position = "none", strip.background = element_rect(colour = "black", fill = "white"),
                  text = element_text(colour = "black"))

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
                sim_data <-dplyr::filter(output)%>%dplyr::group_by(nLoci,sim,Indv)%>%dplyr::summarise(Pure1_sd=sd(Pure1),Pure1=mean(Pure1),
                                                                       Pure2_sd=sd(Pure2),Pure2=mean(Pure2),
                                                                       F1_sd=sd(F1),F1=mean(F1),
                                                                       F2_sd=sd(F2),F2=mean(F2),
                                                                       BC1_sd=sd(BC1),BC1=mean(BC1),
                                                                       BC2_sd=sd(BC2),BC2=mean(BC2))%>%dplyr::ungroup()%>%data.frame()

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
                ggplot(boxdata, aes(x = nLoci, y = value, fill = sim)) +
                  geom_boxplot(alpha = 0.8, outlier.size = 0) +
                  facet_wrap(~class, nrow = 3, scales = "free_y") +
                  labs(y = "PofZ Score", x = "Panel Size (Loci)") +
                  scale_fill_manual(values = c("grey75", "grey75", "grey75")) +
                  theme(panel.background = element_rect(fill = "white", colour = "black"),
                    plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),
                    legend.position = "none", strip.background = element_rect(colour = "black", fill = "white"),
                    text = element_text(colour = "black"))

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
                  ggplot(sim_means2, aes(x = factor(nLoci), y = hybrid, fill = sim)) +
                  geom_boxplot(alpha = 0.8, outlier.size = 0) +
                  facet_wrap(~hclass, nrow = 3, scales = "free_y") +
                  labs(y = "PofZ Score", x = "Panel Size (Loci)") +
                  scale_fill_manual(values=c("grey75","grey75","grey75"))+
                  theme(panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),
                    legend.position = "none", strip.background = element_rect(colour = "black", fill = "white"), text = element_text(colour = "black"))

                #Save plot and data
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

                  #combined hybrid probabilities
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
                  FinalData <- data.frame(ProbOutput%>%dplyr::group_by(nLoci,level,class)%>%dplyr::summarise(mprob = mean(prob,na.rm=T),
                                                                            sdprob = sd(prob,na.rm=T))%>%dplyr::ungroup())
                  FinalData$class <- factor(FinalData$class, levels=c("Pure1","Pure2","F1","F2","BC1","BC2")) # NH class

                  # set plotting levels
                  FinalData$group <- "Pure"
                  FinalData[which(FinalData$class %in% c("BC1","BC2")),"group"] <- "Back-cross"
                  FinalData[which(FinalData$class %in% c("F1","F2")),"group"] <- "Generational hybrids"

                  FinalData$group <-  factor(FinalData$group,levels=c("Pure","Generational hybrids","Back-cross"))
                  FinalData$class <- factor(FinalData$class,levels=c("Pure1","Pure2","F1","F2","BC1","BC2"))

                  #plot the class groupings
                  Efficiency_AllClass_LinePlot <-
                  ggplot(FinalData, aes(x = level, y = mprob, col = class)) +
                  geom_line(lwd = 1.25) +
                  facet_grid( group ~ nLoci, scales = "free_y") +
                  theme(panel.background = element_rect(fill = "white", colour = "black"),
                    plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),
                    legend.position="bottom",strip.background = element_rect(fill="white",colour = "black"), text = element_text(colour = "black")) +
                  scale_color_brewer(palette = "Dark2") +
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
                    ggplot(FinalData2) +
                    geom_line(aes(x=level,y=mprob,col=class),lwd=1.25)+
                    geom_line(aes(x=level,y=mprob+sdprob,col=class),lty=2)+
                    geom_line(aes(x=level,y=mprob-sdprob,col=class),lty=2)+
                    facet_wrap(~nLoci, ncol = 3)+
                    theme(panel.background = element_rect(fill = "white", colour = "black"),
                      plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),
                      legend.position="bottom",strip.background = element_rect(fill="white",colour = "black"), text = element_text(colour = "black")) +
                    scale_color_brewer(palette = "Dark2")+
                    labs(x="Critical PofZ Threshold",y=expression("Efficiency "%+-%"sd"),col="Genotype Frequency Class") +
                    ylim(0, 1)


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
                      facet_wrap(~class,nrow=3,scales="free_y")+
                      scale_color_brewer(palette = "Dark2")+
                      theme(panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"),
                        panel.grid.major = element_line(colour = "grey90"), legend.position = "bottom", strip.background = element_rect(fill = "white", colour = "black"), text = element_text(colour = "black"))+
                      labs(x="Critical PofZ Threshold",y=expression("Efficiency "%+-%"sd"),col="Panel Size (Loci)")+
                      geom_vline(xintercept = Thresholds, lty=2) +
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
                      facet_wrap(~class,nrow=3,scales="free_y")+
                      scale_color_brewer(palette = "Dark2")+
                      theme(panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"),
                        panel.grid.major = element_line(colour = "grey90"),legend.position="bottom",strip.background = element_rect(fill="white",colour = "black"), text = element_text(colour = "black"))+
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
                      facet_grid(~class)+
                      scale_color_brewer(palette = "Dark2")+
                      theme(panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"),
                        panel.grid.major = element_line(colour = "grey90"),legend.position="bottom",strip.background = element_rect(fill="white",colour = "black"), text = element_text(colour = "black")) +
                      labs(x="Critical PofZ Threshold",y=expression("Efficiency "%+-%"sd"),col="Panel Size (Loci)")+
                      ylim(0, 1) +
                      geom_vline(xintercept = Thresholds, lty=2)

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
                      facet_grid(~class)+
                      scale_color_brewer(palette = "Dark2")+
                      theme(panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"),
                        panel.grid.major = element_line(colour = "grey90"),legend.position="bottom",strip.background = element_rect(fill="white",colour = "black"), text = element_text(colour = "black")) +
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
                    facet_grid(~level)+
                    labs(x="Panel Size (Loci)",y=expression("Efficiency "%+-%"sd"),col="Genotype Frequency Class",group="")+scale_color_brewer(palette = "Dark2")+
                    theme(panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"),
                      panel.grid.major = element_line(colour = "grey90"),legend.position="bottom",strip.background = element_rect(fill="white",colour = "black"), text = element_text(colour = "black"))


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
                    facet_grid(~level)+
                    labs(x="Panel Size (Loci)",y=expression("Efficiency "%+-%"sd"),col="Genotype Frequency Class",group="")+
                    scale_color_brewer(palette = "Dark2")+
                    theme(panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"),
                      panel.grid.major = element_line(colour = "grey90"),legend.position="bottom",strip.background = element_rect(fill="white",colour = "black"), text = element_text(colour = "black"))

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

                    classnames <- c("Pure1","Pure2","F1","F2","BC1","BC2")

                    tempType2 <- sim_means

                    #Create long form data for dplyr loop
                    PofZVector_Type2 <- rep(50:99/100,each=nrow(tempType2)) #vector of PofZs

                    tempType2Long <- do.call("rbind", replicate(length(50:99), tempType2, simplify = FALSE))
                    tempType2Long$group <- as.character(PofZVector_Type2)
                    tempType2Long$pofz <- PofZVector_Type2
                    tempType2Long$samplesize=paste(samplesize,collapse=",") #wild card variables to be incorperated into the do function directly
                    tempType2Long$classnames=paste(classnames,collapse=",")

                    missout <- tempType2Long%>%dplyr::group_by(nLoci,sim,pofz)%>%
                        dplyr::do(type2function(.))%>%
                        dplyr::ungroup()%>%data.frame()


                    #calcluate the means among simulations
                  miss_mean <- missout%>%dplyr::group_by(nLoci,pofz,class)%>%
                      dplyr::summarise(mprobP1 = mean(mclass_P1,na.rm=T),
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
                      sdprobBC2 = sd(mclass_BC2,na.rm=T))%>%
                      dplyr::ungroup()%>%data.frame()

                    miss_mean[is.na(miss_mean)]=NA #replace NaN's with NAs

                    colnames(miss_mean)[grep("pofz",colnames(miss_mean))] = "level"

                    #merge with the other data
                    FinalData3 <- merge(miss_mean,FinalData,by=c("nLoci","level","class"))

                    PlotData <- reshape2::melt(FinalData3[c("nLoci","level","class","mprobP1","mprobP2","mprobF1","mprobF2","mprobBC1","mprobBC2")],
                      id.vars=c("nLoci","level","class"))

                    PlotDatasd <- reshape2::melt(FinalData3[c("nLoci","level","class","sdprobP1","sdprobP2","sdprobF1","sdprobF2","sdprobBC1","sdprobBC2")],
                       id.vars=c("nLoci","level","class"))

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
                      temp.plot <- ggplot(dplyr::filter(PlotData,class==i))+
                      geom_line(aes(x=level,y=value,col=variable),lwd=1.25)+
                      geom_line(aes(x=level,y=value+sd,col=variable),lty=2)+
                      geom_line(aes(x=level,y=value-sd,col=variable),lty=2)+
                      facet_grid(~nLoci,scales="free_y")+
                      scale_color_brewer(palette = "Dark2")+
                      theme(legend.position="bottom",strip.background = element_rect(fill="white",colour = "black"))+
                      labs(x="Critical PofZ Threshold",y=paste0("Proportion ",i," misassigned \u00B1 sd"),col="Genotype Frequency Class") +
                        expand_limits(y = 0)

                      # temp.plot

                      if(filetag!=""){ggsave(paste0(dir,"Figures and Data/pdf/",filetag,"_",i,"_TypeII_nloci.pdf"),temp.plot,height = 6,width = 8)} else
                        {ggsave(paste0(dir,paste0("Figures and Data/pdf/",i,"_TypeII_nloci.pdf")),temp.plot,height = 6,width = 8)}

                      if(filetag!=""){ggsave(paste0(dir,"Figures and Data/jpg/",filetag,"_",i,"_TypeII_nloci.jpg"),temp.plot,height = 6,width = 8)} else
                      {ggsave(paste0(dir,paste0("Figures and Data/jpg/",i,"_TypeII_nloci.jpg")),temp.plot,height = 6,width = 8)}

                      }


    #### CHECK HERE - this doesn't work yet  -- fix the plot names.
        ## the plot names are:
#     accuracy_boxplot
# accuracy_lineplotSD
# accuracy_lineplot
# accuracy_lineplot_ClassFacet_SD
# accuracy_lineplot_ClassFacet
# Accuracy_ByThreshold_LinePlot_AllClass
# Accuracy_ByThreshold_LinePlot_PureHyb
# TypeI_BoxPlot
# typeI_lineplot
# meanPofZ_AllClasses_BoxPlot
# meanPofZ_PureHyb_BoxPlot
# Efficiency_AllClass_LinePlot
# Efficiency_PureHyb_LinePlot
# Efficiency_AllClass_LinePlot_ClassFacet_thresh
# Efficiency_AllClass_LinePlot_ClassFacet_NOthresh
# Efficiency_PureHyb_LinePlot_ClassFacet_thresh
# Efficiency_PureHyb_LinePlot_ClassFacet_NOthresh
# Efficiency_ByThreshold_LinePlot_AllClass
# Efficiency_ByThreshold_LinePlot_PureHyb



    ## Create summary booklet

    # if(filetag!=""){
    #   pdf(file = paste0(dir,"Figures and Data/pdf/",filetag,"_OutputBooklet.pdf"))
    #   print(p1);print(h1)
    #   print(p3);print(h3)
    #   print(p4);print(h4)
    #   print(p5);print(h5)
    #   for (i in unique(PlotData$class)){
    #     temp.plot <- ggplot(filter(PlotData,class==i))+
    #       geom_line(aes(x=level,y=value,col=variable),lwd=1.25)+
    #       geom_line(aes(x=level,y=value+sd,col=variable),lty=2)+
    #       geom_line(aes(x=level,y=value-sd,col=variable),lty=2)+
    #       facet_grid(~nLoci,scales="free_y")+
    #       theme_bw()+scale_color_brewer(palette = "Dark2")+
    #       theme(legend.position="bottom",strip.background = element_rect(fill="white",colour = "black"))+
    #       labs(x="Probability threshold",y=paste0("Proportion ",i," misassigned \u00B1 sd"),col="Classification")
    #       print(temp.plot)
    #       }
    #   dev.off()
    # } else {
    #   pdf(file = paste0(dir,"Figures and Data/pdf/",filetag,"_OutputBooklet.pdf"))
    #   print(p1);print(h1)
    #   print(p3);print(h3)
    #   print(p4);print(h4)
    #   print(p5);print(h5)
    #   for (i in unique(PlotData$class)){
    #     temp.plot <- ggplot(filter(PlotData,class==i))+
    #       geom_line(aes(x=level,y=value,col=variable),lwd=1.25)+
    #       geom_line(aes(x=level,y=value+sd,col=variable),lty=2)+
    #       geom_line(aes(x=level,y=value-sd,col=variable),lty=2)+
    #       facet_grid(~nLoci,scales="free_y")+
    #       theme_bw()+scale_color_brewer(palette = "Dark2")+
    #       theme(legend.position="bottom",strip.background = element_rect(fill="white",colour = "black"))+
    #       labs(x="Probability threshold",y=paste0("Proportion ",i," misassigned \u00B1 sd"),col="Classification")
    #     print(temp.plot)
    #   }
    #   dev.off()
    # }

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
