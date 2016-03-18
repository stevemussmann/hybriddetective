#' @name hybridpower
#' @title Power evaluation for NewHybrids analysis of simulated datasets
#' @description Evaluates the accuracy with which NewHybrids assigns individuals of known hybrid class to the correct hybrid class in simulated datasets at varying levels of stringency (PofZ). The code will write graphical and numerical results to the directory provided by the user.
#' @param dir path directory which holds the output from different runs through New Hybrids (e.g. 3 simulations with 3 replicate runs each through NH) note that this directory should only hold the output folders.
#' @param filetag A name tag which will be added to the outputs
#' @param Threshold A threshold which will be added to the plots showing the assignment success for different levels of probability of a given class estimated by NewHybrids. Default is (NULL) so if nothing is specified it will not add this to the output plots (success ~ threshold by class)
#' @param Thresholds A vector of thresholds which will be added to the plots showing the assignment success for different levels of probability of a given class estimated by NewHybrids. Default is (NULL) so if nothing is specified it will not add this to the output plots (success ~ threshold by class).
#' @param samplesize is the number of fish per NH class. By (default: NULL) this data will be extracted from the "*individuals.txt" output from parallelnewhybrids. This can also explicitly defined as a vector (6 values corresponding to # in P1,P2,F1,F2,BC1,BC2) or a path to the *_Individuals.txt output from.
#' @param CT convergence threshold (default: 0.1) denoting what an acceptable proportion of each individual of P1 & P2  can be classified as "F2".
#' @param CTI proportion of individuals (default: 0.5) within a class (P1 and P2) which are permitted to fail exceed CT.
#' @rdname hybridpower
#' @import ggplot2
#' @import magrittr
#' @importFrom dplyr filter summarise ungroup group_by
#' @importFrom grid arrow unit
#' @importFrom stringr str_extract
#' @importFrom reshape2 melt
#' @importFrom scales alpha
#' @export



hybridpower <-function(dir,filetag="",Threshold=NULL,Thresholds=c(0.5,0.6,0.7,0.8,0.9),samplesize=NULL,CT=0.1,CTI=0.5) {

  #set directory for which holds the New Hybrids output folders
  filedir <- dir
  lfiles <- setdiff(list.files(dir),"Figures and Data") #ignores Figures folder in case this is run more than once
  if(length(which(list.files(dir)=="Figures and Data"))==0){dir.create(paste0(dir,"Figures and Data"))} # if there isn't a 'Figures and Data' folder for output create one
  if(length(which(list.files(paste0(dir,"Figures and Data"))=="pdf"))==0){dir.create(paste0(dir,"Figures and Data/pdf"))} #create a folder for pdfs
  if(length(which(list.files(paste0(dir,"Figures and Data"))=="jpg"))==0){dir.create(paste0(dir,"Figures and Data/jpg"))} #create a folder for jpgs
  if(length(which(list.files(paste0(dir,"Figures and Data"))=="data"))==0){dir.create(paste0(dir,"Figures and Data/data"))} #create a folder for data

  #Convergence checker
  arethereproblems = "No"

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


    ## average and SD the  replicate runs of each simulation in New Hybrids
    sim_data <- as.data.frame(output%>%group_by(sim,Indv)%>%
                              summarise(Pure1_sd=sd(Pure1),Pure1=mean(Pure1),
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
    for (i in unique(sim_means$class))
    {
      temp <- filter(sim_means,class==i)
      tout <- data.frame(sim=temp$sim,Indv=temp$Indv,class=i,value=temp[,i])
      boxdata <- rbind(boxdata,tout)
    }

    # Create pot
     p1 <- ggplot(boxdata,aes(x=sim,y=value))+geom_boxplot(alpha=0.8,outlier.size = 0)+theme_bw()+facet_wrap(~class,nrow=3,scales="free_y")+
      labs(y="Class probability",x="Number of loci")+
      theme(strip.background = element_rect(fill="white"),legend.position="none")

    #save plot + Data
    if(filetag!=""){ggsave(paste0(dir,"Figures and Data/pdf/",filetag,"_AssignmentSuccess~simulation-nSNPs_p1.pdf"),p1,height = 8,width = 10)}else
    {ggsave(paste0(dir,"Figures and Data/pdf/AssignmentSuccess~simulation-nSNPs_p1.pdf"),p1,height = 8,width = 10)}

    if(filetag!=""){ggsave(paste0(dir,"Figures and Data/jpg/",filetag,"_AssignmentSuccess~simulation-nSNPs_p1.jpg"),p1,height = 8,width = 10)}else
    {ggsave(paste0(dir,"Figures and Data/jpg/AssignmentSuccess~simulation-nSNPs_p1.jpg"),p1,height = 8,width = 10)}

    if(filetag!=""){write.csv(boxdata, paste0(dir,"Figures and Data/data/",filetag,"_AssignmentSuccess~simulation-nSNPs_OUTPUT.csv"))}else
    {write.csv(boxdata, paste0(dir,"Figures and Data/data/AssignmentSuccess~simulation-nSNPs_OUTPUT.csv"))}

    #Make similar boxplot with summed propbability of being a hybrid
    sim_means2 <- sim_means
    sim_means2$hybrid <- rowSums(sim_means2[,c("F1","F2","BC1","BC2")])
    sim_means2[which(sim_means2$class =="Pure1"),"hybrid"]=sim_means2[which(sim_means2$class =="Pure1"),"Pure1"] #add values of the Pure
    sim_means2[which(sim_means2$class =="Pure2"),"hybrid"]=sim_means2[which(sim_means2$class =="Pure2"),"Pure2"] #add values of the Pure

    sim_means2$hclass <- "Hybrid"
    sim_means2[which(sim_means$class=="Pure1"),"hclass"] <- "Pure1"
    sim_means2[which(sim_means$class=="Pure2"),"hclass"] <- "Pure2"

    sim_means2$hclass <- factor(sim_means2$hclass,levels=c("Pure1","Pure2","Hybrid"))

    h1 <- ggplot(sim_means2,aes(x=sim,y=hybrid))+geom_boxplot(alpha=0.8,outlier.size = 0)+theme_bw()+facet_wrap(~hclass,nrow=3)+
      labs(y="Class probability",x="Simulation")+theme(strip.background = element_rect(fill="white"))+scale_y_continuous(limits=c(0.9,1))

     #Save plot
    if(filetag!=""){ggsave(paste0(dir,"Figures and Data/pdf/",filetag,"_AssignmentSuccess~simulation-nSNPs_Hybrid_h1.pdf"),h1,height = 8,width = 8)}else
    {ggsave(paste0(dir,"Figures and Data/pdf/AssignmentSuccess~simulation-nSNPs_Hybrid_h1.pdf"),h1,height = 8,width = 8)}

    if(filetag!=""){ggsave(paste0(dir,"Figures and Data/jpg/",filetag,"_AssignmentSuccess~simulation-nSNPs_Hybrid_h1.jpg"),h1,height = 8,width = 8)}else
    {ggsave(paste0(dir,"Figures and Data/jpg/AssignmentSuccess~simulation-nSNPs_Hybrid_h1.jpg"),h1,height = 8,width = 8)}

    if(filetag!=""){write.csv(sim_means2, paste0(dir,"Figures and Data/data/",filetag,"_AssignmentSuccess~simulation-nSNPs_Hybrid_OUTPUT.csv"))}else
    {write.csv(sim_means2, paste0(dir,"Figures and Data/data/AssignmentSuccess~simulation-nSNPs_Hybrid_OUTPUT.csv"))}

## create the New Hybrids plot

      NH_melt <- melt(data = sim_means[,-grep("class",colnames(sim_means))], id.vars = c("Indv","sim")) ## melt the data to allow the data to be stacked by indivudal
      colnames(NH_melt) <- c("Indv","sim", "PopProb", "CumProb") ## rename so that its prettier

      NH_melt$sim <- gsub("S","Simulation ",NH_melt$sim) # fix the simulation labels

      ## Plot colours
      col.vec <- c("red", "blue", "grey", "green", "black", "yellow", "brown")


      ## make a nice pretty little plot - and save that bad boy
      #f.name <- paste0(res.folder, res.name, ".consensus.pdf")
      #pdf(file = f.name)

      p2 <- ggplot(NH_melt, aes(x = Indv, y=CumProb, fill=PopProb))+
      geom_bar(stat="identity", position = "stack") +
        scale_fill_manual(values=col.vec)+
        labs(y="Cumulative Probability",x="Individual",fill="")+
        scale_y_continuous(limits = c(0, 1), expand=c(0, 0)) +
        scale_x_continuous(expand=c(0,0)) +
        facet_grid(sim~.)+
        theme(axis.text.x = element_text(colour = "black"), axis.title.x = element_blank(),
              axis.text.y = element_text(colour = "black"), axis.title.y = element_text(size = 15, colour = "black"),
              strip.background = element_rect(fill="white",colour = "black"),panel.margin=unit(1, "lines"))

      #save the figure
      if(filetag!=""){ggsave(paste0(dir,"Figures and Data/pdf/",filetag,"_NewHybridsPlot~simulation_p2.pdf"),p2,height = 8,width = 10)} else
      {ggsave(paste0(dir,"Figures and Data/pdf/NewHybridsPlot~simulation_p2.pdf"),p2,height = 8,width = 10)}


## Look at assignment success as a function of threshold probability
classvec2 <- rep(c("Pure1","Pure2","F1","F2","BC1","BC2"),times=samplesize)

    ProbOutput <- NULL
    for(i in unique(sim_means$sim)){
      tempsub <- filter(sim_means,sim==i)
        for(q in 50:99/100){ # probability of 50 - 99%

          p1.p <- length(which(tempsub[which(classvec2=="Pure1"),"Pure1"] > q))/samplesize[1]
          p2.p <- length(which(tempsub[which(classvec2=="Pure2"),"Pure2"] > q))/samplesize[2]
          F1.p <- length(which(tempsub[which(classvec2=="F1"),"F1"] > q))/samplesize[3]
          F2.p <- length(which(tempsub[which(classvec2=="F2"),"F2"] > q))/samplesize[4]
          BC1.p <- length(which(tempsub[which(classvec2=="BC1"),"BC1"] > q))/samplesize[5]
          BC2.p <- length(which(tempsub[which(classvec2=="BC2"),"BC2"] > q))/samplesize[6]
          tempout <- data.frame(sim=i,level=q,prob=c(p1.p, p2.p,F1.p,F2.p,BC1.p,BC2.p),
                                class=c("Pure1","Pure2","F1","F2","BC1","BC2"))
          ProbOutput <- rbind(ProbOutput,tempout)

        }
    }

    #combined hybrid probabilities
    ProbOutput2 <- NULL
    for(i in unique(sim_means$sim)){
      tempsub <- filter(sim_means,sim==i)
      tempsub$phyb <- rowSums(tempsub[,c("F1","F2","BC1","BC2")])
      for(q in 50:99/100){ # probability of 50 - 99%

        p1.p <- length(which(tempsub[which(classvec2=="Pure1"),"Pure1"] > q))/samplesize[1]
        p2.p <- length(which(tempsub[which(classvec2=="Pure2"),"Pure2"] > q))/samplesize[2]
        Hybrid <- length(which(tempsub[which(classvec2%in%c("F1","F2","BC1","BC2")),"phyb"]>q))/sum(samplesize[3:6])
        tempout <- data.frame(sim=i,level=q,prob=c(p1.p, p2.p,Hybrid),
                              class=c("Pure1","Pure2","Hybrid"))
        ProbOutput2 <- rbind(ProbOutput2,tempout)

      }
    }

# get the mean and standard error for the estimates of assignment succes based on NH probabilty among simulations
      FinalData <- data.frame(ProbOutput%>%group_by(level,class)%>%summarise(mprob = mean(prob,na.rm=T),
                                                                            sdprob = sd(prob,na.rm=T))%>%ungroup())
      FinalData$class <- factor(FinalData$class, levels=c("Pure1","Pure2","F1","F2","BC1","BC2")) # NH class

    # set plotting levels
      FinalData$group <- "Pure"
      FinalData[which(FinalData$class %in% c("BC1","BC2")),"group"] <- "Back-cross"
      FinalData[which(FinalData$class %in% c("F1","F2")),"group"] <- "Generational hybrids"
      FinalData$group <-  factor(FinalData$group,levels=c("Pure","Generational hybrids","Back-cross"))

      #plot the class groupings
      p3 <- ggplot(FinalData,aes(x=level,y=mprob,col=class))+geom_line(lwd=1.25)+theme_bw()+facet_grid(group~.,scales="free_y")+
        theme(legend.position="bottom",strip.background = element_rect(fill="white",colour = "black"))+
        scale_color_brewer(palette = "Dark2")+
        labs(x="Probability threshold",y="Assignment success",col="Classification")

      if(filetag!=""){ggsave(paste0(dir,"Figures and Data/pdf/",filetag,"_AssinmentSuccess~level-class_p3.pdf"),p3,height = 10,width = 8)} else
      {ggsave(paste0(dir,"Figures and Data/pdf/AssinmentSuccess~level-class_p3.pdf"),p3,height = 10,width = 8)}

      if(filetag!=""){ggsave(paste0(dir,"Figures and Data/jpg/",filetag,"_AssinmentSuccess~level-class_p3.jpg"),p3,height = 10,width = 8)} else
      {ggsave(paste0(dir,"Figures and Data/jpg/AssinmentSuccess~level-class_p3.jpg"),p3,height = 10,width = 8)}

      #ComboHybrids
      FinalData2 <- data.frame(ProbOutput2%>%group_by(level,class)%>%summarise(mprob = mean(prob,na.rm=T),
                                                                               sdprob = sd(prob,na.rm=T))%>%ungroup())
      FinalData2$class <- factor(FinalData2$class, levels=c("Pure1","Pure2","Hybrid")) # set plotting levels

      h3 <- ggplot(FinalData2)+
        geom_line(aes(x=level,y=mprob,col=class),lwd=1.25)+
        geom_line(aes(x=level,y=mprob+sdprob,col=class),lty=2)+
        geom_line(aes(x=level,y=mprob-sdprob,col=class),lty=2)+
        theme_bw()+
        theme(legend.position="bottom",strip.background = element_rect(fill="white",colour = "black"))+
        scale_color_brewer(palette = "Dark2")+
        labs(x="Probability threshold",y=expression("Assignment success "%+-%"sd"),col="Classification")

      if(filetag!=""){ggsave(paste0(dir,"Figures and Data/pdf/",filetag,"_AssinmentSuccess~level-class_Hybrid_h3.pdf"),h3,height = 8,width = 8)} else
      {ggsave(paste0(dir,"Figures and Data/pdf/AssinmentSuccess~level-class_Hybrid_h3.pdf"),h3,height = 8,width = 8)}

      if(filetag!=""){ggsave(paste0(dir,"Figures and Data/jpg/",filetag,"_AssinmentSuccess~level-class_Hybrid_h3.jpg"),h3,height = 8,width = 8)} else
      {ggsave(paste0(dir,"Figures and Data/jpg/AssinmentSuccess~level-class_Hybrid_h3.jpg"),h3,height = 8,width = 8)}

    #plot if no threshold specified
      if(length(Threshold)==0){
      p4 <- ggplot(data=FinalData)+geom_line(aes(x=level,y=mprob),lwd=1.25)+
        geom_line(aes(x=level,y=mprob+sdprob),lty=2)+geom_line(aes(x=level,y=mprob-sdprob),lty=2)+
        theme_bw()+facet_wrap(~class,nrow=3,scales="free_y")+theme(strip.background = element_rect(fill="white",colour = "black"))+
        labs(x="Probability threshold",y=expression("Assignment success "%+-%"sd"));p4
      }

    #plot if threshold specified
       if(length(Threshold)!=0){
        #calcualte the lower assignment success for each class
        lowerlim <- data.frame(FinalData%>%group_by(class)%>%
                                 summarise(lim=min(mprob))%>%ungroup())

        threshlimits <- filter(FinalData,level==Threshold)
        threshlimits$lower <- 0.5
        threshlimits$threshold <- Threshold
        threshlimits$x <- (threshlimits$threshold-0.5)/2
        threshlimits <- merge(threshlimits,lowerlim,by = "class")
        threshlimits$lim <- threshlimits$lim-threshlimits$sdprob

        p4 <- ggplot(data=FinalData)+
          geom_line(aes(x=level,y=mprob),lwd=1.25)+
          geom_line(aes(x=level,y=mprob+sdprob),lty=2)+
          geom_line(aes(x=level,y=mprob-sdprob),lty=2)+
          theme_bw()+facet_wrap(~class,nrow=3,scales="free_y")+
          geom_segment(data=threshlimits,aes(x=lower,xend=level,y=mprob,yend=mprob),col="black")+
          geom_segment(data=threshlimits,aes(x=level,xend=level,y=lim,yend=mprob),col="black")+
          scale_y_continuous(expand=c(0.03,0))+scale_x_continuous(expand=c(0.04,0))+
          labs(x="Probability threshold",y=expression("Assignment success "%+-%"sd"))+
          theme(strip.background = element_rect(fill="white",colour = "black"))+
          geom_text(data=threshlimits,aes(x=0.55,y = mprob-(mprob-lim)*0.1,label=round(mprob,2)));p4
       }

      #Save plot
      if(filetag!=""){ggsave(paste0(dir,"Figures and Data/pdf/",filetag,"_AssignmentSuccess~level-error_p4.pdf"),p4,height = 10,width = 8)} else
      {ggsave(paste0(dir,"Figures and Data/pdf/AssignmentSuccess~level-error_p4.pdf"),p4,height = 10,width = 8)}

      if(filetag!=""){ggsave(paste0(dir,"Figures and Data/jpg/",filetag,"_AssignmentSuccess~level-error_p4.jpg"),p4,height = 10,width = 8)} else
      {ggsave(paste0(dir,"Figures and Data/jpg/AssignmentSuccess~level-error_p4.jpg"),p4,height = 10,width = 8)}

      ### Combo Hybrids

      #plot if no threshold specified
      if(length(Threshold)==0){
        h4 <- ggplot(data=FinalData2)+geom_line(aes(x=level,y=mprob),lwd=1.25)+
          geom_line(aes(x=level,y=mprob+sdprob),lty=2)+geom_line(aes(x=level,y=mprob-sdprob),lty=2)+
          theme_bw()+facet_wrap(~class,nrow=3,scales="free_y")+theme(strip.background = element_rect(fill="white",colour = "black"))+
          labs(x="Probability threshold",y="Assignment success ± sd");h4
      }

      #plot if threshold specified
      if(length(Threshold)!=0){
        #calcualte the lower assignment success for each class
        lowerlim <- data.frame(FinalData2%>%group_by(class)%>%
                                 summarise(lim=min(mprob))%>%ungroup())

        threshlimits <- filter(FinalData2,level==Threshold)
        threshlimits$lower <- 0.5
        threshlimits$threshold <- Threshold
        threshlimits$x <- (threshlimits$threshold-0.5)/2
        threshlimits <- merge(threshlimits,lowerlim,by = "class")
        threshlimits$lim <- threshlimits$lim-threshlimits$sdprob

        h4 <- ggplot(data=FinalData2)+
          geom_line(aes(x=level,y=mprob),lwd=1.25)+
          geom_line(aes(x=level,y=mprob+sdprob),lty=2)+
          geom_line(aes(x=level,y=mprob-sdprob),lty=2)+
          theme_bw()+facet_wrap(~class,nrow=3,scales="free_y")+
          geom_segment(data=threshlimits,aes(x=lower,xend=level,y=mprob,yend=mprob),col="black")+
          geom_segment(data=threshlimits,aes(x=level,xend=level,y=lim,yend=mprob),col="black")+
          scale_y_continuous(expand=c(0.03,0))+scale_x_continuous(expand=c(0,0))+
          labs(x="Probability threshold",y="Assignment success ± sd")+
          theme(strip.background = element_rect(fill="white",colour = "black"))+
          geom_text(data=threshlimits,aes(x=0.55,y = mprob-(mprob-lim)*0.1,label=round(mprob,3)));h4
      }

      #Save plot
      if(filetag!=""){ggsave(paste0(dir,"Figures and Data/pdf/",filetag,"_AssignmentSuccess~level-error_Hybrid_h4.pdf"),h4,height = 10,width = 8)} else
      {ggsave(paste0(dir,"Figures and Data/pdf/AssignmentSuccess~level-error_Hybrid_h4.pdf"),h4,height = 10,width = 8)}

      if(filetag!=""){ggsave(paste0(dir,"Figures and Data/jpg/",filetag,"_AssignmentSuccess~level-error_Hybrid_h4.jpg"),h4,height = 10,width = 8)} else
      {ggsave(paste0(dir,"Figures and Data/jpg/AssignmentSuccess~level-error_Hybrid_h4.jgp"),h4,height = 10,width = 8)}

      ## mean plot ----------

      #facet labels
      FinalData$threshold <- paste0(FinalData$level*100,"%")

      p5 <- ggplot(filter(FinalData,level %in% Thresholds),aes(x=threshold,y=mprob,col=class,group=class))+
        geom_point(size=2.5)+geom_path(lwd=0.9)+
        geom_errorbar(aes(ymin=mprob-sdprob,ymax=mprob+sdprob),width=0.1)+
        theme_bw()+
        labs(x="Probability threshold",y=expression("Assignment success "%+-%"sd"),col="Classification",group="")+scale_color_brewer(palette = "Dark2")+
        theme(legend.position="bottom",strip.background = element_rect(fill="white",colour = "black"))

      #Save plot
      if(filetag!=""){ggsave(paste0(dir,"Figures and Data/pdf/",filetag,"_AssignmentSuccess~z-loci_p5.pdf"),p5,height = 8,width = 10)} else
      {ggsave(paste0(dir,"Figures and Data/pdf/AssignmentSuccess~~z-loci_p5.pdf"),p5,height = 8,width = 10)}

      if(filetag!=""){ggsave(paste0(dir,"Figures and Data/jpg/",filetag,"_AssignmentSuccess~z-loci_p5.jpg"),p5,height = 8,width = 10)} else
      {ggsave(paste0(dir,"Figures and Data/jpg/AssignmentSuccess~~z-loci_p5.jpg"),p5,height = 8,width = 10)}


      #Combined Hybrids
      FinalData2$threshold <- paste0(FinalData2$level*100,"%")

      h5 <- ggplot(filter(FinalData2,level %in% Thresholds),aes(x=threshold,y=mprob,col=class,group=class))+
        geom_point(size=2.5)+geom_path(lwd=0.9)+
        geom_errorbar(aes(ymin=mprob-sdprob,ymax=mprob+sdprob),width=0.1)+
        theme_bw()+
        labs(x="Probability threshold",y=expression("Assignment success "%+-%"sd"),col="Classification",group="")+scale_color_brewer(palette = "Dark2")+
        theme(legend.position="bottom",strip.background = element_rect(fill="white",colour = "black"))

      #Save plot
      if(filetag!=""){ggsave(paste0(dir,"Figures and Data/pdf/",filetag,"_AssignmentSuccess~z-loci_Hybrid_h5.pdf"),h5,height = 8,width = 10)} else
      {ggsave(paste0(dir,"Figures and Data/pdf/AssignmentSuccess~~z-loci_Hybrid_h5.pdf"),h5,height = 9,width = 10)}

      if(filetag!=""){ggsave(paste0(dir,"Figures and Data/jpg/",filetag,"_AssignmentSuccess~z-loci_Hybrid_h5.jpg"),h5,height = 8,width = 10)} else
      {ggsave(paste0(dir,"Figures and Data/jpg/AssignmentSuccess~~z-loci_Hybrid_h5.jpg"),h5,height = 8,width = 10)}


      ## Misclassification 'type II' error ------------
      classnames <- c("Pure1","Pure2","F1","F2","BC1","BC2")
      missout <- NULL

        for(i in unique(sim_means$sim)){
          tempsub <- filter(sim_means,sim==i)
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

                tempout <- data.frame(sim=i,level=q,class=z,
                                      mclass_P1=temp4[which(temp4$Var1 == "Pure1"),"Freq"],
                                      mclass_P2=temp4[which(temp4$Var1 == "Pure2"),"Freq"],
                                      mclass_F1=temp4[which(temp4$Var1 == "F1"),"Freq"],
                                      mclass_F2=temp4[which(temp4$Var1 == "F2"),"Freq"],
                                      mclass_BC1=temp4[which(temp4$Var1 == "BC1"),"Freq"],
                                      mclass_BC2=temp4[which(temp4$Var1 == "BC2"),"Freq"])
              } else
              {tempout <- data.frame(sim=i,level=q,class=z,
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

      #calcluate the means among simulations
      miss_mean <- data.frame(missout%>%group_by(level,class)
                              %>%summarise(mprobP1 = mean(mclass_P1,na.rm=T),
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

    #replace NaN's with NAs
    miss_mean[is.na(miss_mean)]=NA



    #merge with the other data
    FinalData3 <- merge(miss_mean,FinalData,by=c("level","class"))

    PlotData <- melt(FinalData3[c("level","class","mprobP1","mprobP2","mprobF1","mprobF2","mprobBC1","mprobBC2")],id.vars=c("level","class"))
    PlotDatasd <- melt(FinalData3[c("level","class","sdprobP1","sdprobP2","sdprobF1","sdprobF2","sdprobBC1","sdprobBC2")],id.vars=c("level","class"))
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
      plotlab=paste("Proportion",i,"missassigned",sep=" ")
      temp.plot <- ggplot(filter(PlotData,class==i))+
        geom_line(aes(x=level,y=value,col=variable),lwd=1.25)+
        geom_line(aes(x=level,y=value+sd,col=variable),lty=2)+
        geom_line(aes(x=level,y=value-sd,col=variable),lty=2)+
        theme_bw()+scale_color_brewer(palette = "Dark2")+
        theme(legend.position="bottom",strip.background = element_rect(fill="white",colour = "black"))+
        labs(x="Probability threshold",y=paste0("Proportion ",i," misassigned \u00B1 sd"),col="Classification")

      if(filetag!=""){ggsave(paste0(dir,"Figures and Data/pdf/",filetag,"_",i,"_MissAssignment~z-nloci.pdf"),temp.plot,height = 6,width = 8)} else
      {ggsave(paste0(dir,paste0("Figures and Data/pdf/",i,"_MissAssignment~z-nloci.pdf"),temp.plot,height = 6,width = 8))}

      if(filetag!=""){ggsave(paste0(dir,"Figures and Data/jpg/",filetag,"_",i,"_MissAssignment~z-nloci.jpg"),temp.plot,height = 6,width = 8)} else
      {ggsave(paste0(dir,paste0("Figures and Data/jpg/",i,"_MissAssignment~z-nloci.jpg"),temp.plot,height = 6,width = 8))}

    }

  #save a complete booklet of all figures as a pdf
    if(filetag!=""){
      pdf(file = paste0(dir,"Figures and Data/pdf/",filetag,"_OutputBooklet.pdf"))
      print(p1);print(h1)
      print(p2)
      print(p3);print(h3)
      print(p4);print(h4)
      print(p5);print(h5)
      for (i in unique(PlotData$class)){
        temp.plot <- ggplot(filter(PlotData,class==i))+
          geom_line(aes(x=level,y=value,col=variable),lwd=1.25)+
          geom_line(aes(x=level,y=value+sd,col=variable),lty=2)+
          geom_line(aes(x=level,y=value-sd,col=variable),lty=2)+
          theme_bw()+scale_color_brewer(palette = "Dark2")+
          theme(legend.position="bottom",strip.background = element_rect(fill="white",colour = "black"))+
          labs(x="Probability threshold",y=paste0("Proportion ",i," misassigned \u00B1 sd"),col="Classification")
        print(temp.plot)
      }
      dev.off()
    } else {
      pdf(file = paste0(dir,"Figures and Data/pdf/",filetag,"_OutputBooklet.pdf"))
      print(p1);print(h1)
      print(p2)
      print(p3);print(h3)
      print(p4);print(h4)
      print(p5);print(h5)
      for (i in unique(PlotData$class)){
        temp.plot <- ggplot(filter(PlotData,class==i))+
          geom_line(aes(x=level,y=value,col=variable),lwd=1.25)+
          geom_line(aes(x=level,y=value+sd,col=variable),lty=2)+
          geom_line(aes(x=level,y=value-sd,col=variable),lty=2)+
          theme_bw()+scale_color_brewer(palette = "Dark2")+
          theme(legend.position="bottom",strip.background = element_rect(fill="white",colour = "black"))+
          labs(x="Probability threshold",y=paste0("Proportion ",i," misassigned \u00B1 sd"),col="Classification")
        print(temp.plot)
      }
      dev.off()
    }

    ## clean workspace
    rm(list=setdiff(ls(), c("p1","p2","p3","p4","p5","h1","h3","h4","h5",
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