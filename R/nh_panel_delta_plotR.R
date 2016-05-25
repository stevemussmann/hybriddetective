#' @name nh_panel_delta_plotR
#' @title Plot change in individual assignment probability at different panel sizes
#' @description Evaluates the accuracy with which NewHybrids assigns individuals of known hybrid class to the correct hybrid class in simulated datasets at varying levels of stringency (PofZ). The code will write graphical and numerical results to the directory provided by the user.
#' @param GPD Filepath to the directory which holds the output from different runs through New Hybrids (e.g. 3 simulations with 3 replicate runs each through NH) note that this directory should only hold the output folders.
#' @param return.workspace A logical query of whether the ggplot2 object should be returned to the workspace for further editing or formatting. Can take values of "TRUE" or "FALSE". NOTE: if return.workspace = TRUE, the user must assign a variable to the function (i.e. mydata <- nh_ind_panel_delta(...))
#' @param save.plot A logical query of whether the plot should be saved as a file in the working directory. NOTE: The file will be saved "delta_plot.xxx", and so must be renamed if the function is to be run more than once to prevent overwriting.
#' @param plot.filetype If save.plot = TRUE, plot.filetype specifies the filetype to save the plot as. The default is "png", with the option of "pdf".
#' @rdname nh_panel_delta_plotR
#' @import ggplot2
#' @import magrittr
#' @importFrom dplyr filter summarise ungroup group_by
#' @importFrom stringr str_extract
#' @importFrom reshape2 melt
#' @export



nh_panel_delta_plotR <- function(GPD, return.workspace, save.plot = FALSE, plot.filetype = "png"){
samplesize = NULL

#set directory for which holds the New Hybrids output folders
  filedir <- GPD
  dirt <- GPD
  lfiles <- setdiff(list.files(dirt),c("Figures and Data", "NewHybrids Plots")) #ignores Figures folder in case this is run more than once and in case plots made

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
              numLociExt <- stringr::str_extract(string = LandAfile, pattern = paste0("from ", "[:digit:]{1,5}", " loci"))
              numLociWorking <- gsub(x = numLociExt, pattern = "from ", replacement = "")
              numLociWorking <- as.numeric(gsub(x = numLociWorking, pattern = " loci", replacement = ""))

              #identify the simulation and repeat info
                S_ident <- gsub("_","",stringr::str_extract(pzfile,paste0("_S","[:digit:]{1}","R","[:digit:]{1}","_")))
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



              output <- rbind(output,tempfile)

          }#end of for loop




    ## average and SD the  replicate runs of each simulation in New Hybrids. Filter is just a holder for the dplyr:: call
      sim_data <- as.data.frame(dplyr::filter(output)%>%group_by(nLoci,sim,Indv)%>%dplyr::summarise(Pure1_sd=sd(Pure1),Pure1=mean(Pure1),
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
      }
    }

    boxdata$nLoci=factor(boxdata$nLoci)

    #Combined loci
    sim_means2 <- sim_means
    sim_means2$hybrid <- rowSums(sim_means2[,c("F1","F2","BC1","BC2")])
    sim_means2[which(sim_means2$class =="Pure1"),"hybrid"]=sim_means2[which(sim_means2$class =="Pure1"),"Pure1"] #add values of the Pure
    sim_means2[which(sim_means2$class =="Pure2"),"hybrid"]=sim_means2[which(sim_means2$class =="Pure2"),"Pure2"] #add values of the Pure

    sim_means2$hclass <- "Hybrid"
    sim_means2[which(sim_means$class=="Pure1"),"hclass"] <- "Pure1"
    sim_means2[which(sim_means$class=="Pure2"),"hclass"] <- "Pure2"

    sim_means2$hclass <- factor(sim_means2$hclass,levels=c("Pure1","Pure2","Hybrid"))


  plot.data <- sim_means2

# plot.data <- read.csv(GPD, header = TRUE)



plot.x.dim = 8
plot.y.dim = 10
x.axis.title.fontsize = 20
x.axis.label.fontsize = 12
x.axis.title = "Panel Size"


tmelt <- melt(data = plot.data, id.vars = c("nLoci", "sim", "Indv"), measure.vars = c("Pure1", "Pure2", "F1", "F2", "BC1", "BC2"))


delta.plot <- ggplot(tmelt, aes(x = factor(nLoci), y = as.factor(Indv))) + geom_tile(aes(fill = variable, alpha = value, height =2)) + scale_fill_manual(values = to.colour) + labs(x = x.axis.title) +
facet_grid(.~ variable) + theme(panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"),
  axis.title.y = element_blank(), axis.text.y = element_blank(), legend.position = "none", axis.ticks.y = element_blank(), panel.grid = element_blank(),
  axis.title.x = element_text(colour = "black", size = x.axis.title.fontsize), axis.text.x = element_text(colour = "black", size = x.axis.label.fontsize))

if(return.workspace == TRUE){
  return(delta.plot)
}


if(save.plot == TRUE){

  if(plot.filetype == "png"){ggsave(delta.plot, filename = "delta_plot.png", height = plot.y.dim, width = plot.y.dim, units = "in", dpi = 600)}
  if(plot.filetype == "pdf"){ggsave(delta.plot, filename = "delta_plot.png", height = plot.y.dim, width = plot.y.dim, units = "in", dpi = 600)}

}



} ### END FUNCTION

