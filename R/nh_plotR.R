#' @name nh_plotR
#' @title NewHybrids Plot
#'
#' @description \code{nh_plotR} plots the cumulative probabilities of assignment for each individual
#' @param NHResults A file path to the NewHybrids (PofZ file) result to be plotted
#' @param ColourVector A logical which can be used to reverse the Pure1/BC1, Pure2/BC2 colours if NEWHYBRIDS switches which population is called Population1 and Population2 between analyses of the same data (allows you to make your plots the same colours)
#' @export
#' @import ggplot2
#' @importFrom reshape2 melt


 nh_plotR <- function(NHResults, ColourVector = 1){ ## this function plots the Q-value output from New Hybrids


   nh_output <- read.table(file = NHResults, header = TRUE)[,-2]

  ## rename the columns so looks better/easier to interpret
   if(ColourVector == 1){
  colnames(nh_output) <- c("Indv", "Pure1", "Pure2","F1", "F2", "BC1", "BC2")
   }

   if(ColourVector == 2){
  colnames(nh_output) <- c("Indv", "Pure2", "Pure1","F1", "F2", "BC2", "BC1")
   }


  nh_output$Indv <-  as.factor(nh_output$Indv) # make the individuals factors to block out the data

  NH_melt <- melt(data = nh_output, id.vars = "Indv") ## melt the data to allow the data to be stacked by indivudal
  colnames(NH_melt) <- c("Indv", "PopProb", "CumProb") ## rename so that its prettier
  NH_melt$PopProb <- factor(x = NH_melt$PopProb, levels = c("Pure1", "Pure2", "F1", "F2", "BC1", "BC2"))
  ## lets give the plot some pretty colours

  col.vec <- c("red", "blue", "grey", "green", "black", "yellow", "brown")


  ## to be used later if decide to break the data or if there are different numbers of individuals in each population type
  # break.by <- nrow(nh_output)/6
  # break.vec <- c(break.by, (break.by*2), (break.by*3), (break.by*4), (break.by*5), (break.by*6))

  ## make a nice pretty little plot
  pretty.plot <- ggplot(NH_melt, aes(x = Indv, y=CumProb, fill=PopProb))

    pretty.plot+geom_bar(stat="identity", position = "stack") + scale_fill_manual(values=col.vec)+ylab("Cumulative Probability")+xlab("Individual") +
      scale_y_continuous(limits = c(0, 1.05), expand=c(0, 0)) +
        theme(axis.title.x = element_text(size = 20,colour = "black"), axis.text.x = element_blank(),
          axis.text.y = element_text(size = 20, colour = "black"), axis.title.y = element_text(size = 20, colour = "black"),
          panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"))


}