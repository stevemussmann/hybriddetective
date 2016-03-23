#' @name n class
#' @title Evalute the number of simulated individuals per NewHybrid class#'
#' @description \code{n_class} will provide a count of each class according to vector of sample names specified in the *_individuals.txt output from  \code{nh_analysis_data_generatoR}
#' @param x vector of sample IDs provided by the *_Individuals.txt output from \code{nh_analysis_data_generatoR}. This must be a vector and not a df.
#' @rdname n_class
#' @export


n_class <- function(x)
  {

  if(length(x)==1){dat <- read.table(x,stringsAsFactors = F);dat <- as.character(dat[,1])}else
    {dat <- as.character(x)}

   #  character vector

  #Identify the specific classes
  Pures <- dat[grep("Pure",dat)]
  Pures <- unlist(lapply(Pures,function(x)unlist(strsplit(x,"_"))[2]))
  BC <- dat[grep("BC",dat)]
  BC <- unlist(lapply(BC,function(x)unlist(strsplit(x,"_"))[2]))
  F1 <- dat[grep("F1",dat)]
  F2 <- dat[grep("F2",dat)]
  #Create dataframe which details the number in each class
  Output <- data.frame(Class=c("P1","P2","F1","F2","BC1","BC2"),
                       n=c(table(factor(Pures, levels=unique(Pures)))[1],
                           table(factor(Pures, levels=unique(Pures)))[2],
                           length(F1),
                           length(F2),
                           table(factor(BC, levels=unique(BC)))[1],
                           table(factor(BC, levels=unique(BC)))[2]))

  #return output
  return(Output)

  }