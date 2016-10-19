#' @name type2function
#' @title Function to evaluate simulation type II error.
#' @description Subordinate function to hybridpowercomp used in dplyr::do(.) when evaluating simulation type II error.
#' @importFrom dplyr filter
#' @rdname type2function
#' @export
#'

type2function <- function(x){

  x <- as.data.frame(ungroup(x))
  samplesize <- as.numeric(unlist(strsplit(x[1,"samplesize"],",")))
  classnames <- unlist(strsplit(x[1,"classnames"],","))
  x <- x[,-grep("samplesize|classnames",colnames(x))] #remove extra columns not needed
  x$missclass <- classnames[apply(x[,classnames],1,which.max)] # what is the class of the highest NH probability
  x$missval <- apply(x[,classnames],1,function(x){as.numeric(x[which.max(x)])})
  x1 <- dplyr::filter(x,class!=missclass & missval>=unique(pofz)) #dataset with missclassifications

  xdummy <- data.frame(Var1 = c("Pure1","Pure2","F1","F2","BC1","BC2"),dummy=NA) # dataframe for dummy values

  tempstore <- NULL #storage of the output
  for (z in classnames){
    x2 <- filter(x1,class == z)
    if(nrow(x2)>0){
      x3 <- as.data.frame(table(x2$missclass)/samplesize[which(classnames==z)]) # percentage of samples miss classed to a given class of a given type of class (i)

      x4 <-  merge(xdummy,x3,by="Var1",all.x = TRUE)

      xout <- data.frame(class=z,
                         mclass_P1=x4[which(x4$Var1 == "Pure1"),"Freq"],
                         mclass_P2=x4[which(x4$Var1 == "Pure2"),"Freq"],
                         mclass_F1=x4[which(x4$Var1 == "F1"),"Freq"],
                         mclass_F2=x4[which(x4$Var1 == "F2"),"Freq"],
                         mclass_BC1=x4[which(x4$Var1 == "BC1"),"Freq"],
                         mclass_BC2=x4[which(x4$Var1 == "BC2"),"Freq"])
    } else
    {xout <- data.frame(class=z,
                        mclass_P1=NA,
                        mclass_P2=NA,
                        mclass_F1=NA,
                        mclass_F2=NA,
                        mclass_BC1=NA,
                        mclass_BC2=NA)} #end else

    tempstore <- rbind(tempstore,xout)

  } #end zloop
  return(tempstore)

} #end of TypeII_Function