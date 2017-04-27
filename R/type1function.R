#' @name type1function
#' @title Function to evaluate type I error.
#' @description Subordinate function to hybridpowercomp used in dplyr::do(.) when evaluating simulation type I error.
#' @param x input dataframe created during hybridpowercomp
#' @importFrom dplyr filter
#' @rdname type1function
#' @export
#'

type1function <- function(x){

  x <- as.data.frame(dplyr::ungroup(x))
  x <- dplyr::filter(x,known %in% c("P1","P2"))
  x[which(as.character(x[, "known"]) == as.character(x[, "max.class"])),"domatch"]=TRUE
  x$whichmax <- apply(x[,5:10],1,function(x){x[which.max(x)]})
  x[which(as.numeric(x$whichmax) < as.numeric(x$pofz)),"isgood"]=FALSE

  good.res <- x[x$max.class %in% c("P1","P2") & x$isgood, ]
  bad.res <- x[!x$max.class %in% c("P1","P2") & x$isgood, ]

  x$pos <- FALSE
  x[x$max.class %in% c("P1","P2") & x$isgood,"pos"] <- TRUE

  numsim <- length(unique(x$sim))

  temp.hold.matrix <- matrix(nrow = 2, ncol = numsim)
  colnames(temp.hold.matrix) <- unique(x$sim)

  good.table <- table(good.res$sim)
  good.temp <- t(data.frame(good.table))
  colnames(good.temp) <- good.temp[1,]

  bad.table <- table(bad.res$sim)
  bad.temp <- t(data.frame(bad.table))
  colnames(bad.temp) <- bad.temp[1,]

  temp.hold.matrix <- matrix(nrow = 2, ncol = numsim)
  temp.hold.matrix[1, ] <- 0

  temp.hold.matrix[2, ] <- 0

  match.cols.bad <- which(colnames(good.temp)%in%colnames(bad.temp))
  match.cols.good <- which(colnames(good.temp)%in%unique(x$sim))

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

  return(data.frame(Prop= typeI))

}
