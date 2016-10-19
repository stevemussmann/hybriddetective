#' @name accuracyfunction
#' @title Function to evaluate simulation accuracy.
#' @description Subordinate function to hybridpowercomp used in dplyr::do(.) when evaluating simulation accuracy.
#' @rdname accuracyfunction
#' @export
#'

accuracyfunction <- function(x){
  x = data.frame(ungroup(x))
  x[which(as.character(x[, "known"]) == as.character(x[, "max.class"])), "domatch"] = TRUE
  x$whichmax <- apply(x[, 5:10], 1, function(x){x[which.max(x)]})
  x[which(as.numeric(x$whichmax) < as.numeric(x$pofz)), "isgood"] = FALSE

  ## filter to retain only instances where assignment matches known category, and the PofZ is greater than the critical
  x1 <- x[x$isgood & x$domatch, ]
  x1$known <- factor(x = x1$known, levels = c("P1", "P2", "F1", "F2", "BC1", "BC2")) ## change factor levels

  ## filter so the assignment is retained (above PofZ), no matter if it is correct (matches known) or not
  x2 <- x[x$isgood, ]

  ## NEED to look at the number of individuals assigned to each category by: simulation, replicate, and number of loci
  ### WANT to average WITHIN simulation
  ## table() will create counts of individuals
  ## CREATE tables that are simulations by replicates, within each genotype category, within each number of loci
  ## table[simulation, replicate, category, loci]

  sumtable <- table(x1$sim,x1$rep) ### CORRECT TABLE - assigned correctly at critical PofZ (numerator)
  sumtable2 <- table(x2$sim,x2$rep) ### Assigned TABLE - assigned to group at critical PofZ (denominator)

  num.sims <- length(unique(x1$sim)) ### how man simulations were conducted - allows to be variable, and will be populated from the data provided

  means <- apply(X = sumtable, MARGIN = 1, FUN = mean, na.rm = TRUE)/apply(X = sumtable2, MARGIN = 1, FUN = mean, na.rm = TRUE)

  return(data.frame(simulation=as.character(rownames(sumtable)), means = means))

}