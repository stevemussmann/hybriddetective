#' @name nh_build_example_results
#' @title Creates a folder containing example NewHybrids analyses which can be evaluated using hybridpowercomp
#' @description Creates a folder containing example NewHybrids analyses which can be evaluated using hybridpowercomp
#' @param dir File path to the directory in write the NewHybrids results (in individual folders as returned by parallelNH_XX). The default is NULL, in which case the folder will be written to the working directory
#' @param remove_example A logical indicating if the example data should be removed from the users hard drive. Default = FALSE. NOTE: If a dir was previously specified, it must be specified again
#' @rdname nh_build_example_results
#' @export


nh_build_example_results <- function(dir = NULL, remove_example = FALSE){

  ## If the directory is not speficied, set dir to the working directory
  if(length(dir) < 1){
    dir = getwd()
  }

  if(remove_example == FALSE){

  ## load in the data
  inds <- data("inds")
  sim1_rep1 <- data("sim1_rep1")
  sim1_rep2 <- data("sim1_rep2")
  sim1_rep3 <- data("sim1_rep3")
  sim2_rep1 <- data("sim2_rep1")
  sim2_rep2 <- data("sim2_rep2")
  sim2_rep3 <- data("sim2_rep3")

  ## Create Results folder
  dir.create(paste0(dir, "nh.results"))
  ## create folders for example data in results folder
  dir.create(paste0(dir, "nh.results/sim1_rep1"))
  dir.create(paste0(dir, "nh.results/sim1_rep2"))
  dir.create(paste0(dir, "nh.results/sim1_rep3"))
  dir.create(paste0(dir, "nh.results/sim2_rep1"))
  dir.create(paste0(dir, "nh.results/sim2_rep2"))
  dir.create(paste0(dir, "nh.results/sim2_rep3"))


  write(inds, paste0(dir, "nh.results/sim1_rep1/example_individuals.txt"))
  write(inds, paste0(dir, "nh.results/sim1_rep2/example_individuals.txt"))
  write(inds, paste0(dir, "nh.results/sim1_rep3/example_individuals.txt"))
  write(inds, paste0(dir, "nh.results/sim2_rep1/example_individuals.txt"))
  write(inds, paste0(dir, "nh.results/sim2_rep2/example_individuals.txt"))
  write(inds, paste0(dir, "nh.results/sim2_rep3/example_individuals.txt"))


  write.table(sim1_rep1, paste0(dir, "nh.results/sim1_rep1/example_sim1_rep1_PofZ.txt"),row.names = FALSE, quote = FALSE, col.names = FALSE)
  write.table(sim1_rep2, paste0(dir, "nh.results/sim1_rep2/example_sim1_rep2_PofZ.txt"),row.names = FALSE, quote = FALSE, col.names = FALSE)
  write.table(sim1_rep3, paste0(dir, "nh.results/sim1_rep3/example_sim1_rep3_PofZ.txt"),row.names = FALSE, quote = FALSE, col.names = FALSE)
  write.table(sim2_rep1, paste0(dir, "nh.results/sim2_rep1/example_sim2_rep1_PofZ.txt"),row.names = FALSE, quote = FALSE, col.names = FALSE)
  write.table(sim2_rep2, paste0(dir, "nh.results/sim2_rep2/example_sim2_rep2_PofZ.txt"),row.names = FALSE, quote = FALSE, col.names = FALSE)
  write.table(sim2_rep3, paste0(dir, "nh.results/sim2_rep3/example_sim2_rep3_PofZ.txt"),row.names = FALSE, quote = FALSE, col.names = FALSE)
  }

  if(remove_example == TRUE){

    unlink(paste0(dir, "nh.results"), recursive = TRUE)

  }


}