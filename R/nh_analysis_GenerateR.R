#' @name nh_analysis_GenerateR
#' @title NewHybrids analysis file maker
#'
#' @description \code{nh_analysis_GenerateR} Merges simulated genotypes with the genotypes of unknown/experimental individuals, producing a file to be analyzed by NewHybrids. Will also output a dataframe containing the names of the individuals (including those that were simulated) in the NewHybrids formatted file.
#' @param ReferencePopsData A file path to a NewHybrids formatted file containing genotypes from the simulated ancestral populations. This can be the result of any of the freqbasedsim functions, or a file created using the function genepop_XXXX from the package genepopedit
#' @param UnknownIndivs A dataframe of the genotypes of the individuals to be analyzed for possible hybrid ancestry. This can either be a genepop format file, or a NewHybrids format file. Note - the number of loci in ReferencePopsData must equal the number in UnkownIndivs
#' @param outputName outputName an optioanal character vector to be applied as the name of the output. The default is NULL, in which case the output name is constructed from the name of the input, with the suffix _SiRj_NH added where i is the number of simulations corresponding to the output, and j is the number of replicates of the ith simulation. NH refers to the fact that the output is in NewHybrids format
#' @param indiv.file
#' @param sample.size sample.size an integer number of simulated individuals to be created for each of the six hybrid classes
#'    (e.g. 200 * each of Pure1, Pure2, F1, F2, BC1 and BC2 = 1200 total simulated individuals); default is 200
#' @param NumSims an integer number of simulated datasets to be created; default is 1
#' @param NumReps an integer number of replicates to be created for each of the n simulated datasets specified
#'    in NumSims; default is 1
#' @param cats.include Optional character vector list denoting which hybrid categories should be included in the output; default is all categories.
#' @export
#' @importFrom tidyr separate
#' @importFrom stringr str_split str_extract str_detect
#' @import plyr


nh_analysis_GenerateR <- function(ReferencePopsData, UnknownIndivs, sim.pops.include = c("Pure1", "Pure2"), output.name){
  ### read in teh simulated data
  sim.file <- read.table(ReferencePopsData, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

  path.start <- getwd()

  ## check if the simulated data is GENEPOP or NewHybrids format. This will make a difference.

  header.sim <- sim.file[1,] ## if it is genepop, it will have a single entry in the first position

  if(str_detect(string = header.sim, pattern = "NumIndivs")==FALSE){
    cats <- genepopedit::genepop_detective(GenePop = ReferencePopsData) ### get the names of the populations -- not sure if strictly needed - ask Ryan if can use numeric pop ID
    writeLines("GENEPOP format detected for SIMULATED DATA. Assuming hybrid category order = Pure 1, Pure 2, F1, F2, Back Cross to Pure 1, Back Cross to Pure 2") ### warn that assuming this order

    pop.no.convert <- c("Pure1", "Pure2", "F1", "F2", "BC1", "BC2") ### make a dataframe that matches up to the order of hybrid categories assumed

    inds.get <- which(pop.no.convert %in% sim.pops.include) ### numeric value of which pops assumed match those requested

    genepopedit::subset_genepop(GenePop = ReferencePopsData, keep = TRUE, sPop = inds.get, path = paste0(path.start, "/", "sim.subset.txt")) ## subset

    sim.inds.include <- genepopedit::genepop_flatten(GenePop = paste0(path.start, "/", "sim.subset.txt")) ### read back in and flatten
    sim.inds.include <- sim.inds.include[,-c(2,3)]
    file.remove(paste0(path.start, "/", "sim.subset.txt"))  ### remove the file that was made by subset_genepop
    sim.inds.include.vector <- sim.inds.include[,1] ### get a vector of individual IDs

    sim.inds.Loci <- colnames(sim.inds.include)

      }

  ## if the input file is NewHybrids format, it should have two items in the first row
  if(str_detect(string = header.sim, pattern = "NumIndivs")==TRUE){

    sim.file <- read.table(ReferencePopsData, header = FALSE, skip = 4, stringsAsFactors = FALSE) ## read it in, but skip the first 4 rows because these are not needed - makes a flattened DF
    sim.inds.Loci <- sim.file[1,] ### the first row will ahve the loci names,
    sim.file <- sim.file[-1,] ## remove the loci names
    colnames(sim.file) <- sim.inds.Loci ## add them back in as column names

    NHResultsDir_Split <- unlist(str_split(string = ReferencePopsData, pattern = "/")) ### need to get the directory in which the file is so can get the idnvidual file to get the number of inds in each cat
    NHResultsDir_Split <- NHResultsDir_Split[-grep(x = NHResultsDir_Split, pattern = ".txt")]
    NHResultsDir <- paste0(paste(NHResultsDir_Split, collapse = "/"), "/")
    get.files.list <- list.files(NHResultsDir)
    indiv.file <- read.table(paste0(NHResultsDir, "/", get.files.list[grep(x = get.files.list, pattern = "individuals")])) ## read in the individual file

    Output <- n_class(x = paste0(NHResultsDir, "/", get.files.list[grep(x = get.files.list, pattern = "individuals")])) ## get the # of inds in each cat


    ### Need to determine the range of rows that represent each hybrid category, the subset the requested individuals
    Pure1 <- Output[1,2]
    Pure2 <- Output[2,2]
    F1 <- Output[3,2]
    F2 <- Output[4,2]
    BC1 <- Output[5,2]
    BC2 <- Output[6,2]

    Pure1.inds <- 1:Pure1
    Pure2.inds <- (Pure1 + 1):(Pure1 + Pure2)
    F1.inds <- (Pure1 + Pure2 + 1):(Pure1 + Pure2 + F1)
    F2.inds <- (Pure1 + Pure2 + F1 + 1):(Pure1 + Pure2 + F1 + F2)
    BC1.inds <- (Pure1 + Pure2 + F1 + F2 + 1):(Pure1 + Pure2 + F1 + F2 + BC1)
    BC2.inds <- (Pure1 + Pure2 + F1 + F2 + BC1 + 1):sum(Output$n)

    pop.location.vec <- list(Pure1.inds, Pure2.inds, F1.inds, F2.inds, BC1.inds, BC2.inds)

    Output$Class <- c("Pure1", "Pure2", "F1", "F2", "BC1", "BC2")


    inds.get <- which(Output$Class %in% sim.pops.include)
    inds.get.subset.vec <- unlist(pop.location.vec[inds.get])

    sim.inds.include <- sim.file[inds.get.subset.vec,]
    sim.inds.include.vector <- indiv.file[inds.get.subset.vec, 1]

      } ## END IF Simulated data is NH format

  ### end of input section for simulated data


  ### meow read in the unknown/experimental data

  ## as was done for the simulated data, need to check if entry is a NewHybrids or GENEPOP format file
  unknown.file <- read.table(UnknownIndivs, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

  header.unknown <- unknown.file[1,]
  if(str_detect(string = header.unknown, pattern = "NumIndivs")==FALSE){ ### if a GenePop format file then will have a single entry in the first row
    unknown.indivs.exist <- genepopedit::genepop_detective(GenePop = UnknownIndivs, variable = "Inds") ## get a list of individuals
    pops.exist <- genepopedit::genepop_detective(GenePop = UnknownIndivs) ##
    ag.frame <- data.frame(Exits=pops.exist, ag.to = rep("Pop1", times = length(pops.exist)))

    genepopedit::subset_genepop_aggregate(GenePop = UnknownIndivs, keep = TRUE, agPopFrame = ag.frame, path = paste0(path.start, "/", "unknown.agged.txt"))

    unknown.flattened <- genepopedit::genepop_flatten(GenePop = paste0(path.start, "/", "unknown.agged.txt"))
    unknown.flattened <- unknown.flattened[,-c(2,3)]
    unknown.inds.include <- unknown.flattened
    unknown.Loci <- colnames(unknown.flattened)

    file.remove(paste0(path.start, "/", "unknown.agged.txt"))

}




  #### if it is a NewHybrids format file
  if(str_detect(string = header.unknown, pattern = "NumIndivs")==TRUE){

    unknown.file <- read.table(UnknownIndivs, header = FALSE, skip = 4, stringsAsFactors = FALSE) ## skip the first 4 lines, will build these after anyways
    unknown.Loci <- unknown.file[1,] ## the loci are in the first row
    unknown.file <- unknown.file[-1,] ### remove the first row, these are the loci names - not needed here
    colnames(unknown.file) <- unknown.Loci ## now make them the column names
    unknown.inds.include <- unknown.file ### data to include
    ## if the data are read in as a NH file, then there should be an associated individual file - modify the path to the NH file to get the individual file
    NHResultsDir_Split <- unlist(str_split(string = ReferencePopsData, pattern = "/"))
    NHResultsDir_Split <- NHResultsDir_Split[-grep(x = NHResultsDir_Split, pattern = ".txt")]
    NHResultsDir <- paste0(paste(NHResultsDir_Split, collapse = "/"), "/")
    get.files.list <- list.files(NHResultsDir)
    unknown.indivs.exist <- read.table(paste0(NHResultsDir, "/", get.files.list[grep(x = get.files.list, pattern = "individuals")])) ### hold the individual file to appened to teh simulated individuals

    Output <- n_class(x = paste0(NHResultsDir, "/", get.files.list[grep(x = get.files.list, pattern = "individuals")])) ## also want to have the numbers of individuals in each population

    }



### error check that the simulated individuals and the unknown individuals have the same number of alleles - if not, fail and return error message

  if(length(setdiff(unknown.Loci[-1], sim.inds.Loci[-1])) > 0){stop("The Simulated and Unknown datasets must contain the same marker names.")}



  ###
  indivs.in.dataset <- c(sim.inds.include.vector, unknown.indivs.exist)
  insertNumIndivs <- paste("NumIndivs", length(indivs.in.dataset))

  insertNumLoci <- paste("NumLoci", length(sim.inds.Loci[-1])) ## will probably have to be -1

  ### hard coded stuff
  insertDigits <- "Digits 3"
  insertFormat <- "Format Lumped"


  LociNames <- paste(sim.inds.Loci[-1], collapse = " ")
  insertLociName <- paste("LocusNames", LociNames)

  insert.meta.data <- c(insertNumIndivs, insertNumLoci, insertDigits, insertFormat, insertLociName)

    sim.unknown.combined <- rbind(sim.inds.include[,-1], unknown.inds.include[,-1])
    sim.ind.renameforNH <- c(1:nrow(sim.unknown.combined))
    sim.unknown.combined <- data.frame(sim.ind.renameforNH, sim.unknown.combined)
    sim.unknown.output <- do.call(paste, c(data.frame(sim.unknown.combined[,]), sep = " "))

    data.out <- c(insert.meta.data, sim.unknown.output)

    write(x = data.out, file = output.name)

}