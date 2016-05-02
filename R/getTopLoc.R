#' @name getTopLoc
#' @title Creates a panel of the top n unlinked loci, and exports the list of loci
#'
#' @description \code{getTopLoc} Extracts the genotypes of individuals at the top n (by Fst) unlinked loci. The resultant file can then be used to simulate the panel efficacy using freqbasedsim. Alternatively, multiple panel sizes can be created by using truncated vectors of the exported loci names in the function subset_genepop from the genepopedit package (www.github.com/rystanley)
#' @param GPD A file path to the GENEPOP format file you wish to create your panel from
#' @param LDpop A string which denotes which of the two populations you wish to calculate linkage disequilibrium in. The options are "Pop1" or "Pop2", or "Both" if the LD is to be calculated based on both populations.
#' @param panel.size An integer number of loci to include in the panel
#' @param r2.threshold The minimum r^2 threshold to consider a pair of loci to be in LD
#' @param ld.window Number of adjacent SNPs to compare each SNP against for LD - default is NULL, which translates to a window size of 99999, which essentially asks to compare each SNP against all others
#' @param allocate.PGD.RAM An integer value in GB to specify the maximum amount of RAM to allocate to PGDspider. The default is 1 GB, which should be sufficient for most analyses.
#' @param where.PLINK A file path to the PLINK installation folder
#' @param where.PGDspider A file path to the PGDspider installation folder
#' @export
#' @importFrom genepopedit genepop_detective subset_genepop genepop_fstat
#' @importFrom hierfstat read.fstat wc
#' @importFrom stringr str_split str_extract
#' @import plyr


getTopLoc <- function(GPD, LDpop = "Pop1", panel.size, where.PLINK, where.PGDspider, return.environment = TRUE, save.LociandIndividuals = FALSE){

  writeLines("Reading Data")

  path.start <- getwd()  ### where to write the files created by genepopedit to

  pops.exist <- genepopedit::genepop_detective(GPD) ## see what populations are in the file

  ##panel.size must be an integer

  if(panel.size%%1 != 0){
    httr::BROWSE("https://en.wikipedia.org/wiki/Integer")
    stop("Panel size must be an integer.")
  }

  ## Function must have exactly two populations to work - fail if not
  if(length(pops.exist) != 2){
    stop("File must contain two populations. The package genepopedit offers editing functions")
  }

  #Variable checks
  if(allocate.PGD.RAM%%1 != 0){
    stop("Please specify an integer GB value to allocate to PGDspider.")
  }

  allocate.PGD.RAM <- allocate.PGD.RAM*1024

  if(r2.threshold < 0 | r2.threshold > 1){
    stop("r^2 threshold must be a value between 0 and 1")
  }

  if(r2.threshold<0.2){
    writeLines("Linkage detection threshold is low (<0.2). Linkage will be classified at a higher frequency than default PLINK selection parameters.")
  }

  if(length(ld.window)==0){
    ld.window = 99999 ### sets the LD window to essentially check every SNP pairwise
  }

  if(ld.window < 0){
    stop("LD window must be non-negative")
  }

  writeLines("Creating training and working datasets")
  ### subsample to get a training and simulated dataset
  inds.sub <- genepopedit::genepop_sample(GenePop = GPD, nsample = 0.5)

  genepopedit::subset_genepop_individual(GenePop = GPD, indiv = inds.sub, keep = TRUE, path = paste0(path.start, "/", "GPD_for_sim.txt"))
  genepopedit::subset_genepop_individual(GenePop = GPD, indiv = inds.sub, keep = FALSE, path = paste0(path.start, "/", "GPD_for_GET_TOP_LOC.txt"))

  sim.path <- paste0(path.start, "/", "GPD_for_sim.txt")
  GPD.Top <- paste0(path.start, "/", "GPD_for_GET_TOP_LOC.txt")



  ### some people might want to get LD for both pops at once - sure why not give them the option


  if(LDpop != "Both"){

  ## makes a dataframe to match whatever name is given to the populations to Pop1 and Pop2
popLDsubsetDF <- data.frame(op=pops.exist, rename=c("Pop1", "Pop2")) ## make a DF to pick the right pop from the file

subPOP <- popLDsubsetDF[which(popLDsubsetDF[,2] == LDpop),1] ### get the name of the pop to be subsetted

subPOP <- as.character(subPOP)

## subset out the population in which LD is to be calculated - this will make a file, which will be deleted after
genepopedit::subset_genepop(GenePop = GPD.Top, sPop = subPOP, keep = TRUE, path = paste0(path.start, "/", "subset_for_LD.txt"))
## remember path to the file created by subset_genepop
sub_data_path <- paste0(path.start, "/", "subset_for_LD.txt")


  }


  if(LDpop == "Both"){

    popLDsubsetDF <- data.frame(op=pops.exist, rename=c("Pop1", "Pop1")) ## make a both the same

    genepopedit::subset_genepop_aggregate(GenePop = GPD.Top, agPopFrame = popLDsubsetDF, path = paste0(path.start, "/", "subset_for_LD.txt"))
    sub_data_path <- paste0(path.start, "/", "subset_for_LD.txt")
    ## now rename

    genepopedit::subset_genepop_rename(GenePop = sub_data_path, path = sub_data_path, nameframe = popLDsubsetDF)

  }


writeLines("Calculating Fst")
### change the format of the original file to FSTAT so can get Fst

# modify the path to play nice with spaces

path.start.PGD <- gsub(x = path.start, pattern = " ", replacement = "\\")
where.PGDspider.PGD <- gsub(x = where.PGDspider, pattern = " ", replacement = "\\ ", fixed = TRUE)

GP_FSTAT_SPID <- "# spid-file generated: Fri Apr 08 10:53:23 ADT 2016

# GENEPOP Parser questions
PARSER_FORMAT=GENEPOP

# Enter the size of the repeated motif (same for all loci: one number; different: comma separated list (e.g.: 2,2,3,2):
GENEPOP_PARSER_REPEAT_SIZE_QUESTION=
# Select the type of the data:
GENEPOP_PARSER_DATA_TYPE_QUESTION=SNP
# How are Microsat alleles coded?
GENEPOP_PARSER_MICROSAT_CODING_QUESTION=REPEATS

# FSTAT Writer questions
WRITER_FORMAT=FSTAT

# Specify which data type should be included in the FSTAT file  (FSTAT can only analyze one data type per file):
FSTAT_WRITER_DATA_TYPE_QUESTION=SNP
# Save label file
FSTAT_WRITER_LABEL_FILE_QUESTION=
# Do you want to save an additional file with labels (population names)?
FSTAT_WRITER_INCLUDE_LABEL_QUESTION=false
# Specify the locus/locus combination you want to write to the FSTAT file:
FSTAT_WRITER_LOCUS_COMBINATION_QUESTION=
"

write(x = GP_FSTAT_SPID, file = paste0(path.start, "/", "GP_FSTAT.spid"))
remember.FSTAT.spidpath.WD <- paste0(path.start, "/", "GP_FSTAT.spid")


### move spid file to the PGDspider folder
file.copy(from = paste0(path.start, "/GP_FSTAT.spid"), to = where.PGDspider, overwrite = TRUE)
remember.FSTAT.spidpath.PGD <- paste0(where.PGDspider, "/", "GP_FSTAT.spid")
## move the input file as well to the same location as PGDspider - this makes this step so much easier
file.copy(from <- GPD.Top, to = where.PGDspider, overwrite = TRUE)
remember.GPD.Top.PGD <- paste0(where.PGDspider, "GPD_for_GET_TOP_LOC.txt")


#### OS X and LINUX CALL

if(Sys.info()["sysname"] != "Windows"){

### create a string to call PGDspider
input.file.call <- "-inputfile GPD_for_GET_TOP_LOC.txt"
execute.SPIDER <- paste0("java -Xmx", allocate.PGD.RAM, "m -Xms512m -jar PGDSpider2-cli.jar")
spid.call <- "-spid GP_FSTAT.spid"
input.format <- "-inputformat GENEPOP"
output.format <- "-outputformat FSTAT"
goto.spider <- paste0("cd ", where.PGDspider.PGD, "; ", execute.SPIDER)
output.file.path <- "-outputfile for_FST.txt"
## string to run
run.PGDspider <- paste0(goto.spider, " ", input.file.call, " ", input.format, " ", output.file.path, " ", output.format, " ", spid.call)


### run PGDspider through system
system(run.PGDspider)

remember.for_FST.PGD <- paste0(where.PGDspider, "for_FST.txt")

} # End MAC LINUX IF


#### Windows call

if(Sys.info()["sysname"] == "Windows"){

### create a string to call PGDspider
input.file.call <- "-inputfile GPD_for_GET_TOP_LOC.txt"
execute.SPIDER <- paste0("java -Xmx", allocate.PGD.RAM, "m -Xms512m -jar PGDSpider2-cli.jar")
spid.call <- "-spid GP_FSTAT.spid"
input.format <- "-inputformat GENEPOP"
output.format <- "-outputformat FSTAT"
goto.spider <- paste0("cd ", where.PGDspider.PGD, " && ", execute.SPIDER)
output.file.path <- "-outputfile for_FST.txt"
## string to run
run.PGDspider <- paste0(goto.spider, " ", input.file.call, " ", input.format, " ", output.file.path, " ", output.format, " ", spid.call)


### run PGDspider through system
shell(run.PGDspider)

remember.for_FST.PGD <- paste0(where.PGDspider, "for_FST.txt")

} # End WINDOWS IF



### move the FSTAT format file back to the working directory
file.copy(from = paste0(where.PGDspider, "/for_FST.txt"), to = path.start, overwrite = TRUE)
## remember the path of the file created by genepop_fstat
fst_data_path <- paste0(path.start, "/", "for_FST.txt")


### read in the FSTAT formatted file
for.fst <- hierfstat::read.fstat("for_FST.txt")
## calculate Fst
FST.dat <- suppressWarnings(hierfstat::wc(for.fst))
### get the Fst values
FSTs <- FST.dat$per.loc$FST
### create a dataframe that is the names of the Loci, and their corresponding Fst
FST.df <- data.frame(colnames(for.fst)[-1], FSTs)
names(FST.df)[1] <- "loci"
## reorder the dataframe from highest to lowest Fst
FST.df <- FST.df[base::order(FST.df$FSTs, decreasing = TRUE),]
# head(FST.df)

writeLines("Calculating Linkage")
### convert file to .ped and .map using PGD spider


## create spid file

spidTop <- "# spid-file generated: Thu Mar 10 09:40:22 AST 2016

# GENEPOP Parser questions
PARSER_FORMAT=GENEPOP

# Enter the size of the repeated motif (same for all loci: one number; different: comma separated list (e.g.: 2,2,3,2):
GENEPOP_PARSER_REPEAT_SIZE_QUESTION=
# Select the type of the data:
GENEPOP_PARSER_DATA_TYPE_QUESTION=SNP
# How are Microsat alleles coded?
GENEPOP_PARSER_MICROSAT_CODING_QUESTION=REPEATS

# PED Writer questions
WRITER_FORMAT=PED

# Save MAP file"

map.loc <- paste0("PED_WRITER_MAP_FILE_QUESTION= ", "PGDtest")

spidBottom <- "# Replacement character for allele encoded as 0 (0 encodes for missing data in PED):
PED_WRITER_ZERO_CHAR_QUESTION=
# Specify the locus/locus combination you want to write to the PED file:
PED_WRITER_LOCUS_COMBINATION_QUESTION=
# Do you want to save an additional MAP file with loci information?
PED_WRITER_MAP_QUESTION=true"

spid.file <- c(spidTop, map.loc, spidBottom)

## write spid file
write(x = spid.file, file = paste0(path.start, "/", "hyb.spid"))
remember.LD.spidpath.WD <- paste0(path.start, "/", "hyb.spid")


### move spid file to the PGDspider folder
file.copy(from = paste0(path.start, "/hyb.spid"), to = where.PGDspider, overwrite = TRUE)
remember.LD.spidpath.PGD <- paste0(where.PGDspider, "/", "hyb.spid")
## move the input file as well to the same location as PGDspider - this makes this step so much easier
file.copy(from <- sub_data_path, to = where.PGDspider, overwrite = TRUE)
remember.sub_data_path.PGD <- paste0(where.PGDspider, "subset_for_LD.txt")


#### OS X LINUX call

if(Sys.info()["sysname"] != "Windows"){
### create a string to call PGDspider
input.file.call <- paste0("-inputfile subset_for_LD.txt")
execute.SPIDER <- paste0("java -Xmx", allocate.PGD.RAM, "m -Xms512m -jar PGDSpider2-cli.jar")
spid.call <- "-spid hyb.spid"
input.format <- "-inputformat GENEPOP"
output.format <- "-outputformat PED"
goto.spider <- paste0("cd ", where.PGDspider.PGD, "; ", execute.SPIDER)
output.file.path <- paste0("-outputfile PGDtest.ped")
## string to run
run.PGDspider <- paste0(goto.spider, " ", input.file.call, " ", input.format, " ", output.file.path, " ", output.format, " ", spid.call)

### run PGDspider through system
system(run.PGDspider)
remember.PEDpath.PGD <- paste0(where.PGDspider, "PGDtest.ped")
remember.MAPpath.PGD <- paste0(where.PGDspider, "PGDtest.map")

} # END OSX LINUX if


#### Windows call

if(Sys.info()["sysname"] == "Windows"){
### create a string to call PGDspider
input.file.call <- paste0("-inputfile subset_for_LD.txt")
execute.SPIDER <- paste0("java -Xmx", allocate.PGD.RAM, "m -Xms512m -jar PGDSpider2-cli.jar")
spid.call <- "-spid hyb.spid"
input.format <- "-inputformat GENEPOP"
output.format <- "-outputformat PED"
goto.spider <- paste0("cd ", where.PGDspider.PGD, " && ", execute.SPIDER)
output.file.path <- paste0("-outputfile PGDtest.ped")
## string to run
run.PGDspider <- paste0(goto.spider, " ", input.file.call, " ", input.format, " ", output.file.path, " ", output.format, " ", spid.call)

### run PGDspider through system
shell(run.PGDspider)
remember.PEDpath.PGD <- paste0(where.PGDspider, "PGDtest.ped")
remember.MAPpath.PGD <- paste0(where.PGDspider, "PGDtest.map")

} # END OSX LINUX if


## move the created ped and map files to the PLINK folder

ped.path <- paste0(where.PGDspider, "/", "PGDtest.ped")
map.path <- paste0(where.PGDspider, "/", "PGDtest.map")

file.copy(from = ped.path, to = where.PLINK, overwrite = TRUE)
file.copy(from = map.path, to = where.PLINK, overwrite = TRUE)

plink_ped_path <- paste0(where.PLINK, "/", "PGDtest.ped")
plink_map_path <- paste0(where.PLINK, "/", "PGDtest.map")


### prepare a string to call PLINK
### modify PLINK path so it plays nice with system
where.PLINK.go <- gsub(x = where.PLINK, pattern = " ", replacement = "\\ ", fixed = TRUE)
go.to.PLINK <- paste0("cd ", where.PLINK.go)


### OSX LINUX PLINK call
if(Sys.info()["sysname"] != "Windows"){
execute.PLINK <- paste0(go.to.PLINK, "; ", "./plink --file PGDtest --r2 --ld-window-r2 ", r2.threshold, " --ld-window ", ld.window, " --noweb")
### run PLINK through system
system(execute.PLINK)
}

### Windows PLINK CALL
if(Sys.info()["sysname"] == "Windows"){
execute.PLINK <- paste0(go.to.PLINK, " && ", "plink --file PGDtest --r2 --ld-window-r2 ", r2.threshold, " --ld-window ", ld.window, " --noweb")
### run PLINK through system
shell(execute.PLINK)
}


## copy the LD file created by PLINK to the working directory
file.copy(from = paste0(where.PLINK, "plink.ld"), to = path.start)

## the format of the LD file is a bit messed up - it is not a regular matrix - have to modify the file a bit to get it to read in properly

## rename the file and change it to a txt
file.rename(from = paste0(path.start, "/", "plink.ld"), to = paste0(path.start, "/", "plink.txt"))
remember.PLINKLDFile.WDpath <- paste0(path.start, "/", "plink.txt")

## read the file in as characters - there are an uneven number of spaces between data, so have to remove those before can read as a table
LDs <- readChar(con = "plink.txt", file.info("plink.txt")$size) ## read in as a long character string
LDs_reform <- gsub(x = LDs, pattern = "\\ +", replacement = " ") ### remove more than 1 space
## save as a new .txt file, which is now formatted to be read into R easily
write(x = LDs_reform, file = "LDsReform.txt")
remember.LDReformat.path <- paste0(path.start, "/LDsReform.txt")
## read back in
Linked <- read.table("LDsReform.txt", header = TRUE)
## keep only the needed columns
Linked <- Linked[c("SNP_A", "SNP_B", "R2")]
### turn both columns into a single vector of unduplicated loci names
ld.unique <- c(as.character(Linked$SNP_A), as.character(Linked$SNP_B))
ld.unique <- as.character(ld.unique[which(duplicated(ld.unique)==FALSE)])

if(length(ld.unique > 1)){

  FST.df2 <- FST.df[which(FST.df$loci %in% ld.unique),]
  FST.ld.ordered <- as.character(FST.df2[order(FST.df2$FSTs,decreasing=T),"loci"])

  holdlist <- list()
  for(i in FST.ld.ordered){
    hold <- c(i,as.character(Linked[which(Linked[,1]%in%i),2]),
              as.character(Linked[which(Linked[,2]%in%i),1]))

    hold2 <- NULL
    for(k in 1:length(hold)){hold2 <- c(hold2,which(FST.ld.ordered == hold[k]))}
    holdlist[[i]] <- hold2
  }

  linked.ranks.df2 <- plyr::rbind.fill(lapply(holdlist,function(y){as.data.frame(t(y),stringsAsFactors=FALSE)}))

  to.keep <- FST.ld.ordered[genepopedit::Optimfunc(linked.ranks.df2)]
  to.drop <- setdiff(as.character(FST.df2$loci),to.keep)

  FST.df2$loci <- as.character(FST.df2$loci)
  FST.df$loci <- as.character(FST.df$loci)

  to.keep <- FST.df[-which(FST.df$loci %in% to.drop),]

  writeLines("Writing output")
  Unlinked.panel <- FST.df[which(FST.df$loci %in% to.keep$loci),]

  your.panel <- FST.df
  your.panel_un <- Unlinked.panel
  your.panel <- your.panel[order(your.panel$FSTs,decreasing=TRUE),]
  your.panel_unlinked <- your.panel_un[order(your.panel_un$FSTs,decreasing = TRUE),]
  your.panel_unlinked <- your.panel_unlinked[1:panel.size,] #return the panel size
  Linked.df <- Linked
  writeLines("Process Completed.")
}

if(length(ld.unique) < 1){ # if not linked loci
  writeLines("Writing output")
  writeLines("No linked loci detected, all loci returned.")
  your.panel <- FST.df
  your.panel <- your.panel[order(your.panel$FSTs,decreasing=TRUE),]
  your.panel_unlinked <- your.panel[1:panel.size,]#return the panel size
  Linked.df <- NA
  }


 if(nrow(your.panel_unlinked) <  panel.size){

   writeLines(paste0("Desired panel size ", panel.size,  " is larger than total number of unlinked loci ", length(FST.order.vec[-to.cut.out]), ". ", "Returning ALL unlinked loci."))
   panel.size = nrow(your.panel_unlinked)

 }


 GPsplit <- c(stringr::str_split(string = GPD, pattern = "/"))

  outNameHold <- stringr::str_extract(GPsplit, paste0("[:word:]{3,}", ".txt"))
  GPD.path <- gsub(x = GPD, pattern = outNameHold, replacement = "")
  outNameHold <- gsub(x = outNameHold, pattern = ".txt", replacement = "")

  outName <- paste0(GPD.path, outNameHold, "_", panel.size, "_Loci_Panel.txt")

 genepopedit::subset_genepop(GenePop = sim.path,
                             subs = as.character(your.panel_unlinked$loci),
                             keep = TRUE, path = outName)

 new.sim.GP <- data.table::fread(outName,
                                 header = FALSE, sep = "\t",
                                 stringsAsFactors = FALSE)

 if(save.LociandIndividuals == TRUE){
 write(x = your.panel_unlinked$loci, file = paste0(GPD.path, paste0("Top_", panel.size, "_loci.txt")))
 write(x = inds.sub, file = paste0(GPD.path, "_individuals_for_Simulation.txt"))
}


 if(return.environment == TRUE){

   nh_topLoci <- list()
   nh_topLoci$One <- new.sim.GP
   nh_topLoci$Two <- your.panel_unlinked$loci
   nh_topLoci$Three <- inds.sub

   names(nh_topLoci) = c(paste0(outNameHold, "_", panel.size, "_Loci_Panel"),
      paste0("Top_", panel.size, "_loci"),
      paste0(outNameHold, "_IndividualsInSim"))

 }


file.remove(sim.path)
file.remove(GPD.Top)
file.remove(sub_data_path)
file.remove(remember.FSTAT.spidpath.WD)
file.remove(remember.FSTAT.spidpath.PGD)
file.remove(remember.GPD.Top.PGD)
file.remove(remember.for_FST.PGD)
file.remove(fst_data_path)
file.remove(remember.LD.spidpath.WD)
file.remove(remember.LD.spidpath.PGD)
file.remove(remember.sub_data_path.PGD)
file.remove(remember.PEDpath.PGD)
file.remove(remember.MAPpath.PGD)
file.remove(plink_ped_path)
file.remove(plink_map_path)
file.remove(remember.PLINKLDFile.WDpath)
file.remove(remember.LDReformat.path)

writeLines("Process Completed.")

 }

