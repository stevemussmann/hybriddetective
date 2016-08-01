#' @name nh_analysis_simulateR_generateR
#' @title NewHybrids data generator
#'
#' @description \code{nh_analysis_data_generateR} Quickly creates simulated reference populations from user provided data, and merges these with the genotypes of unknown/experimental individuals, producing a file to be analyzed by NewHybrids. Will also output a dataframe containing the names of the individuals (including those that were simulated) in the NewHybrids formatted file. The simulation function is the same as in freqbasedsim
#' @param ReferencePopsData A file path to a GenePop formatted file containing genotypes from two (2) ancestral populations. This is the data from which the simulated hybrids will be constructed
#' @param UnknownIndivs A dataframe of genotypes of unknown/experimental individuals to be analyzed for hybrid category generated using the function *genepop_flatten* from the package *genepopedit*. Individuals or loci can be subsetted from this flattened data set if desired, or manipulations can conducted before flattening using the appropriate functions from *genepopedit*. Note - the number of loci in ReferencePopsData must equal the number in UnkownIndivs
#' @param outputName outputName A character vector to be applied as the name of the output.
#' @param pop.groups Optional character vector denoting how the two ancestral populations should be named; default is PopA and PopB
#' @param sample.size an integer number of simulated individuals to be created for each of the six hybrid classes (viz. Pure1, Pure2, F1, F2, BC1, BC2). The default is 200 (viz. 200 * each of Pure1, Pure2, F1, F2, BC1 and BC2 = 1200 total simulated individuals)
#' @param NumSims an integer number of simulated datasets to be created. The default is 1
#' @param NumReps  integer number of replicates of each of the NumSims simulated dataset to be created. The default is 1
#' @param cats.include Optional character vector list denoting which hybrid categories should be included in the output; default is all categories.
#' @export
#' @importFrom tidyr separate
#' @importFrom stringr str_split str_extract str_detect
#' @import plyr




nh_analysis_simulateR_generateR <- function(ReferencePopsData, UnknownIndvs, outputName = NULL, pop.groups = c("Pure1", "Pure2"), sample.size = 200, NumSims = 1, NumReps = 1, cats.include = c("PopA", "PopB", "F1", "F2", "BCA", "BCB")){

  if(length(out.name) == 0){out.name=ReferencePopsData}

  ## Top part of this function has been direcrtly copied from the simulate hybrids function


  ## Get to work!
      ## Load the data - this will create the centred Pure1 and Pure2 individuals against which the UNKNOWN INDIVIUDUALS will be compared
      ReferencePops <-  read.table(ReferencePopsData, header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)

      ## get the name of the data being simulated - will be used to name the output if outputname isn't specified
      GPsplit <- c(stringr::str_split(string = ReferencePopsData, pattern = "/")) ## split by backslash
      outNameHold <- stringr::str_extract(GPsplit, paste0("[:word:]{3,}", ".txt")) ## get the name to call the simulated files from the file name entered - hold here to affix Reps and Sims later

      ## NewHybrids requires the number of individuals being analyzed be specified
      NumIndivsSIM <- sample.size*2 # <- this is just the number to be simulated/group * number of groups [6]


    ## Stacks version information
        stacks.version <- ReferencePops[1,] # this could be blank or any other source. ## this was duplicated from another function - not sure if needed

    #Remove first label of the stacks version
        ReferencePops <- as.vector(ReferencePops)
        ReferencePops <- ReferencePops[-1,]

    #Add an index column to ReferencePops and format as a dataframe
        ReferencePops <- data.frame(data=ReferencePops,ind=1:length(ReferencePops))
        ReferencePops$data <- as.character(ReferencePops$data)

    #ID the rows which flag the Populations
        Pops  <-  which(ReferencePops$data == "Pop" | ReferencePops$data =="pop" | ReferencePops$data == "POP")
        npops  <-  1:length(Pops)

    ## Seperate the data into the column headers and the rest
        ColumnData <- ReferencePops[1:(Pops[1]-1),"data"]  ### SNP Names
        NumLoci <- length(ColumnData) ### NewHybrids Requires the number of LOCI be specified

        snpData <- ReferencePops[Pops[1]:NROW(ReferencePops),]  ### Genotypes - this is where the magic starts

    #Get a datafile with just the snp data no pops
        tempPops <- which(snpData$data=="Pop"| snpData$data =="pop" | snpData$data == "POP")
        snpData <- snpData[-tempPops,]

    #Seperate the snpdata
    #First we pull out the population data which follows "TEXT ,  "
        temp <- separate(snpData,data,into=c("Pops","snps"),sep=",")
        temp$snps <- substring(temp$snps,3) # delete the extra spaces at the beginning
        temp2 <- data.frame(do.call(rbind, str_extract_all(temp$snps, "[0-9]{3}")))

    ## Going to have to break the two alleles of the SNPS apart - this will thus double the number of columns
    ## SO <- will want to have SNP_A and SNP_A2
        ColumnData2 <- ColumnData ## Duplicatet the SNP names
        ColumnData2 <- paste(ColumnData2, "2", sep = ".")    ## add .2 to each duplicated name

    ## can't just append the duplicated names to the end of the original names - have to intersperse them
        places = rep(1:length(ColumnData)*2) ### creates a list of even numbers 2X as long as the number of columns
  ##i.e. the lenght of the original plus the duplicated names
  ## - this will also mark the position to insert the duplicates
        ColumnData.Dup = rep(NA, times = length(ColumnData)*2) ### make an object to feed names into
      for(i in 1:length(ColumnData)){ ### for loop to add original, then duplicated name
          a = places[i]-1 ## Original names go first, and are in the odd positions
          b = places[i] ### Duplicated names go second and are in the even posoitons
          Col.name.orig = ColumnData[i] ## Get the name in the ith position
          Col.name.plus2 = ColumnData2[i] ## get the name in the ith positon
          ColumnData.Dup[a] = Col.name.orig ## add the original name to the new vector
          ColumnData.Dup[b] = Col.name.plus2 ## add the duplicate name to the new vector
              } ## End of Loop


    #Contingency to see if R read in the top line as the "stacks version" -- modified to deal with the duplicated SNP names
    if (length(temp2)!=length(ColumnData.Dup)){colnames(temp2) <- c(stacks.version, paste(stacks.version, "2", sep = "."),ColumnData.Dup)}
    if (length(temp2)==length(ColumnData.Dup)){colnames(temp2) <- ColumnData.Dup}
    #if (length(temp2)/2!=length(ColumnData)){stacks.version="No stacks version specified"}

    ## Get the Alpha names
    NamePops=temp[,1] # Names of each

    if(length(pop.groups) == 0){ ### If unique grouping IDs ≠ number of "Pop" user must give vector of groupings
                                ### equal to number of "Pop" or else the function will fail
    NameExtract=stringr::str_extract(NamePops, "[A-z]{3,}" ) ### if looking at higher order grouping (i.e. pops in
        # regions) can have more unique coding than "Pop" - will want to remove original names so can
        ## keep track of which unique groupings cross. i.e. Cross by "Pop", but remember ID of parents
          } ## End of IF statement


  # extract the text from the individuals names to denote population
  ## Now add the population tags using npops (number of populations and Pops for the inter differences)
    tPops <- c(Pops,NROW(ReferencePops))
      PopIDs <- NULL
          for (i in 2:length(tPops)){
            hold <- tPops[i]-tPops[i-1]-1
            if(i==length(tPops)){hold=hold+1}
            pophold <- rep(npops[i-1],hold)
            PopIDs <- c(PopIDs,pophold)
          } ## end of loop

    temp2$Pop <- PopIDs;


     if(length(pop.groups)!=0){
     hold.names=stringr::str_extract(NamePops, "[A-z]{3}" ) ## This may need to be improved in published version
        for(i in 1:length(unique(PopIDs))){
          u.ID.no <- unique(PopIDs)[i]
          to <- min(which(PopIDs==u.ID.no))
          from <- max(which(PopIDs==u.ID.no))
      hold.names[to:from] = paste(pop.groups[i], hold.names[to:from], sep=".")
    }
    NameExtract <- hold.names
     }

     ## get the nubmer of indivudals within each "Pop" grouping --- the Number of individuals in the two ancesntal populations need not be the same as the nubmer of
        ## individuals to be simulated
    PopLengths <- table(temp2$Pop)

     ## Need to be able to tell what row each individual is in, and what population it is
        ind.vector = c(1:nrow(temp)) ### make a vector that is the number of individuals
        ind.matrix = data.frame(temp2$Pop, ind.vector) ## add populatuions to that

      temp.split <- split(x = temp2, f = temp2$Pop)

      pop.recall <- NULL
          for(i in 1:length(temp.split)){
            popn <- paste(pop.groups[i], "pop", sep = "_")
            temp.split.hold = temp.split[[i]]
            temp.split.hold = temp.split.hold[-which(names(temp.split.hold) == "Pop")]
            assign(x = popn, value = temp.split.hold)
            pop.recall <- c(pop.recall, popn)
                } ## End of loop


      mat.name.recall <- NULL
      for(i in 1:length(pop.recall)){
          temp.mat <- data.frame(matrix(vector(), 2, length(temp2)/2))
          pop.get <- get(pop.recall[i])

            temp.mat.hold <- NULL
              for(k in 1:nrow(pop.get)){
                ind.hold <- pop.get[k,]
                temp.mat[1,] <- t(t(ind.hold[c(T,F)]))
                temp.mat[2,] <-  t(t(ind.hold[c(F,T)]))
                temp.mat.hold <- rbind(temp.mat.hold, temp.mat)
                    } ## End of K loop

            mat.out.name <- paste(pop.recall[i],"matrix", sep = "_")
            assign(x = mat.out.name, value = temp.mat.hold)
            mat.name.recall <- c(mat.name.recall, mat.out.name)
                } ## End of I loop

      ##### Loooooooooooop tha Sims! #####

               for(sim in 1:NumSims){

                    ### MAKE PURE CROSS - centre the data -
                    pure.name.recall <- NULL

                      for(k in 1:length(pop.groups)){
                          pop1 <- get(mat.name.recall[k])
                          pop2 <- get(mat.name.recall[k])

                          off.interspersed.out <- NULL
                          for(i in 1:sample.size){

                            hold.off.pop1 <- apply(pop1, FUN = sample, 2, 1)
                            hold.off.pop2 <- apply(pop2, FUN = sample, 2, 1)
                            hold.off.interspersed <- data.frame(c(rbind(hold.off.pop1, hold.off.pop2)))

                            off.interspersed.out <- rbind(off.interspersed.out, t(hold.off.interspersed))
                                  } ## End of I loop

                          pure.name <- paste("Pure", pop.groups[k], sep = "_")
                          pure.name.recall <- c(pure.name.recall, pure.name)

                          assign(x = pure.name, value = off.interspersed.out)
                              } ## End of K loop


              inv.pure.name.recall <- NULL
              for(i in 1:length(pop.recall)){
                temp.mat <- data.frame(matrix(vector(), 2, length(temp2)/2))
                pop.get <- get(pure.name.recall[i])

                temp.mat.hold <- NULL
                  for(k in 1:nrow(pop.get)){
                    ind.hold <- pop.get[k,]
                    temp.mat[1,] <- t(t(ind.hold[c(T,F)]))
                    temp.mat[2,] <-  t(t(ind.hold[c(F,T)]))
                    temp.mat.hold <- rbind(temp.mat.hold, temp.mat)
                      } ## End of K loop

                inv.pure.out.name <- paste(pure.name.recall[i],"inv", sep = "_")
                assign(x = inv.pure.out.name, value = temp.mat.hold)
                inv.pure.name.recall <- c(inv.pure.name.recall, inv.pure.out.name)
                      } ## end of i loop

      ### MAKE F1 CROSS

          pop1 <- get(inv.pure.name.recall[1])
          pop2 <- get(inv.pure.name.recall[2])
          F1.out <- NULL
              for(i in 1:sample.size){

                  hold.off.pop1 <- apply(pop1, FUN = sample, 2, 1)
                  hold.off.pop2 <- apply(pop2, FUN = sample, 2, 1)
                  hold.off.interspersed <- data.frame(c(rbind(hold.off.pop1, hold.off.pop2)))

                  F1.out <- rbind(F1.out, t(hold.off.interspersed))
                      } ## END of loop


          temp.mat <- data.frame(matrix(vector(), 2, length(temp2)/2))
          pop.get <- F1.out

      inv.F1 <- NULL
      for(k in 1:nrow(pop.get)){
        ind.hold <- pop.get[k,]
        temp.mat[1,] <- t(t(ind.hold[c(T,F)]))
        temp.mat[2,] <-  t(t(ind.hold[c(F,T)]))
        inv.F1 <- rbind(inv.F1, temp.mat)
        } ## END K LOOP


        ### MAKE F2 CROSS
          pop1 <- inv.F1
          pop2 <- inv.F1
          F2.out <- NULL
            for(i in 1:sample.size){

                hold.off.pop1 <- apply(pop1, FUN = sample, 2, 1)
                hold.off.pop2 <- apply(pop2, FUN = sample, 2, 1)
                hold.off.interspersed <- data.frame(c(rbind(hold.off.pop1, hold.off.pop2)))

                F2.out <- rbind(F2.out, t(hold.off.interspersed))
                } ## END I LOOP


  ### MAKE Back CROSS


  BC.name.recall <- NULL
  for(k in 1:length(pop.groups)){

    pop1 <- get(inv.pure.name.recall[k])
    pop2 <- inv.F1

  off.interspersed.out <- NULL
  for(i in 1:sample.size){

    hold.off.pop1 <- apply(pop1, FUN = sample, 2, 1)
    hold.off.pop2 <- apply(pop2, FUN = sample, 2, 1)
    hold.off.interspersed <- data.frame(c(rbind(hold.off.pop1, hold.off.pop2)))

    off.interspersed.out <- rbind(off.interspersed.out, t(hold.off.interspersed))
      } ## END I LOOP

    BC.name <- paste("BC", pop.groups[k], sep = "_")
    BC.name.recall <- c(BC.name.recall, BC.name)

    assign(x = BC.name, value = off.interspersed.out)
      } ## END K LOOP


  for(i in 1:length(pure.name.recall)){
    off.name <- paste(pure.name.recall[i], c(1:sample.size), sep="_")
    hold.dat <- get(pure.name.recall[i])
    hold.dat <- data.frame(off.name, hold.dat)
    ColumnData.Dup.insert = c("ID", ColumnData.Dup)
    colnames(hold.dat) = ColumnData.Dup.insert
    assign(x = pure.name.recall[i], value = hold.dat)
  } ## END I LOOP


  f1.off.name <- paste("F1", c(1:sample.size), sep = "_")
  F1.out <- data.frame(f1.off.name, F1.out)
  colnames(F1.out) <- c("ID", ColumnData.Dup)

  f2.off.name <- paste("F2", c(1:sample.size), sep = "_")
  F2.out <-  data.frame(f2.off.name, F2.out)
  colnames(F2.out) <- c("ID", ColumnData.Dup)

  BC.name.recall
  for(i in 1:length(BC.name.recall)){
    off.name <- paste(BC.name.recall[i], c(1:sample.size), sep="_")
    hold.dat <- get(BC.name.recall[i])
    hold.dat <- data.frame(off.name, hold.dat)
    ColumnData.Dup.insert = c("ID", ColumnData.Dup)
    colnames(hold.dat) = ColumnData.Dup.insert
    assign(x = BC.name.recall[i], value = hold.dat)
  }


  for(b in 1:length(pure.name.recall)){

    fam.to.bind.name <- pure.name.recall[b]
    fam.to.bind <- get(fam.to.bind.name)
    indiv.hold <- fam.to.bind[,1]
    loci.bind <- which(stringr::str_detect(string = colnames(fam.to.bind), pattern = "\\.2")==TRUE)
      col.out <- NULL

      for(k in 1:length(loci.bind)){
        place.1 <- (loci.bind[k]-1)
        place.2 <- loci.bind[k]
        hold.col <- paste0(fam.to.bind[,place.1], fam.to.bind[,place.2])
        col.out <- cbind(col.out, hold.col)

          }

    fam.reord <- cbind(indiv.hold,col.out)
    colnames(fam.reord) <- c(colnames(fam.to.bind[1]), colnames(fam.to.bind[c((loci.bind-1))]))
    assign(x = fam.to.bind.name, fam.reord)


  for(b in 1:length(pure.name.recall)){

    fam.to.remove.untyped.name <- pure.name.recall[b]

    fam.to.remove.untyped <- get(fam.to.remove.untyped.name)
    fam.to.remove.untyped[which(stringr::str_detect(string = fam.to.remove.untyped, pattern = "000")==TRUE)] = "000000"
    assign(x = fam.to.remove.untyped.name, value = fam.to.remove.untyped)
  }

      }

  ## F1

  fam.to.bind.name <- "F1.out"

  fam.to.bind <- get(fam.to.bind.name)
  indiv.hold <- fam.to.bind[,1]
  loci.bind <- which(stringr::str_detect(string = colnames(fam.to.bind), pattern = "\\.2")==TRUE)
    col.out <- NULL

    for(k in 1:length(loci.bind)){
      place.1 <- (loci.bind[k]-1)
      place.2 <- loci.bind[k]
      hold.col <- paste0(fam.to.bind[,place.1], fam.to.bind[,place.2])
      col.out <- cbind(col.out, hold.col)

        }

  fam.reord <- cbind(indiv.hold,col.out)
  colnames(fam.reord) <- c(colnames(fam.to.bind[1]), colnames(fam.to.bind[c((loci.bind-1))]))
  assign(x = fam.to.bind.name, fam.reord)


  fam.to.remove.untyped.name <- "F1.out"

  fam.to.remove.untyped <- get(fam.to.remove.untyped.name)
  fam.to.remove.untyped[which(stringr::str_detect(string = fam.to.remove.untyped, pattern = "000")==TRUE)] = "000000"
  assign(x = fam.to.remove.untyped.name, value = fam.to.remove.untyped)



  fam.to.bind.name <- "F2.out"

  fam.to.bind <- get(fam.to.bind.name)
  indiv.hold <- fam.to.bind[,1]
  loci.bind <- which(stringr::str_detect(string = colnames(fam.to.bind), pattern = "\\.2")==TRUE)
    col.out <- NULL

    for(k in 1:length(loci.bind)){
      place.1 <- (loci.bind[k]-1)
      place.2 <- loci.bind[k]
      hold.col <- paste0(fam.to.bind[,place.1], fam.to.bind[,place.2])
      col.out <- cbind(col.out, hold.col)

        }

  fam.reord <- cbind(indiv.hold,col.out)
  colnames(fam.reord) <- c(colnames(fam.to.bind[1]), colnames(fam.to.bind[c((loci.bind-1))]))
  assign(x = fam.to.bind.name, fam.reord)


  fam.to.remove.untyped.name <- "F2.out"

  fam.to.remove.untyped <- get(fam.to.remove.untyped.name)
  fam.to.remove.untyped[which(stringr::str_detect(string = fam.to.remove.untyped, pattern = "000")==TRUE)] = "000000"
  assign(x = fam.to.remove.untyped.name, value = fam.to.remove.untyped)


  for(b in 1:length(BC.name.recall)){

    fam.to.bind.name <- BC.name.recall[b]

    fam.to.bind <- get(fam.to.bind.name)
    indiv.hold <- fam.to.bind[,1]
    loci.bind <- which(stringr::str_detect(string = colnames(fam.to.bind), pattern = "\\.2")==TRUE)
      col.out <- NULL

      for(s in 1:length(loci.bind)){
        place.1 <- (loci.bind[s]-1)
        place.2 <- loci.bind[s]
        hold.col <- paste0(fam.to.bind[,place.1], fam.to.bind[,place.2])
        col.out <- cbind(col.out, hold.col)

          }

    fam.reord <- cbind(indiv.hold,col.out)
    colnames(fam.reord) <- c(colnames(fam.to.bind[1]), colnames(fam.to.bind[c((loci.bind-1))]))
    assign(x = fam.to.bind.name, fam.reord)

  }

  for(b in 1:length(BC.name.recall)){

    fam.to.remove.untyped.name <- BC.name.recall[b]

    fam.to.remove.untyped <- get(fam.to.remove.untyped.name)
    fam.to.remove.untyped[which(stringr::str_detect(string = fam.to.remove.untyped, pattern = "000")==TRUE)] = "000000"
    assign(x = fam.to.remove.untyped.name, value = fam.to.remove.untyped)
  }

  #Now recompile the NewHybrids


  pop.names <- c(pure.name.recall, "F1.out", "F2.out", BC.name.recall)


  sim.pop.name.reference <- c("PopA", "PopB", "F1", "F2", "BCA", "BCB")
  pop.names.use <- pop.names[which(sim.pop.name.reference %in% cats.include)]


  ### make vector of names of simulated individuals
  popvecout <- NULL
  for(i in 1:length(pop.names.use)){

    pvecmake <- paste0(pop.names.use[i], "_", c(1:nrow(get(pop.names.use[i]))))

    popvecout <- c(popvecout, pvecmake)

  }


  ### compile simulated data
  sim.out <- NULL
  for(i in 1:length(pop.names.use)){
    hold.pop <- get(pop.names.use[i])
    sim.out <- rbind(sim.out, hold.pop)
  }


  sim.out[,1] <- c(1:nrow(sim.out))

  Loci.sim <- do.call(paste, c(data.frame(sim.out[,]), sep = " "))
  cd2 <- paste(ColumnData, collapse = " ")

  NumIndivs.SIM <- length(Loci.sim)
  NumIndivs.UNKOWN <- nrow(unknown.pop)
  NumIndivs <- NumIndivs.SIM + NumIndivs.UNKOWN


  unknown.indv.vec <- droplevels(unknown.pop[,1])

  ## make output individual vector

  output.indv.vec <- c(popvecout, as.character(unknown.indv.vec))
  output.indv.df <- data.frame(c(1:length(output.indv.vec)), output.indv.vec)
  names(output.indv.df) <- c("number", "ID")

  unknown.pop.use <- unknown.pop
  unknown.pop.use[,1] <- c((length(Loci.sim)+1):NumIndivs)
  Loci.Unknown <- do.call(paste, c(data.frame(unknown.pop.use[,]), sep = " "))

  ##### here is where the UnknownIndvs are added

  insertNumIndivs <- paste("NumIndivs", NumIndivs)
  insertNumLoci <- paste("NumLoci", NumLoci)
  insertYourDigits <- "Digits 3"
  insertFormat <- "Format Lumped"
  insertLociName <- paste("LocusNames", cd2)

  Loci.out <- c(insertNumIndivs, insertNumLoci, insertYourDigits, insertFormat, insertLociName,   Loci.sim, Loci.Unknown)
  #write.table(x = Loci.out, file = "NSDrop1.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

   outNameGive <- out.name
  outNameGive <- paste0(outNameGive, "_S", sim)

  for(reps in 1:NumReps){

    outNameGive <- paste0(outNameGive, "R", reps, "_analdata.txt")
    write.table(x = Loci.out, file = outNameGive, row.names = FALSE, col.names = FALSE, quote = FALSE)
    outNameGivereRef <- paste0(out.name, "_IDref.txt")
    write.table(x = output.indv.df, file = outNameGivereRef, row.names = FALSE, col.names = FALSE, quote = FALSE)
      } ## END REPS LOOP

    } ### End of sim loop

  } ## End of function

