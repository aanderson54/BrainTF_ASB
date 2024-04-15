#!/usr/bin/env R
#Construct_fasta_variant_sequences.R

################################################################################
#This script is used to convert a phased VCF to a table that can more easily be
#parsed and used by a later script, Construct_fasta_variant_sequences.R
#
#Run under R 3.4.3, but can be used under other versions so long as the
#     parallel package is available.
#
#This program takes the following arguments:
# theVariants: Path to the parsed phased VCF table produced by VCF_to_VCFTable.R.
# OUTDIR: Directory to which sequences will be written.
# genomeFile: path to a fasta of the genome being used for this build.
# hangoverLength: how far to extend beyond the last variant when constructing
#     variant regions. This number also denotes how far away two variants
#     can be to be incorporated into the same variant region,
#     and whether or not two nearby variants of different phases will be
#     removed from consideration or not. Recommended that this be equal to
#     read length +1
#
################################################################################


################################################################################
################################################################################
#Load Libraries and define options.
################################################################################
################################################################################

library(parallel)
options("scipen"=100)


################################################################################
################################################################################
#Define Functions
################################################################################
################################################################################


#In order to run this program in parallel, it works best to write the entire
#thing in a very large main function, and put within that function all code
#required to run the program within that function.  This is because using the
#parallel package starts independent sessions that don't remember any
#functions or libraries called outside of that session.

mainFunction <- function(theChromosome, theVariants, outputSequences, problemSequences, problemVariants, genomeFile, hangoverLength) {

  #Load in two relevant packages used in the programs below.
  library(Rsamtools)
  library(GenomicRanges)

  ##############################################################################
  #The following function merges overlapping regions in a pre-formatted file
  #this facilitates getting genomic regions which will contain all of the variants
  #within {hangoverLength}bp of one another.
  ##############################################################################
  mergeRegions <- function(regions, hangoverLength) {
    if(nrow(regions)==1) {
      finalRegions <- regions
      return(finalRegions)
    }
    startLine <- 1
    nextLine <- startLine+1
    finalRegions <- c()
    while(startLine<=nrow(regions)) {
      if(startLine==nrow(regions) || nextLine >= nrow(regions)) {
        regions[startLine,"stop"] <- regions[nextLine,"stop"]
        finalRegions <- rbind(finalRegions, regions[startLine,])
        startLine <- nextLine+1
        nextLine <- startLine+1
        break()
      }
      if(as.numeric(regions[startLine,"stop"])>=(as.numeric(regions[nextLine,"start"])+(hangoverLength-1))) {
        regions[startLine,"stop"] <- regions[nextLine,"stop"]
        nextLine <- nextLine+1
      }
      else {
        finalRegions <- rbind(finalRegions, regions[startLine,])
        startLine <- nextLine
        nextLine <- startLine+1
      }
    }
    return(finalRegions)
  }

  ##############################################################################
  #The following function is used to introduce relevant variants into each region
  #If an error occurs at some step, it returns "NA", otherwise it returns the
  #modified sequences.
  ##############################################################################
  introduceVariants <- function(relevantVars, theSequence) {
    patSequence <- strsplit(theSequence, split="")[[1]]
    matSequence <- strsplit(theSequence, split="")[[1]]
    patTitleAddition <- ""
    matTitleAddition <- ""

    if(length(grep("/", relevantVars[,"Genotype"], fixed=T))>0 || length(grep(".", relevantVars[,"Genotype"], fixed=T))>0 ) { return(NA) }

    else {
      for (m in 1:nrow(relevantVars)) {
        theAlt <- strsplit(relevantVars[m,"Alt"], split=",")[[1]]
        allAlleles <- c(relevantVars[m,"Ref"], theAlt)
        theGenotype <- strsplit(relevantVars[m,"Genotype"], split="|", fixed=T)[[1]]
        theGenotype <- as.numeric(theGenotype)+1

        patReplace <- allAlleles[theGenotype[1]]
        matReplace <- allAlleles[theGenotype[2]]

        startLoc <- as.numeric(relevantVars[m,2])-as.numeric(theseRanges[j,3])+1
        endLoc <- startLoc+nchar(allAlleles[1])-1
        if(endLoc>length(patSequence)) { return(NA) }
        if(paste(toupper(patSequence[startLoc:endLoc]), collapse="")==toupper(allAlleles[1]) && paste(toupper(matSequence[startLoc:endLoc]), collapse="")==toupper(allAlleles[1])) {
          patSequence[startLoc:endLoc] <- rep("", endLoc-startLoc+1)
          patSequence[startLoc] <- patReplace
          matSequence[startLoc:endLoc] <- rep("", endLoc-startLoc+1)
          matSequence[startLoc] <- matReplace

          if(patReplace==matReplace) {theState <- "HOM"}
          else {theState <- "HET"}
          patTitleAddition <- paste(patTitleAddition, relevantVars[m,2], ",", patReplace, ",", theState, ";", sep="")
          matTitleAddition <- paste(matTitleAddition, relevantVars[m,2], ",", matReplace, ",", theState, ";", sep="")

        }
        else{ return(NA) }
      }
    }
    patSequence <- paste(patSequence, collapse="")
    matSequence <- paste(matSequence, collapse="")
    return(c(patSequence, matSequence, patTitleAddition, matTitleAddition))
  }

  ##############################################################################
  #Set up appropriate files for output sequences and any problem variant/sequences.
  #These problems are of two varieties: 1) Unphased variants and those within range
  #of them; 2) two variants (phased or not) which come from different phase sets and
  #are in range of one another.
  #Also, identify the genome
  ##############################################################################
  outputTitles <- paste(outputSequences, theChromosome, "_Titles.txt", sep="")
  outputSequences <- paste(outputSequences, theChromosome, ".fasta", sep="")
  problemSequences <- paste(problemSequences, theChromosome, ".fasta", sep="")
  problemVariants <- paste(problemVariants, theChromosome, ".txt", sep="")
  genome <- FaFile(file=genomeFile)


  #Restrict the variants only to the particular chromosome under consideration.
  theVariantsChrom <- theVariants[theVariants[,"Chrom"]==theChromosome,]

  #Identify all unique Phase Sets, and handle them separately.
  thePhaseSets <- unique(theVariantsChrom[,"PhaseSet"])

  ##############################################################################
  #Loop through each of the phase sets and identify all groups of variants such that
  #each variant is within 100bp of the next variant in the phase set; Watch for
  #situtions wherein a variant is unphased.  Also check in each case if a
  #variant from another phase set is within {hangoverLength-1}bp of a variant of interest.
  ##############################################################################
  for (i in 1:length(thePhaseSets)) {
    if(i%%1000==0) { print(paste(theChromosome, i, length(thePhaseSets), sep="   "))}
    theseVariants <- theVariantsChrom[theVariantsChrom[,"PhaseSet"]==thePhaseSets[i],]
    theseRanges <- c()

    #Set up the range for each variant.
    for (j in 1:nrow(theseVariants)) {
      start <- as.numeric(theseVariants[j,2])-(hangoverLength-1)
      stop <- as.numeric(theseVariants[j,2])+(hangoverLength-1)
      theseRanges <- rbind(theseRanges, c(theseVariants[j,1], "+", start, stop))
    }
    colnames(theseRanges) <- c("chrom", "strand", "start", "stop")

    #Merge the ranges together, where possible.
    theseRanges <- mergeRegions(theseRanges, hangoverLength)
    theseRanges_GR <- makeGRangesFromDataFrame(as.data.frame(theseRanges))

    #for each merged range, get the genomic sequence.
    theSequences <- getSeq(genome, theseRanges_GR, as.character=T)
    theSequences <- as.character(theSequences)

    #For each of the resulting sequences, construct a title for the fasta.
    for (j in 1:length(theSequences)) {
      patTitle <- paste(">", theseRanges[j,1], ":", theseRanges[j,3], "-", theseRanges[j,4], ":P:", sep="")
      matTitle <- paste(">", theseRanges[j,1], ":", theseRanges[j,3], "-", theseRanges[j,4], ":M:", sep="")

      relevantVars <- theVariantsChrom[as.numeric(theVariantsChrom[,2])>=as.numeric(theseRanges[j,3]),]
      relevantVars <- relevantVars[as.numeric(relevantVars[,2])<= as.numeric(theseRanges[j,4]),]

      #If, among the variants within the range of this sequence, variants come
      #from more than one phase set, write this as a problem sequence.
      if(length(unique(relevantVars[,"PhaseSet"]))>1) {
        write.table(cbind(relevantVars, rep(j, nrow(relevantVars)), rep("InRangeOfAlternatePhaseSet", nrow(relevantVars))), problemVariants, append=T, sep="\t", quote=F, row.names=F, col.names=F)
        write(paste(patTitle, theSequences[j], matTitle, theSequences[j], sep="\n"), problemSequences, append=T, sep="\t")
        next()
      }

      #Introduce the variants into the sequence of interest.
      returnedSequence <- NA
      returnedSequence <- introduceVariants(relevantVars, as.character(theSequences[j]))

      #Write the variant sequences.
      if(is.na(returnedSequence[1])) {
        write.table(cbind(relevantVars, rep(j, nrow(relevantVars)), rep("UnresolvedGenotype", nrow(relevantVars))), problemVariants, append=T, sep="\t", quote=F, row.names=F, col.names=F)
        write(paste(patTitle, theSequences[j], matTitle, theSequences[j], sep="\n"), problemSequences, append=T, sep="\t")
      }
      if(!is.na(returnedSequence[1])) {
        write(paste(patTitle, returnedSequence[1], matTitle, returnedSequence[2], sep="\n"), outputSequences, append=T, sep="\n")
        write(paste(paste(patTitle, returnedSequence[3], sep=""), paste(matTitle, returnedSequence[4], sep=""), sep="\n"), outputTitles, append=T, sep="\n")
      }

    }

  }

}



################################################################################
################################################################################
#Begin Script
################################################################################
################################################################################

args=commandArgs(trailingOnly=T)
theVariants <- args[1]
OUTDIR <- args[2]
genomeFile <- args[3]
hangoverLength <- as.numeric(as.character(args[4]))


################################################################################
#Identify the variants of interest and assign appropriate column names.
#Detect the number of cores and start a cluster.
################################################################################
theVariants <- read.table(theVariants, header=F, sep="\t", stringsAsFactors=F)
colnames(theVariants) <- c("Chrom", "Position", "ID", "Ref", "Alt", "PhaseSet", "Genotype")


numCores <- detectCores()
theClust <- makeCluster(numCores)
theChromosomes <- unique(theVariants[,"Chrom"])


problemVariants <- paste(OUTDIR, "/ProblemVariants_", sep="")
problemSequences <- paste(OUTDIR, "/ProblemSequences_", sep="")
outputSequences <- paste(OUTDIR, "/AllVar_Sequences_", sep="")

################################################################################
#Apply the Main function across all chromosomes.
################################################################################
clusterApply(theClust, x=theChromosomes, fun <- mainFunction, theVariants=theVariants, outputSequences=outputSequences, problemSequences=problemSequences, problemVariants=problemVariants, genomeFile=genomeFile, hangoverLength=hangoverLength)
