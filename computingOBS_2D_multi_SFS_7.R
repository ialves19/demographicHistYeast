#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

########################
########################
##
##
## This R script requires: 
##  1. genotype m, strains/clade and sample files to be in the same dir
##  2. the model file is in the output folder
##
##  This script is called using: computingOBS_2D_multi_SFS.7.R <inDir> <outDir> <strains/clade> <GTm> <model> <samples>
##
##  Isabel Alves - March 2024
##
########################
########################

#library(tidyverse, lib.loc = "/shared/home/ialves/R/x86_64-redhat-linux-gnu-library/4.3")
library(tidyverse)
#library(ggpubr)
#library(ggplot2)
set.seed(988599694)

# arguments 
inDir <- args[1]
outDir <- args[2]
dfStrains <- args[3]
GTmName <- args[4]
modelFile <- args[5]
samplesF <- args[6]
boot <- args[7]

# To debug
# inDir <- "/shared/home/ialves/demoHist_yeast3039/04-analysis/fastsimcoal/SFS_134samples"
# outDir <- "/shared/home/ialves/demoHist_yeast3039/04-analysis/fastsimcoal/fiveClades"
# dfStrains <- "df_Clade_Strains.txt"
# GTmName <- "samples27.coded.AA.GT"
# modelFile <- "model_OneDomest.txt"
# samplesF <- "samples27.txt"

# ----------------
# opening files
# ----------------
df_clades <- read.table(paste0(inDir, "/", dfStrains), header = F, sep = "\t")
model <- scan(paste0(outDir, "/", modelFile), nlines = 1, what = character())
struModel <- scan(paste0(outDir, "/", modelFile), skip = 1, nlines = 1, what = numeric())
df_zig <- read.table(paste0(outDir, "/", modelFile), skip = 2, header = F, sep = "\t")
GTmatrix <- read.table(paste0(inDir, "/", GTmName), header = F, sep = "\t")
samples <- scan(paste0(inDir, "/", samplesF), what = character())
colnames(GTmatrix) <- samples

# ----------------
# creating model specific folder
# ----------------
if(boot == F) {
  sfsDir <- paste0(outDir, "/", model, "_dSFS")
} else {
  nb.boot <-  as.numeric(unlist(strsplit(unlist(strsplit(GTmName, split = "/"))[1], split = "_"))[3])
  sfsDir <- paste0(outDir, "/", model, "_", nb.boot, "_dSFS")
  model <- paste0(model, "_", nb.boot)
}

#sfsDir <- paste0(outDir, "/", model, "_dSFS")
if(dir.exists(sfsDir)){
  print("Folder ", sfsDir, " already exists.")
  print(paste0("Generating files in: ", sfsDir))
} else {
  print(paste0("Folder ", sfsDir, " doesn't exist."))
  print(paste0("Creating model folder: ", sfsDir))
  dir.create(sfsDir, showWarnings = T)
}

# ----------------
# transforming clade names 
# ----------------
listt <- lapply(df_clades$V1, function(x) { unlist(strsplit(x, split = "\\."))[1:2]})
cladePerStrainList <- lapply(listt, function(x) { paste0(x, collapse = ".") })
cladePerStrainDF <- do.call(rbind, cladePerStrainList)
df_clades <- df_clades %>% mutate(Clade=cladePerStrainDF) %>% select(-V1)
# converting it into a list
listSamplesPerCladeIDsClean <- list()
for(clade in c(unique(df_clades$Clade))) {
  listSamplesPerCladeIDsClean[[clade]] <- unlist(df_clades %>% filter(Clade == clade) %>% select(V2))
}

# copying GTm
GT_copy <- GTmatrix
# ----------------
# transforming zygosity into a list
# ----------------
zigosityStatus <- as.list(df_zig$V2)
names(zigosityStatus) <- df_zig$V1
cladeNames <- names(zigosityStatus)
# ----------------
# Creating pop IDs
# ----------------
# this takes into account whether or not we have empty populations, ie ghosts
# it should count starting from 1
# the empty deme before +1 is zero-based
emptyDemes <- which(struModel == 0)
nbPopTLPfile <- length(struModel)
popIDs <- rep(1, nbPopTLPfile)
names(popIDs) <- rep("Empty", nbPopTLPfile)
popIDs[emptyDemes] <- 0
names(popIDs)[popIDs == 1] <- cladeNames

# ----------------
# Generating 2D-SFS
# ----------------
pairsClades <- combn(length(cladeNames),2, replace=F)

popTable <- data.frame(cbind(seq(0,(length(cladeNames)-1), by=1), cladeNames))
colnames(popTable) <- c("popID", "CladeNames")
nbSamples <- length(samples)
write.table(popTable, file=paste0(sfsDir, "/popLabels_", nbSamples,"samples_", model,".txt"), row.names = F, col.names = T, quote = F)

sfs_2d_list <- list()
# compute 2D-dSFS
for(cladePair in 1:ncol(pairsClades)) {
  
  #cladePair <- 1
  cladeOne_name <-  cladeNames[pairsClades[1,cladePair]] 
  cladeTwo_name <-  cladeNames[pairsClades[2,cladePair]]
  cat("Computing 2D-SFS for the pair:", cladeOne_name, cladeTwo_name, sep = "\n")
  
  strainsCladeOne <- unlist(listSamplesPerCladeIDsClean[[cladeNames[pairsClades[1, cladePair]]]])
  strainsCladeTwo <- unlist(listSamplesPerCladeIDsClean[[cladeNames[pairsClades[2, cladePair]]]])
  
  sampleSetOne <- colnames(GT_copy[,colnames(GT_copy) %in% strainsCladeOne])
  sampleSetTwo <- colnames(GT_copy[,colnames(GT_copy) %in% strainsCladeTwo])
  
  subGTmOne <- GT_copy[,colnames(GT_copy) %in% strainsCladeOne]
  subGTmTwo <- GT_copy[,colnames(GT_copy) %in% strainsCladeTwo]
  
  if(zigosityStatus[[cladeNames[pairsClades[1,cladePair]]]] == "Homozygous") {
    subGTmOne[subGTmOne == 1] <- sample(c(0,2), 1)  
    subGTmOne[subGTmOne == 2] <- 1
    
    entriesOne <- rep(0, length(sampleSetOne)+1)
    names(entriesOne) <- as.character(0:length(sampleSetOne))
    #entries[names(entries) %in% names(table(apply(subGTm, 1, sum)))] <- table(apply(subGTm, 1, sum))
  } else {
    
    entriesOne <- rep(0, (length(sampleSetOne)*2)+1)
    names(entriesOne) <- as.character(0:(length(sampleSetOne)*2))
    
  }
  
  if(zigosityStatus[[cladeNames[pairsClades[2,cladePair]]]] == "Homozygous") {
    
    subGTmTwo[subGTmTwo == 1] <- sample(c(0,2), 1)  
    subGTmTwo[subGTmTwo == 2] <- 1
    
    entriesTwo <- rep(0, length(sampleSetTwo)+1)
    names(entriesTwo) <- as.character(0:length(sampleSetTwo))
    #entries[names(entries) %in% names(table(apply(subGTm, 1, sum)))] <- table(apply(subGTm, 1, sum))
  } else {
    
    entriesTwo <- rep(0, (length(sampleSetTwo)*2)+1)
    names(entriesTwo) <- as.character(0:(length(sampleSetTwo)*2))
    
  }
  
  countsSetOne <- apply(subGTmOne, 1, sum)
  countsSetTwo <- apply(subGTmTwo, 1, sum)
  
  SFS_2D <- matrix(rep(0, length(entriesOne)*length(entriesTwo)), ncol = length(entriesTwo), nrow=length(entriesOne))
  
  colnames(SFS_2D) <- paste0("d", (which(names(popIDs) == cladeTwo_name)-1), "_", names(entriesTwo))
  rownames(SFS_2D) <- paste0("d", (which(names(popIDs) == cladeOne_name)-1), "_", names(entriesOne))
  
  
  for(snp in 1:length(countsSetOne)) { #length(countsSetOne)
    SFS_2D[as.numeric(countsSetOne[snp])+1,as.numeric(countsSetTwo[snp])+1] <- SFS_2D[as.numeric(countsSetOne[snp])+1,as.numeric(countsSetTwo[snp])+1]+1 
  }
  
  cat("1 observation", file = paste0(sfsDir, "/", model, "_","jointDAFpop", which(names(popIDs) == cladeTwo_name)-1, "_", which(names(popIDs) == cladeOne_name)-1,".tmp"), sep = "\n")
  write.table(t(SFS_2D), file=paste0(sfsDir, "/", model, "_","jointDAFpop", which(names(popIDs) == cladeTwo_name)-1, "_", which(names(popIDs) == cladeOne_name)-1,".tmp"), row.names = TRUE, col.names = NA, quote = F, sep = "\t", append = T)
  
} 
# --------------------------
# ----------------
# ------


# ----------------
# Generating multi-SFS
# ----------------

# creating lists
strainsClade <- list()
sampleSet <-  list()
subGT <-  list()
entries <-  list()
countsSFS <-list()

# this loop goes over the populations, checks it's zygosity levels and transforms the genotype matrix accordingly. 
# this means that for homozygous strains the genotypes encoded as 2 pass to 1 as they are considered haploid. 
# the outcome of this loop is a vector with the counts of derived alleles per site on the genotype matrix 
for(cladeNb in 1:length(cladeNames)) {
  
  #cladeNb <- 1
  strainsClade[[cladeNb]] <- unlist(listSamplesPerCladeIDsClean[[cladeNames[cladeNb]]])
  
  sampleSet[[cladeNb]] <- colnames(GT_copy[,colnames(GT_copy) %in% strainsClade[[cladeNb]]])
  
  subGT[[cladeNb]] <- GT_copy[,colnames(GT_copy) %in% strainsClade[[cladeNb]]]
  
  if(zigosityStatus[[cladeNames[cladeNb]]] == "Homozygous") {
    subGT[[cladeNb]][subGT[[cladeNb]] == 1] <- sample(c(0,2), 1)  
    subGT[[cladeNb]][subGT[[cladeNb]] == 2] <- 1
    
    entries[[cladeNb]] <- rep(0, length(sampleSet[[cladeNb]])+1)
    names(entries[[cladeNb]]) <- as.character(0:length(sampleSet[[cladeNb]]))
    #entries[names(entries) %in% names(table(apply(subGTm, 1, sum)))] <- table(apply(subGTm, 1, sum))
  } else {
    
    entries[[cladeNb]] <- rep(0, (length(sampleSet[[cladeNb]])*2)+1)
    names(entries[[cladeNb]]) <- as.character(0:(length(sampleSet[[cladeNb]])*2))
    
  }
  
  countsSFS[[cladeNb]] <- apply(subGT[[cladeNb]], 1, sum)
  
} # end of loop counting the number of derived alleles within population

a <- array(rep(0,prod(lengths(entries))), dim = lengths(entries))
countMono <- 0
countSing <- 0
for(snp in 1:length(countsSFS[[cladeNb]])) { #length(countsSetOne)
  #snp <- 1    
  coord <- sapply(countsSFS, `[[`, snp)
  idx <- paste(coord+1, collapse = ",")
  eval(parse(text = paste0("a[", idx, "] <- a[", idx, "]+1")))
  if(sum(coord) == 0) {
    countMono <- countMono+1
    print(paste0(countMono, " monomorphic sites"))
  } else if(sum(coord) == 1) {
    countSing <- countSing+1
    print(paste0(countSing, " singletons"))
    
  }
}
### Retrieving the multi-dimensional SFS with vector of indexes
entriesL <- lengths(entries)
n = length(entriesL)
ini = array(1, n)
fin = entriesL
#b <- array(1:prod(entriesL), dim = entriesL)

v = ini
iterN = prod(fin - ini + 1)
mSFS <- vector(mode = "numeric", length = iterN)

for(iter in 1:iterN) {
  #print(b[v]);
  idx <- paste(v, collapse = ",")
  print(paste(v-1, collapse = ","))
  mSFS[iter] <- eval(parse(text = paste0("a[", idx, "]")))
  #eval(parse(text = paste0("b[", idx, "]")))
  #print(eval(parse(text = paste0("b[", idx, "]"))))
  for(k in n:1) {
    if(v[k] < fin[k]) {
      #print(k);
      #print(v[k]);
      
      v[k] = v[k] + 1;
      break; 
    }
    v[k] = ini[k];
  }
}
cat("1 observations. No. of demes and sample sizes are on next line", file=paste0(sfsDir, "/", model, "_DSFS.tmp"), sep = "\n")
cat(paste(length(cladeNames), paste(lengths(entries)-1, collapse = "\t"), sep = "\t"), file=paste0(sfsDir, "/", model, "_DSFS.tmp"), sep = "\n", append = T)
cat(paste(mSFS, collapse = "\t"), file=paste0(sfsDir, "/", model, "_DSFS.tmp"), sep = "\n", append = T)


