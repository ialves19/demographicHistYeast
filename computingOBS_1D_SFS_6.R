#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

########################
########################
##
##
## This R script requires: 
##  1. the genotype matrix, 
##  2. the sample file, 
##  3. the table with strains/clade, and
##  4. the table with the zygosity 
##
##  to be in the same dir: wrkDir
##  This script is called using: computingOBS_1D_SFS_6.R
##  the seed number is related to the randomly drawn allele at heterozygous
##  sites in the homozygous strains.
##  In case you want to treat the homozygous strains as heterozygous, 
##  change the the table with the zygosity accordingly. 
##
##  Isabel Alves - March 2024
##
########################
########################

#library(tidyverse, lib.loc = "/shared/home/ialves/R/x86_64-redhat-linux-gnu-library/4.3")
library(tidyverse)
#library(ggpubr)
library(ggplot2)
set.seed(988599694)

# arguments 
wrkDir <- args[1]
samplesFile <- args[2]
dfStrains <- args[3]
dfZygosity <- args[4]
GTmName <- args[5]
boot <- args[6]

if(boot == F) {
  outDir <- wrkDir
} else {
  bootDir <- unlist(strsplit(GTmName, split = "/"))[1]
  outDir <- paste0(wrkDir, "/", bootDir)
}

# wrkDir <- "/shared/home/ialves/demoHist_yeast3039/04-analysis/fastsimcoal/SFS_137samples"
# samplesFile <- "samples27.txt"
# dfStrains <- "df_Clade_Strains.txt"
# dfZygosity <- "df_zygosity.txt"
# GTmName <- "samples27.coded.AA.GT"
# ----------------
# opening files
# ----------------
df_clades <- read.table(paste0(wrkDir, "/", dfStrains), header = F, sep = "\t")
df_zig <- read.table(paste0(wrkDir, "/", dfZygosity), header = F, sep = "\t")
GTmatrix <- read.table(paste0(wrkDir, "/", GTmName), header = F, sep = "\t")
samples <- scan(paste0(wrkDir, "/", samplesFile), what = character())
colnames(GTmatrix) <- samples
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

# ----------------
# computing 1D SFS
# ----------------
CladeIDs <- names(zigosityStatus)
sfs_list <- list()
dsfs_1D <- list()
# compute 1D-dSFS
for(clade in CladeIDs) {
  
  #clade <- CladeIDs[5]
  strains <- unlist(listSamplesPerCladeIDsClean[[clade]])
  sampleSet <- colnames(GT_copy[,colnames(GT_copy) %in% strains])
  subGTm <- GT_copy[,colnames(GT_copy) %in% strains]
  if(zigosityStatus[[clade]] == "Homozygous") {
    subGTm[subGTm == 1] <- sample(c(0,2), 1)  
    subGTm[subGTm == 2] <- 1
    
    entries <- rep(0, length(sampleSet)+1)
    names(entries) <- as.character(0:length(sampleSet))
    entries[names(entries) %in% names(table(apply(subGTm, 1, sum)))] <- table(apply(subGTm, 1, sum))
    sfs_list[[clade]] <- entries
    write.table(file=paste0(outDir, "/1D-SFS_", clade, ".sfs"), as.matrix(sfs_list[[clade]]), col.names = F, row.names = T, quote = F)

  } else {
    entries <- rep(0, (length(sampleSet)*2)+1)
    names(entries) <- as.character(0:(length(sampleSet)*2))
    entries[names(entries) %in% names(table(apply(subGTm, 1, sum)))] <- table(apply(subGTm, 1, sum))
    sfs_list[[clade]] <- entries
    write.table(file=paste0(outDir, "/1D-SFS_", clade, ".sfs"), as.matrix(sfs_list[[clade]]), col.names = F, row.names = T, quote = F)
  }
}
