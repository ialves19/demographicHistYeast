#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
########################
########################
##
##
## This R script requires: 
##  1. genotype, chrPos and sample files to be in the same dir
##  2. SNP matrix and operational table (strains info) to be in the same dir
##
##  This script is called using: 
##  sbatch convertREFtoAA_5.sh -g /shared/home/ialves/demoHist_yeast3039/04-analysis/fastsimcoal/SFS_134samples/allChrom.1SNP.25samples.noORF.noMiss.SNVs.demoHist \
##  -s samples25.txt -l list_of_strains_per_Clade.txt -t "onlyNonCoding" -d /shared/home/ialves/demoHist_yeast3039/03-data
##
##  Isabel Alves - March 2024
##
########################
########################

#library(tidyverse, lib.loc = "/shared/home/ialves/R/x86_64-redhat-linux-gnu-library/4.3")
library(tidyverse)
GTMatrix <- args[1]
samplesFile <- args[2]
listOfStrains <- args[3]
typeSNPs <- args[4]
pathToData <- args[5]
  
tmpDir <- unlist(strsplit(GTMatrix, split="/"))
pathToAnalysis <- paste(tmpDir[-length(tmpDir)], collapse="/")
rm(tmpDir)
# projDir <- "/shared/home/ialves/demoHist_yeast3039"
# GTMatrix <- "SFS_137samples/demoInf_27Strains_SNVs"
# samplesFile <- "SFS_137samples/samples27.txt"
# listOfStrains <- "SFS_137samples/list_of_strains_per_Clade.txt"
# --------------
# opening tables
# --------------
GT_pos <- read.table(file=paste0(GTMatrix, ".chrpos"), header = F, sep = "\t")
GT_m <- read.table(file=paste0(GTMatrix, ".coded.GT"), header = F, sep = "\t")
SNVs_AA_notFunctional <- read.csv(file = paste0(pathToData, "/SNPs_", typeSNPs, ".csv"))
openCladeInfo <- read.csv(file=paste0(pathToData, "/operationalTable_Full3039Sace_Clades_Ploidy_Aneuploidy.csv"))
# --------------
# opening samples
# --------------
samples <- scan(file=paste0(pathToAnalysis,"/",samplesFile), what = character())
# give names to the genotypes
colnames(GT_m) <- samples

# correct matrix if additional columns are present
GT_m <- GT_m %>% select(all_of(samples))
# --------------
# creating a commun column between matrices
# --------------
colnames(GT_pos) <- c("Chromosome", "Position")
GT_pos <- GT_pos %>% mutate(SNVid=paste0(Chromosome, "_", Position))

indexALTSNVs <- which(GT_pos$SNVid %in% SNVs_AA_notFunctional$SNVid)
# --------------
# filter GT matrices to keep only the list of SNPs we want 
# --------------
GT_pos <- GT_pos[indexALTSNVs,]
GT_m <- GT_m[indexALTSNVs,]

# --------------
# get info about the strains
# --------------
samplesPerPop <- read.table(file=paste0(pathToAnalysis, "/", listOfStrains), header = T, sep = "\t")
CladeIDs <- names(table(samplesPerPop$CladeID))
samplesPerPop <- samplesPerPop %>% mutate(StandardizedName=Strains)
samplesPerPop <- left_join(samplesPerPop, openCladeInfo, by="StandardizedName")

# excluding aneuploids and keeping only the most common zygosity strains
listSamplesPerCladeIDsClean <- list()
zigosityStatus <- list()

for(clade in CladeIDs) {
  #clade <- CladeIDs[2]
  zigosityStatus[[clade]] <- names(table(samplesPerPop %>% filter(CladeID==clade, Ploidy=="2", Aneuploidy=="euploid") %>% select(Zygosity)))[table(samplesPerPop %>% filter(CladeID==clade, Ploidy=="2", Aneuploidy=="euploid") %>% select(Zygosity)) == max(table(samplesPerPop %>% filter(CladeID==clade, Ploidy=="2", Aneuploidy=="euploid") %>% select(Zygosity)))]
  if(clade == "6. West Europe Wine") {
    zigosityStatus[[clade]] <- "Heterozygous"
  }
  listSamplesPerCladeIDsClean[[clade]] <- as.character(unlist(samplesPerPop %>% filter(CladeID==clade, Ploidy=="2", Aneuploidy=="euploid", Zygosity == zigosityStatus[[clade]]) %>% select(StandardizedName)))
}
# --------------
# # transforming allele counts according to AA
# --------------
SNVs_subset <- left_join(GT_pos, SNVs_AA_notFunctional, by="SNVid")
altSNVs <- which(SNVs_subset$AlignAllele == "ALT")
GT_copy <- GT_m

for(k in altSNVs) {
  GT_copy[k,which(GT_copy[k,] == 2)] <- 4
  GT_copy[k,which(GT_copy[k,] == 0)] <- 2
  GT_copy[k,which(GT_copy[k,] == 4)] <- 0
  
}
# --------------
# saving files
# --------------
# transforming list Samples per Clade into df
df_listSamples <- do.call(rbind,(lapply(listSamplesPerCladeIDsClean, function(x) {data.frame(x)})))
colnames(df_listSamples) <- c("StandardizedName")

# saving samples per clade
write.table(df_listSamples, file=paste0(pathToAnalysis, "/df_Clade_Strains.txt"), sep = "\t", col.names = F, row.names = T, quote = F)
# transforming list with zygosity status into a df
df_zigozity <- do.call(rbind, zigosityStatus)
# saving zygosity status
write.table(df_zigozity, file=paste0(pathToAnalysis, "/df_zygosity.txt"), col.names = F, row.names = T, quote = F, sep = "\t")
# saving GT matrix
nbSamples <- dim(GT_copy)[2]
write.table(GT_copy, file=paste0(pathToAnalysis, "/samples", nbSamples, ".coded.AA.GT"), col.names = F, row.names = F, quote = F, sep = "\t")
# saving matrix pos
tmp_m <- SNVs_subset %>% select(all_of(c("Chromosome","Position.x","REF", "ALT", "AlignAllele", "AlignStatus", "AA")))
write.table(tmp_m, file=paste0(pathToAnalysis, "/samples", nbSamples, ".AA.chrPos"), col.names = F, row.names = F, quote = F, sep = "\t")



