#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
########################
########################
##
##
## This R script requires: 
##
##  Isabel Alves - March 2024
##
########################
########################

#library(tidyverse, lib.loc = "/shared/home/ialves/R/x86_64-redhat-linux-gnu-library/4.3")
library(tidyverse)
projDir <-"/shared/home/ialves/demoHist_yeast3039"
AAmatrix <- "03-data/SNPs_all.csv"
samplesFile <- "04-analysis/ABBA-BABA/samplesToKeep.tmp"
listOfStrains <- "04-analysis/ABBA-BABA/strains_of_cladeList.txt"
vcfName <- "04-analysis/ABBA-BABA/AW_AD_ED.mac1.Dsuite"
# --------------
# opening tables
# --------------
GT_pos <- read.table(file=paste0(projDir,"/", vcfName, ".chrpos"), header = F, sep = "\t")
colnames(GT_pos) <- c("Chromosome", "Position")
GT_pos <- GT_pos %>% mutate(SNVid=paste0(Chromosome, "_", Position))
AAallSNVs <- read.csv(file = paste0(projDir, "/", AAmatrix))

GT_pos <- left_joint(GT_pos,AAallSNVs,by=SNVid)
GT_pos <- GT_pos %>% filter(!is.na(Chromosome.x))

alt_pos <-  GT_pos %>% filter(AlignAllele=="ALT") %>% select(Chromosome, Position.x)
ref_pos <-  GT_pos %>% filter(AlignAllele=="REF") %>% select(Chromosome, Position.x)

write.table(alt_pos, file = paste0(projDir, "/", vcfName, ".allSNVs.alt"), col.names = F,
            row.names = F, quote = F, sep = "\t")
write.table(ref_pos, file = paste0(projDir, "/", vcfName, ".allSNVs.ref"), col.names = F,
            row.names = F, quote = F, sep = "\t")