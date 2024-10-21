#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#myPaths <- .libPaths()
#libs <- c("/shared/home/ialves/R/x86_64-redhat-linux-gnu-library/4.3", .libPaths())
#.libPaths(new = libs)
#myPaths <- c(myPaths, "/shared/home/ialves/R/x86_64-redhat-linux-gnu-library/4.3/tidyverse")
#.libPaths("/shared/home/ialves/R/x86_64-redhat-linux-gnu-library/4.3") 

#library(tidyverse, lib.loc = "/shared/home/ialves/R/x86_64-redhat-linux-gnu-library/4.3")
library(tidyverse)
typeSNPs <- args[1]
#typeSNPs <- "onlyNonCoding" or "all"


wrkDir <- "/shared/home/ialves/demoHist_yeast3039"
# opening file with AA information
AAstateFileName <- "/shared/home/ialves/ancestralState/03-Data/full3039Matrix.AllPositions.DP10.GQ20.Mind20.99pNonMiss.ExcHet99.chromosome"
AA <- data.frame()

#opening files with AA 
for(chr in 1:16) {
  
  openFile <- read.table(file=paste0(AAstateFileName, chr, ".JubinExtracted.SNmono.BiSNPs.AA"), 
                         header = T, sep = "\t")
  AA <- rbind(AA, openFile %>% filter(AA != "?"))
  
}
AA <- AA %>% mutate(SNVid=paste0(Chromosome, "_", Position))

# opening files with annotation information

#open file with annotation
# here we are changing also the ChrID which in the original ann matrix is chrI chrII and so on
openAnno <- read.table(file=paste0(wrkDir, "/03-data/annotation/full3039Matrix.DP10.GQ20.Mind20.99pNonMiss.ExcHet99.ann_4SIFT_SIFTannotations.xls"), sep = "\t", header = T)

romanNb <- as.roman(1:16)
chrIDnew <- rep(NA, nrow(openAnno))


for(chr in 1:16) {
  chrIDnew[which(openAnno$CHROM == paste0("chr", romanNb[chr]))] <- paste0("chromosome", chr)
}

sum(is.na(chrIDnew))
openAnno <- openAnno %>% mutate(Chromosome=chrIDnew)
openAnno <- openAnno %>% mutate(SNVid=paste0(Chromosome, "_", POS))

GT_AA_anno_m <- left_join(AA, openAnno, by="SNVid")

if(typeSNPs == "all") {
  
  write.csv(GT_AA_anno_m, file = paste0(wrkDir, "/03-data/SNPs_", typeSNPs, ".csv"), row.names = F, col.names = T)
} else if(typeSNPs == "onlyNonCoding") {
  SNVs_AA_notFunctional <- GT_AA_anno_m %>% filter(is.na(POS))
  write.csv(SNVs_AA_notFunctional, file = paste0(wrkDir, "/03-data/SNPs_", typeSNPs, ".csv"), row.names = F, col.names = T)
}

