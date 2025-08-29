## Creating sample size table 
# The idea is to keep info on the 
# number of samples and sites used for 
# demographic history inference. 
library(kableExtra)
library(tidyverse)
library(ggplot2)
library(knitr)

SFSfolder <- "/Users/isabel/Dropbox/UnivSTRASBOURG/PROJECTS/demoHist_yeast3039/05-results"
obsSFS <- c("SCALLING_WildEast_haploid/OneDomest_dSFS/multi-SFS/OneDomest_DSFS.obs",
                      "FiveClades/CHINESEwild/OneDomest_dSFS/OneDomest_DSFS.obs",
                      "FiveClades/TAIWANESEwild/OneDomestTaiw_dSFS/OneDomestTaiw_DSFS.obs")
names(obsSFS) <- c("31samples", "27samples", "25samples")

sfsSize <- list()

for(obs in 1:length(obsSFS)){
  #obs <- obsSFS[1]
  fileName <- paste0(SFSfolder, "/", obsSFS[obs])
  tmpSamples <- read.table(file=fileName, skip = 1, nrows = 1)
  tmpObs <- read.table(file=fileName, skip = 2, nrows = 1)
  sfsSize[[names(obsSFS)[obs]]]$length <- length(tmpObs)
  sfsSize[[names(obsSFS)[obs]]]$nonzero <- sum(tmpObs != 0)
}

mSfsSize <- t(matrix(unlist(sfsSize), nrow = 2, byrow = F))
rownames(mSfsSize) <- c("31samples", "27samples", "25samples")
colnames(mSfsSize) <- c("nbSFSentries", "NonzeroEntries")

samplesAndData <-  data.frame(DataName=c("Chi+Taiw", "Chi", "Taiw"), 
                              Samples=c(31,27,25),
                              nbCallableSitesNoMiss=c(8923794,8972807,8995776), 
                              nbSNVsNoMiss=c(308717,290681,153234), nbNonCodingSNVs=c(48155,45085,23232))
samplesAndData <- cbind(samplesAndData, mSfsSize)

kable(samplesAndData, row.names = F, caption = "Dataset") %>% 
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = F)