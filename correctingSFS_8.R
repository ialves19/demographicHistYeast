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

# arguments 
inDir <- args[1]
sites <- args[2]
modelFile <- args[3]

# inDir <- "/shared/home/ialves/demoHist_yeast3039/04-analysis/fastsimcoal/fiveClades"
# sites <- "sites27samples.txt"
# modelFile <- "model_OneDomest.txt"
# ----------------
# opening FILES
# ----------------
model <- scan(paste0(inDir, "/", modelFile), nlines = 1, what = character())
struModel <- scan(paste0(inDir, "/", modelFile), skip = 1, nlines = 1, what = numeric())
df_zig <- read.table(paste0(inDir, "/", modelFile), skip = 2, header = F, sep = "\t")

# ----------------
# creating model specific folder
# ----------------
sfsDir <- paste0(inDir, "/", model, "_dSFS")
runScript <- 0
if(dir.exists(sfsDir)){
  print(paste0("Folder ", sfsDir, " already exists"))
  print(paste0("Correcting SFS in: ", sfsDir))
  runScript <- 1
} else {
  print(paste0("ERROR: Folder ", sfsDir, " doesn't exist"))
}


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
# opening sites info file
# ----------------
sitesList <- list()
sitesTmp <- scan(file=paste0(sfsDir, "/", sites), what = character())
sitesList$mono <- as.numeric(sitesTmp[2])
sitesList$poly <- as.numeric(sitesTmp[4])
sitesList$polyFinal <- as.numeric(sitesTmp[6])
# ----------------
# computing total number of monomorphic sites
# ----------------
propSNPsFinal <- sitesList$polyFinal/sitesList$poly
nbSitesAA_NCD_NoMiss_approx <- sitesList$mono*propSNPsFinal
print(paste0("Number of sites surveyed: ", nbSitesAA_NCD_NoMiss_approx))

if(runScript) {
  # ----------------
  # correcting 2D-SFS
  # ----------------
  pairsClades <- combn(length(cladeNames),2, replace=F)
  for(p in 1:ncol(pairsClades)) {
    
    #p <- 1
    #cladePair <- 1
    cladeOne_name <-  cladeNames[pairsClades[1,p]] 
    cladeTwo_name <-  cladeNames[pairsClades[2,p]]
    cat("Correcting mono 2D-SFS for the pair:", cladeOne_name, cladeTwo_name, sep = "\n")
    
    m_2DSFS <- read.table(paste0(sfsDir, "/", model, "_jointDAFpop", which(names(popIDs) == cladeTwo_name)-1, "_", which(names(popIDs) == cladeOne_name)-1, ".tmp"), header = T, row.names = 1, skip = 1)
    
    finalNbSNPs_sfs <- sum(m_2DSFS)-m_2DSFS[1,1] #here we lost some variable sites in the homozygous strains 
    print("Number of SNPs in the sfs:")
    finalNbSNPs_sfs
    
    
    nbMono_sfs <- nbSitesAA_NCD_NoMiss_approx-finalNbSNPs_sfs+m_2DSFS[1,1]
    print(paste0("Number of monomorphic sites: ", nbMono_sfs))
    
    m_2DSFS[1,1] <- ceiling(nbMono_sfs)
    
    cat("1 observation", file = paste0(sfsDir, "/", model, "_","jointDAFpop", which(names(popIDs) == cladeTwo_name)-1, "_", which(names(popIDs) == cladeOne_name)-1,".obs"), sep = "\n")
    write.table(m_2DSFS, file=paste0(sfsDir, "/", model, "_","jointDAFpop", which(names(popIDs) == cladeTwo_name)-1, "_", which(names(popIDs) == cladeOne_name)-1,".obs"), row.names = TRUE, col.names = NA, quote = F, sep = "\t", append = T)
    
  }
  # ------------------------------
  # ------------------
  # -----------
  # ----------------
  # opening multi-SFSs
  # ----------------
  headermultiSFS <- scan(file=paste0(sfsDir, "/", model, "_DSFS.tmp"), skip = 1, nlines = 1)
  multiSFS <- scan(file=paste0(sfsDir, "/", model, "_DSFS.tmp"), skip = 2, nlines = 1)
  # ----------------
  # correcting multi-SFSs
  # ----------------
  finalNbSNPs_sfs <- sum(multiSFS[-1]) #here we lost some variable sites in the homozygous strains 
  print(paste0("Number of SNPs in the sfs: ", finalNbSNPs_sfs))
  
  nbMono_sfs <- nbSitesAA_NCD_NoMiss_approx-finalNbSNPs_sfs+multiSFS[1]
  print(paste0("Number of monomorphic sites: ", nbMono_sfs))
  
  multiSFS[1] <- ceiling(nbMono_sfs)
  # ----------------
  # exporting multi-SFSs
  # ----------------
  cat("1 observations. No. of demes and sample sizes are on next line", file=paste0(sfsDir, "/", model, "_DSFS.obs"), sep = "\n")
  cat(paste(headermultiSFS, collapse = "\t"), file=paste0(sfsDir, "/", model, "_DSFS.obs"), sep = "\n", append = T)
  cat(paste(multiSFS, collapse = "\t"), file=paste0(sfsDir, "/", model, "_DSFS.obs"), sep = "\n", append = T)
}