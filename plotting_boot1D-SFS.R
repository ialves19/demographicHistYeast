library(tidyverse)

# ----------------
# Variables
# ----------------
modelDir <- "/shared/home/ialves/demoHist_yeast3039/04-analysis/fastsimcoal/fiveClades"
modelFile <- "model_TwoDomest_EurOld_Taiw_stru.txt"
wrkDir <- "/shared/home/ialves/demoHist_yeast3039/04-analysis/fastsimcoal/SFS_134samples"
nbBootstraps <- 100
SizePlot <- 35
# ----------------
# open files
# ----------------
modelPop <- read.table(file=paste0(modelDir,"/", modelFile), header = F, skip = 2, sep = "\t")

# ----------------
# open boot sfs 
# ----------------
v_pop <- modelPop$V1
dsfs_pop <- list()

# ----------------
# getting 1D-sfs across pop
# ----------------
for(p in 1:length(v_pop)) {
  
  #p <- v_pop[1]
  #p <- 2
  sfs1d_list <- list()
  for(boot in 1:nbBootstraps) {
    
    sfs1d_list[[boot]] <- read.table(file=paste0(wrkDir, "/bootstrap_DS_", boot, "/1D-SFS_", v_pop[p], ".sfs"), 
                                     row.names = 1, header = F, sep = " ")
      
  }
  
  sfs1d.df <- do.call(cbind, sfs1d_list)
  colnames(sfs1d.df) <- paste0("boot", 1:nbBootstraps)
  
  sfs1d.df.long <-  sfs1d.df %>% mutate(sfsEntries=row.names(sfs1d.df)) %>% pivot_longer(-sfsEntries,names_to = "boot", values_to="SFS-1D") 
  sfs1d.df.long$sfsEntries <- as.numeric(sfs1d.df.long$sfsEntries)
  dsfs_pop[[v_pop[p]]] <- sfs1d.df.long %>% mutate(Clade=rep(v_pop[p], nrow(sfs1d.df.long)))
  #rm(sfs1d.df)
}

# ----------------
# plotting 
# ----------------
   sfs.plot <- do.call(rbind, dsfs_pop) %>%
     filter(sfsEntries!=0) %>% ggplot(aes(x=sfsEntries, y=`SFS-1D`, colour=boot)) + 
     geom_line() + facet_wrap(~Clade) + scale_color_viridis(discrete = TRUE) + theme_classic() + theme(legend.position="none") +
     xlab("dSFS-entries") + ylab("Counts")
ggsave(file=paste0(wrkDir, "/bootstrap_DS_", boot, "/1D-SFS_allClades.pdf"), sfs.plot,
       height = SizePlot*5, width = SizePlot*8, units = "mm")

