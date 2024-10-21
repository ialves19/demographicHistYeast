wrkDir <- "/Users/isabel/Dropbox/UnivSTRASBOURG/PROJECTS/demoHist_yeast3039/04-analysis/ABBA-BABA/ABBA-BABA_1279strains"

fName_BBAA <- "ABBA-BABA_BBAA.txt"
openDF_BBAA <- read.table(paste0(wrkDir, "/", fName_BBAA), header = T, sep = "\t")

fName_Dmin <- "ABBA-BABA_Dmin.txt"
openDF_Dmin <- read.table(paste0(wrkDir, "/", fName_Dmin), header = T, sep = "\t")

colorCodeDf <- read.csv(file="/Users/isabel/Dropbox/UnivSTRASBOURG/PROJECTS/demoHist_yeast3039/05-results/Victor/CladesColors.csv", header = T)


openDF_BBAA <- openDF_BBAA %>% mutate(P1P2=paste0(P1,", ",P2), standError=Dstatistic/Z.score, Name=P3)
openDF_Dmin <- openDF_Dmin %>% mutate(P1P2=paste0(P1,", ",P2), standError=Dstatistic/Z.score, Name=P3)

mergedDF_BBAA <- left_join(openDF_BBAA, colorCodeDf, by="Name")
mergedDF_BBAA <- mergedDF_BBAA %>% mutate(adjPval=p.adjust(p.value,method="BH"))

mergedDF_BBAA %>% filter(adjPval < 0.01) %>% count()
mergedDF_BBAA %>% count()

mergedDF_BBAA <- mergedDF_BBAA %>% filter(!is.na(P2), !is.na(P3))

# order Clades according to clade number
# ---> !! closer numbers are closer in the NJ tree !! <---
orderClades <- list()
orderClades[["Name"]] <- sort(unique(mergedDF_BBAA$P3))
x <- strsplit(sort(unique(mergedDF_BBAA$P3)), split = "\\.")
x <- as.numeric(lapply(x, `[[`, 1))
orderClades[["Number"]] <- as.numeric(x)
df_orderClades <- as.data.frame(do.call(cbind, orderClades))
df_orderClades$Number <- as.numeric(df_orderClades$Number) 

df_orderClades[order(df_orderClades$Number),]
#####---------------------
#####----------
#####----

cladesOrder <- df_orderClades[order(df_orderClades$Number),1]

colorCode <- colorCodeDf$Color[match(factor(cladesOrder),colorCodeDf$Name)]
mergedDF_BBAA$P3 <- factor(mergedDF_BBAA$P3, levels = cladesOrder, ordered = T)

# plot D stats heatmap
# keep only the maximum values 
maxBBAA_df <-  mergedDF_BBAA %>% group_by(P3, P2) %>% top_n(1, Dstatistic) 

Dstats_m <- as.data.frame(maxBBAA_df$Dstatistic, ncol=1)
colnames(Dstats_m) <- "Dstats"
Dstats_m <- Dstats_m %>% mutate(P3=maxBBAA_df$P3, P2=maxBBAA_df$P2) 
Dstats_m <- Dstats_m %>% pivot_wider(names_from = P3, values_from = Dstats)
Dstats_m <- Dstats_m %>% filter(!is.na(P2))
mDstats <- t(as.matrix(Dstats_m[,-1]))
colnames(mDstats) <- Dstats_m$P2
mDstats <- mDstats[match(cladesOrder,rownames(mDstats)),match(cladesOrder,colnames(mDstats))]
colnames(mDstats) <- cladesOrder
mDstats <-  mDstats[,!colnames(mDstats) == "24. Lab Strains"]
#mDstats[upper.tri(mDstats)] <- t(mDstats)[upper.tri(mDstats)]
#{
my_sample_name <- colorCodeDf %>% filter(Level == "Clade", Name %in% rownames(mDstats)) %>% select(Name)
colnames(my_sample_name) <- "Clade"
my_sample_col <- colorCodeDf %>% filter(Level == "Clade", Name %in% rownames(mDstats)) %>% select(Color)
row.names(my_sample_name) <- unlist(colorCodeDf %>% filter(Level == "Clade", Name %in% rownames(mDstats)) %>% select(Name))
tmp <- my_sample_col$Color
names(tmp) <- my_sample_name$Clade
my_color <- list(Clade=tmp)
p <- pheatmap(t(mDstats), cutree_rows = 4, border_color = "white",
         annotation_colors = my_color, annotation_col = my_sample_name, annotation_legend = F)
#}
#my_hclust_gene <- hclust(dist(mDstats), method = "complete")
###
## p-value matrix
##
pVal_m <- as.data.frame(maxBBAA_df$adjPval, ncol=1)
colnames(pVal_m) <- "Adj.Pvalue"
pVal_m <- pVal_m %>% mutate(P3=maxBBAA_df$P3, P2=maxBBAA_df$P2) 
pVal_m <- pVal_m %>% pivot_wider(names_from = P3, values_from = Adj.Pvalue)
pVal_m <- pVal_m %>% filter(!is.na(P2))
mPValue <- t(as.matrix(pVal_m[,-1]))
colnames(mPValue) <- pVal_m$P2
mPValue <- mPValue[match(cladesOrder,rownames(mPValue)),match(cladesOrder,colnames(mPValue))]
colnames(mPValue) <- cladesOrder
mPValue <-  mPValue[,!colnames(mPValue) == "24. Lab Strains"]

#mPValue[upper.tri(mPValue)] <- t(mPValue)[upper.tri(mPValue)]

alpha <- c(0.01,0.001,0.0001)
mPValueSig <- mPValue
mPValueSig[mPValue > alpha[1]] <- "ns"
mPValueSig[mPValue < alpha[1] & mPValue > alpha[2]] <- "*"
mPValueSig[mPValue < alpha[2] & mPValue > alpha[3]] <- "**"
mPValueSig[mPValue < alpha[3]] <- "***"

pheatmap(t(mDstats), cutree_rows = 4,cutree_cols = 4, border_color = "white",
         display_numbers = t(mPValueSig), fontsize_number = 10, 
         annotation_colors = my_color, annotation_col = my_sample_name, annotation_legend = F)

#####
##
## F4-ratio
##
f4ratio_m <- as.data.frame(maxBBAA_df$f4.ratio, ncol=1)
colnames(f4ratio_m) <- "Dstats"
f4ratio_m <- f4ratio_m %>% mutate(P3=maxBBAA_df$P3, P2=maxBBAA_df$P2) 
f4ratio_m <- f4ratio_m %>% pivot_wider(names_from = P3, values_from = Dstats)
f4ratio_m <- f4ratio_m %>% filter(!is.na(P2))
mf4ratio <- t(as.matrix(f4ratio_m[,-1]))
colnames(mf4ratio) <- f4ratio_m$P2
mf4ratio <- mf4ratio[match(cladesOrder,rownames(mf4ratio)),match(cladesOrder,colnames(mf4ratio))]
colnames(mf4ratio) <- cladesOrder

mf4ratio[upper.tri(mf4ratio)] <- t(mf4ratio)[upper.tri(mf4ratio)]
pheatmap(mf4ratio, cutree_rows = 4, cutree_cols = 4, border_color = "white",
         annotation_colors = my_color, annotation_col = my_sample_name, annotation_legend = F)
