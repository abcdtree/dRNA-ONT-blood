.libPaths(libdir)

"%&%" <- function(a,b) paste(a,b, sep = "")


library(ggplot2)
library(tidyverse)
library(viridis)
library(reshape)

get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

resdir_sel="/data"



outf_pearson = paste(resdir_sel, "/corr_pearson.txt", sep="/")
outf_spearman = paste(resdir_sel, "/corr_spearman.txt", sep="/")

corrp_df=read.table(outf_pearson,head=T)
rownames(corrp_df)=corrp_df[,1]

rownames(corrp_df) <- c("bac1_illumina","bac1_ONT","bac2_illumina","bac2_ONT","bac3_illumina","bac3_ONT",
                        "vir1_illumina","vir1_ONT","vir4_illumina","vir4_ONT","bac4_illumina","bac4_ONT","bac5_illumina","bac5_ONT",
                        "bac6_illumina","bac6_ONT","vir5_illumina","vir5_ONT","vir2_illumina","vir2_ONT","vir3_illumina","vir3_ONT",
                        "vir6_illumina","vir6_ONT")
#colnames(corrp_df)
corrp_m=as.matrix(corrp_df[,-1])
colnames(corrp_m) <- rownames(corrp_df)

corrp_m


corrp_m = corrp_m[order(row.names(corrp_m)),order(colnames(corrp_m))]
corrp_m

type="Pearson"
corr_m=corrp_m

upper_tri <- get_upper_tri(corr_m)
melted_cormat <- melt(upper_tri, na.rm = TRUE)

gpl<-ggplot(data = melted_cormat, aes(X2, X1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", na.value="white",
                       midpoint = 0, limit = c(0,1), space = "Lab", 
                       name=type) +
  theme_minimal()+ 
  theme(panel.background = element_blank(),axis.text.x = element_text(angle = 45, vjust = 1, 
                                                                      size = 6, hjust = 1),axis.text.y = element_text(vjust = 1, 
                                                                                                                      size = 6, hjust = 1))+
  coord_fixed()

gpl<-gpl + 
  geom_text(aes(X2, X1, label = value), color = "black", size = 2) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))

gpl
pdf(work_dir%&%"/Peason_Corr_heatmap-rescale-reorder-transcript-level.pdf")
gpl
dev.off()

corrsp_df=read.table(outf_spearman,head=T)
rownames(corrsp_df)=corrsp_df[,1]
corrsp_m=as.matrix(corrsp_df[,-1])

type="Spearman"
corr_m=corrsp_m

upper_tri <- get_upper_tri(corr_m)
melted_cormat <- melt(upper_tri, na.rm = TRUE)

gpl<-ggplot(data = melted_cormat, aes(X2, X1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", na.value="white",
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name=type) +
  theme_minimal()+ 
  theme(panel.background = element_blank(),axis.text.x = element_text(angle = 45, vjust = 1, 
                                                                      size = 6, hjust = 1),axis.text.y = element_text(vjust = 1, 
                                                                                                                      size = 6, hjust = 1))+
  coord_fixed()

gpl<-gpl + 
  geom_text(aes(X2, X1, label = value), color = "black", size = 2) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))

gpl

pdf(resdir_sel%&%"/Spearman_Corr_heatmap.pdf")
gpl
dev.off()
