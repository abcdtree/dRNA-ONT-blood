.libPaths(libdir)

"%&%" <- function(a,b) paste(a,b, sep = "")


library(ggplot2)
library(tidyverse)
library(viridis)
library(reshape)
library(RColorBrewer)

data <- read.table("/data/category_num.txt", header=TRUE)
coul <- brewer.pal(5, "Set2") 

data$group <- factor(data$group, levels=c("novel_not_in_catalog", "novel_in_catalog", "incomplete_splice_match","full_splice_match","intergenic"))
ggplot(data, aes(x=group, y=number, fill=as.factor(group))) + 
    geom_bar(stat = "identity", width=0.6) +
    scale_fill_manual(values = coul) +
    scale_y_continuous(breaks=seq(0,300,100),limits = c(0, 300))+
    theme_minimal()+ 
    theme(legend.position="none",axis.text.x = element_text(angle = 10, vjust = 0.8, 
                                                                      size = 10, hjust = 0.5),axis.text.y = element_text(vjust = 1, , 
                                                                                                                      size = 10, hjust = 1))+
                                                                                                                      
    xlab("Structural category")+
    ylab("Number of novel isoforms")

data <- read.table("/data/category_num.txt", header=TRUE)
coul <- brewer.pal(5, "Set2") 

data$group <- factor(data$group, levels=c("novel_not_in_catalog", "novel_in_catalog", "incomplete_splice_match","full_splice_match","intergenic"))
ggplot(data, aes(x=group, y=number, fill=as.factor(group))) + 
  geom_bar(stat = "identity", width=0.6) +
  scale_fill_manual(values = coul) +
  theme_minimal()+ 
  scale_y_continuous(breaks=seq(0,600,200),limits = c(0, 600))+
  theme(legend.position="none",axis.text.x = element_text(angle = 10, vjust = 0.8, 
                                                          size = 10, hjust = 0.5),axis.text.y = element_text(vjust = 1, , 
                                                                                                           size = 10, hjust = 1))+
  
  xlab("Structural category")+
  ylab("Number of novel isoforms")
  

data <- read.table("/data/exons_num.txt", header=TRUE)
sorted_data  <- data [order(data$group), ]
barplot(height=sorted_data$number, names=sorted_data$group, 
        col=rgb(0.8,0.1,0.1,0.6),
        xlab="Number of exons", 
        ylab="Number of novel isoforms", 
        ylim=c(0,90),
        names.arg=sorted_data$group,cex.names = 1.0
)


data <- read.table("/data/exons_num.txt", header=TRUE)
sorted_data  <- data [order(data$group), ]
barplot(height=sorted_data$number, names=sorted_data$group, 
        col=rgb(0.8,0.1,0.1,0.6),
        xlab="Number of exons", 
        ylab="Number of novel isoforms", 
        ylim=c(0,200),
        names.arg=sorted_data$group,cex.names = 1.0
)

data <- read.table("/data/chrom_num.txt", header=TRUE)
data$group <- factor(data$group, levels=c("chr1", "chr2", "chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"))
sorted_data  <- data [data$group, ]

barplot(height=sorted_data$number, names=sorted_data$group, 
        col="#69b3a2",
        xlab="chromosome", 
        ylab="Number of novel isoforms", 
        ylim=c(0,80),
        names.arg=sorted_data$group,cex.names = 1.0
)


data <- read.table("/data/chrom_num.txt", header=TRUE)
data$group <- factor(data$group, levels=c("chr1", "chr2", "chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"))
sorted_data  <- data [data$group, ]

barplot(height=sorted_data$number, names=sorted_data$group, 
        col="#69b3a2",
        xlab="chromosome", 
        ylab="Number of novel isoforms", 
        ylim=c(0,100),
        names.arg=sorted_data$group,cex.names =1.0
)



data <- read.table("/data/category_translength.txt", header=TRUE)
coul <- brewer.pal(5, "Set2") 
data$group <- factor(data$group, levels=c("novel_not_in_catalog", "novel_in_catalog", "incomplete_splice_match","full_splice_match","intergenic"))


ggplot(data,aes(x=group, y=number, fill=group)) +
  geom_boxplot() +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  scale_fill_manual(values = coul) +
  theme_minimal()+ 
  scale_y_continuous(breaks=seq(0,10000,2000),limits = c(0, 10000))+
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    axis.text.x = element_text(angle = 10, vjust = 0.8, 
                                                      size = 10, hjust = 0.5),axis.text.y = element_text(vjust = 1, , 
                                                                                                         size = 10, hjust = 1)
  ) +
  #ggtitle("Basic boxplot") +
  xlab("Structural category")+
  ylab("Inferred transcript length")

data <- read.table("/data/category_translength.txt", header=TRUE)
coul <- brewer.pal(5, "Set2") 
data$group <- factor(data$group, levels=c("novel_not_in_catalog", "novel_in_catalog", "incomplete_splice_match","full_splice_match","intergenic"))


ggplot(data,aes(x=group, y=number, fill=group)) +
  geom_boxplot() +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  scale_fill_manual(values = coul) +
  theme_minimal()+ 
  scale_y_continuous(breaks=seq(0,10000,2000),limits = c(0, 10000))+
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    axis.text.x = element_text(angle = 10, vjust = 0.8, 
                               size = 10, hjust = 0.5),axis.text.y = element_text(vjust = 1, , 
                                                                                  size = 10, hjust = 1)
  ) +
  #ggtitle("Basic boxplot") +
  xlab("Structural category")+
  ylab("Inferred transcript length")

