
library(dplyr)
library(fgsea)
library(clusterProfiler)
library(data.table)
library(stringr)
#library(ggplot2)
library(topGO)
library(Rgraphviz)
library(pathview)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)

library(MASS)
library(pracma)

"%&%" <- function(a,b) paste(a,b, sep = "")


rhome1="/data/"
rapidshome = paste(rhome1,"rapids_rna_seq", sep="/");
datahome= rhome1%&%"humanRNA/"
mt = read.table(datahome%&%"chrM_gene_v35.txt", sep="\t", head=F, as.is = T)

mt_name_df=data.frame(strsplit(as.character(mt[,9]), "\\;"))
row_as_char <- as.character(mt_name_df[4, ])
split_result <- strsplit(row_as_char, "gene_name=")
mt_name=unlist(split_result)
mt_name=mt_name[mt_name!=""]

#sorted_genes_list <- gene_df[order(-gene_df$mean), ]$gene
gene_df=read.table("polyA_length_on_gene_overlap_mRNA_rmGlobin.txt",header=T)
gene_df_mit=gene_df[gene_df$gene%in%mt_name,]
gene_df_nuclear=gene_df[!gene_df$gene%in%mt_name,]


############# excluding mitochrondrial genes##############
virus_gene_list=gene_df_nuclear$virus_mean
bac_gene_list=gene_df_nuclear$bacteria_mean
full_gene_list = gene_df_nuclear$mean

outdir="/data/ALL_MIT/"
gene_list=full_gene_list
names(gene_list) = gene_df_nuclear$gene

file_path <- outdir%&%"gene_polyAlen.txt"
write.table(gene_list, file=file_path, row.names=TRUE, col.names=FALSE, quote=FALSE,sep="\t")

gene_list = sort(gene_list, decreasing = TRUE) ##From polyA length largest to smallest
gene_list = gene_list[!duplicated(names(gene_list))]
head(gene_list)

DEG.entrez_id = mapIds(x = org.Hs.eg.db,
                       keys = names(gene_list),
                       keytype = "SYMBOL",
                       column = "ENTREZID")
gene_list_entrezid=gene_list
names(gene_list_entrezid)=DEG.entrez_id

##https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/
gse_go <- gseGO(geneList=gene_list, 
                ont ="ALL", 
                keyType = "SYMBOL", 
                nPerm = 10000, 
                minGSSize = 3, 
                maxGSSize = 800, 
                pvalueCutoff = 0.05, 
                verbose = TRUE, 
                OrgDb = org.Hs.eg.db, 
                pAdjustMethod = "BH")

gse_df=as.data.frame(gse_go)
outfile=paste(outdir,"GO_ALL.csv",sep="")
write.csv(gse_df , outfile, row.names=FALSE, quote = FALSE)

outfile1 <-paste(outdir,"/gsea_GO_ALL.pdf",sep="")
pdf(file=outfile1,width = 10,height = 10)
dotplot(gse_go, showCategory=10, split=".sign") + facet_grid(.~.sign)
dev.off()

outfile1 <-paste(outdir,"/enrich_GO_ALL.pdf",sep="")
pdf(file=outfile1,width = 10,height = 10)
ridgeplot(gse_go) 
 #+ labs(x = "polyA tail lengths")
dev.off()

pdf(outdir%&%"GSEA_enrichment_GO_ALL.pdf") 
for(p in 1:50){
  print(gse_go$Description[p])
  print(gseaplot(gse_go, by = "all", title = gse_go$Description[p], geneSetID = p))
  #dev.off()
}
dev.off() 

outfile0=paste(outdir,"gse_go.rds",sep="")
saveRDS(gse_go,outfile0)


gse_goBP <- gseGO(geneList=gene_list, 
                  ont ="BP", 
                  keyType = "SYMBOL", 
                  nPerm = 10000, 
                  minGSSize = 3, 
                  maxGSSize = 800, 
                  pvalueCutoff = 0.05, 
                  verbose = TRUE, 
                  OrgDb = org.Hs.eg.db, 
                  pAdjustMethod = "BH")

gse_BPdf=as.data.frame(gse_goBP)
outfile=paste(outdir,"GO_BP.csv",sep="")
write.csv(gse_BPdf , outfile, row.names=FALSE, quote = FALSE)

outfile1 <-paste(outdir,"/gsea_GO_BP.pdf",sep="")
pdf(file=outfile1,width = 10,height = 10)
dotplot(gse_goBP, showCategory=10, split=".sign") + facet_grid(.~.sign)
dev.off()

outfile1 <-paste(outdir,"/enrich_GO_BP.pdf",sep="")
pdf(file=outfile1,width = 10,height = 10)
ridgeplot(gse_goBP) #+ labs(x = "enrichment distribution")
dev.off()

pdf(outdir%&%"GSEA_enrichment_GO_BP.pdf") 
for(p in 1:50){
  print(gse_goBP$Description[p])
  print(gseaplot(gse_goBP, by = "all", title = gse_goBP$Description[p], geneSetID = p))
  #dev.off()
}
dev.off() 

outfile0=paste(outdir,"gse_goBP.rds",sep="")
saveRDS(gse_goBP,outfile0)


gse_goCC <- gseGO(geneList=gene_list, 
                  ont ="CC", 
                  keyType = "SYMBOL", 
                  nPerm = 10000, 
                  minGSSize = 3, 
                  maxGSSize = 800, 
                  pvalueCutoff = 0.05, 
                  verbose = TRUE, 
                  OrgDb = org.Hs.eg.db, 
                  pAdjustMethod = "BH")

gse_CCdf=as.data.frame(gse_goCC)
outfile=paste(outdir,"GO_CC.csv",sep="")
write.csv(gse_CCdf , outfile, row.names=FALSE, quote = FALSE)

outfile1 <-paste(outdir,"/gsea_GO_CC.pdf",sep="")
pdf(file=outfile1,width = 10,height = 10)
dotplot(gse_goCC, showCategory=10, split=".sign") + facet_grid(.~.sign)
dev.off()

outfile1 <-paste(outdir,"/enrich_GO_CC.pdf",sep="")
pdf(file=outfile1,width = 10,height = 10)
ridgeplot(gse_goCC) #+ labs(x = "enrichment distribution")
dev.off()

pdf(outdir%&%"GSEA_enrichment_GO_CC.pdf") 
for(p in 1:50){
  print(gse_goCC$Description[p])
  print(gseaplot(gse_goCC, by = "all", title = gse_goCC$Description[p], geneSetID = p))
  #dev.off()
}
dev.off() 

outfile0=paste(outdir,"gse_goCC.rds",sep="")
saveRDS(gse_goCC,outfile0)


gse_goMF <- gseGO(geneList=gene_list, 
                  ont ="MF", 
                  keyType = "SYMBOL", 
                  nPerm = 10000, 
                  minGSSize = 3, 
                  maxGSSize = 800, 
                  pvalueCutoff = 0.05, 
                  verbose = TRUE, 
                  OrgDb = org.Hs.eg.db, 
                  pAdjustMethod = "BH")

gse_MFdf=as.data.frame(gse_goMF)
outfile=paste(outdir,"GO_MF.csv",sep="")
write.csv(gse_MFdf , outfile, row.names=FALSE, quote = FALSE)

outfile1 <-paste(outdir,"/gsea_GO_MF.pdf",sep="")
pdf(file=outfile1,width = 10,height = 10)
dotplot(gse_goMF, showCategory=10, split=".sign") + facet_grid(.~.sign)
dev.off()

outfile1 <-paste(outdir,"/enrich_GO_MF.pdf",sep="")
pdf(file=outfile1,width = 10,height = 10)
ridgeplot(gse_goMF) #+ labs(x = "enrichment distribution")
dev.off()

pdf(outdir%&%"GSEA_enrichment_GO_MF.pdf") 
for(p in 1:50){
  print(gse_goMF$Description[p])
  print(gseaplot(gse_goMF, by = "all", title = gse_goMF$Description[p], geneSetID = p))
  #dev.off()
}
dev.off() 

outfile0=paste(outdir,"gse_goMF.rds",sep="")
saveRDS(gse_goMF,outfile0)

gse_KEGG <- gseKEGG(geneList=gene_list_entrezid, 
                    organism = "hsa", 
                    keyType = "kegg", 
                    nPerm = 10000, 
                    minGSSize = 3, 
                    maxGSSize = 800, 
                    pvalueCutoff = 0.05, 
                    verbose = TRUE, 
                    pAdjustMethod = "BH")


df_combine=cbind(names(gene_list),names(gene_list_entrezid))
colnames(df_combine)=c("symbol","ENTREZID")
outfile_=paste(outdir,"GeneNameMap.csv",sep="")
write.csv(df_combine , outfile_, row.names=FALSE, quote = FALSE)

gse_KEGGdf=as.data.frame(gse_KEGG)
outfile=paste(outdir,"KEGG.csv",sep="")
write.csv(gse_KEGGdf , outfile, row.names=FALSE, quote = FALSE)

outfile1 <-paste(outdir,"/gsea_KEGG.pdf",sep="")
pdf(file=outfile1,width = 10,height = 10)
dotplot(gse_KEGG, showCategory=10, split=".sign") + facet_grid(.~.sign)
dev.off()

outfile1 <-paste(outdir,"/enrich_KEGG.pdf",sep="")
pdf(file=outfile1,width = 10,height = 10)
ridgeplot(gse_KEGG) #+ labs(x = "enrichment distribution")
dev.off()

pdf(outdir%&%"GSEA_enrichment_KEGG.pdf") 
for(p in 1:22){
  print(gse_KEGG$Description[p])
  print(gseaplot(gse_KEGG, by = "all", title = gse_KEGG$Description[p], geneSetID =p))
  #dev.off()
}
dev.off() 

outfile0=paste(outdir,"gse_KEGG.rds",sep="")
saveRDS(gse_KEGG,outfile0)

########################################################################
##########################Check the peaks for distribution ############
#######check the peaks

###For KEGG enrichment####
for(p in 1:dim(gse_KEGGdf)[1]){
  print(gse_KEGGdf$Description[p])
  ptw_check=gse_KEGGdf$Description[p]

  gse_core_sel=gse_KEGGdf[gse_KEGGdf$Description==ptw_check,]$core_enrichment
  gse_sel_list=unlist(strsplit(gse_core_sel,"/"))
  df_combine=as.data.frame(df_combine)
  genename=df_combine[df_combine$ENTREZID%in%gse_sel_list,]$symbol
  
  full_gene_df=as.data.frame(gene_list)
  genename_sel=full_gene_df[genename,]
  names(genename_sel)=genename
  
  # Find local maxima
  find_peaks <- function(density_data) {
    y <- density_data$y
    peaks <- which(diff(sign(diff(y))) == -2) + 1
    return(density_data$x[peaks])
  }
  
  density_data <- density(genename_sel)
  
  # Find local maxima
  peaks <- findpeaks(density_data$y, nups = 1, ndowns = 1, npeaks = 2)
  peak_positions <- density_data$x[peaks[,2]]
  print(peak_positions)
  
  # Segment data based on proximity to peaks
  assign_peak <- function(data, peaks) {
    sapply(data, function(x) {
      peak <- which.min(abs(peaks - x))
      return(peak)
    })
  }
  
  peak_assignments <- assign_peak(genename_sel, peak_positions)
  print(peak_assignments)
  
  # Create a data frame
  gene_data <- data.frame(gene = names(genename_sel), polyA_len = genename_sel, peak = peak_assignments)
  
  file_path <- outdir%&%"KEGG/"%&%gse_KEGGdf$Description[p]%&%"_KEGG_PEAK.txt"
  write.table(gene_data, file=file_path, row.names=FALSE, col.names=c("GENE","PolyA_Len","PEAK_Group"), quote=FALSE,sep="\t")
  
  # Split the data into groups based on peaks
  group1 <- subset(gene_data, peak == 1)
  group2 <- subset(gene_data, peak == 2)
  
  # Example enrichment analysis
  enrichment_analysis <- function(group) {
    return(table(group$gene))
  }
  
  enrichment1 <- enrichment_analysis(group1)
  enrichment2 <- enrichment_analysis(group2)
  
  print("Enrichment for Peak 1:")
  print(enrichment1)
  
  print("Enrichment for Peak 2:")
  print(enrichment2)
}


###For GO-BP enrichment####
for(p in 1:dim(gse_BPdf)[1]){
  print(gse_BPdf$Description[p])
  ptw_check=gse_BPdf$Description[p]
  
  gse_core_sel=gse_BPdf[gse_BPdf$Description==ptw_check,]$core_enrichment
  gse_sel_list=unlist(strsplit(gse_core_sel,"/"))
  df_combine=as.data.frame(df_combine)
  genename=df_combine[df_combine$symbol%in%gse_sel_list,]$symbol
  
  full_gene_df=as.data.frame(gene_list)
  genename_sel=full_gene_df[genename,]
  names(genename_sel)=genename
  
  # Find local maxima
  find_peaks <- function(density_data) {
    y <- density_data$y
    peaks <- which(diff(sign(diff(y))) == -2) + 1
    return(density_data$x[peaks])
  }
  
  density_data <- density(genename_sel)
  
  # Find local maxima
  peaks <- findpeaks(density_data$y, nups = 1, ndowns = 1, npeaks = 2)
  peak_positions <- density_data$x[peaks[,2]]
  print(peak_positions)
  
  # Segment data based on proximity to peaks
  assign_peak <- function(data, peaks) {
    sapply(data, function(x) {
      peak <- which.min(abs(peaks - x))
      return(peak)
    })
  }
  
  peak_assignments <- assign_peak(genename_sel, peak_positions)
  print(peak_assignments)
  
  # Create a data frame
  gene_data <- data.frame(gene = names(genename_sel), polyA_len = genename_sel, peak = peak_assignments)
  
  file_path <- outdir%&%"GO_BP/"%&%gse_BPdf$Description[p]%&%"_GOBP_PEAK.txt"
  write.table(gene_data, file=file_path, row.names=FALSE, col.names=c("GENE","PolyA_Len","PEAK_Group"), quote=FALSE,sep="\t")
  
  # Split the data into groups based on peaks
  group1 <- subset(gene_data, peak == 1)
  group2 <- subset(gene_data, peak == 2)
  
  # Example enrichment analysis
  enrichment_analysis <- function(group) {
    return(table(group$gene))
  }
  
  enrichment1 <- enrichment_analysis(group1)
  enrichment2 <- enrichment_analysis(group2)
  
  print("Enrichment for Peak 1:")
  print(enrichment1)
  
  print("Enrichment for Peak 2:")
  print(enrichment2)
}


###For GO-CC enrichment####
for(p in 1:dim(gse_CCdf)[1]){
  print(gse_CCdf$Description[p])
  ptw_check=gse_CCdf$Description[p]
  
  gse_core_sel=gse_CCdf[gse_CCdf$Description==ptw_check,]$core_enrichment
  gse_sel_list=unlist(strsplit(gse_core_sel,"/"))
  df_combine=as.data.frame(df_combine)
  genename=df_combine[df_combine$symbol%in%gse_sel_list,]$symbol
  
  full_gene_df=as.data.frame(gene_list)
  genename_sel=full_gene_df[genename,]
  names(genename_sel)=genename
  
  # Find local maxima
  find_peaks <- function(density_data) {
    y <- density_data$y
    peaks <- which(diff(sign(diff(y))) == -2) + 1
    return(density_data$x[peaks])
  }
  
  density_data <- density(genename_sel)
  
  # Find local maxima
  peaks <- findpeaks(density_data$y, nups = 1, ndowns = 1, npeaks = 2)
  peak_positions <- density_data$x[peaks[,2]]
  print(peak_positions)
  
  # Segment data based on proximity to peaks
  assign_peak <- function(data, peaks) {
    sapply(data, function(x) {
      peak <- which.min(abs(peaks - x))
      return(peak)
    })
  }
  
  peak_assignments <- assign_peak(genename_sel, peak_positions)
  print(peak_assignments)
  
  # Create a data frame
  gene_data <- data.frame(gene = names(genename_sel), polyA_len = genename_sel, peak = peak_assignments)
  
  file_path <- outdir%&%"GO_CC/"%&%gse_CCdf$Description[p]%&%"_GOCC_PEAK.txt"
  write.table(gene_data, file=file_path, row.names=FALSE, col.names=c("GENE","PolyA_Len","PEAK_Group"), quote=FALSE,sep="\t")
  
  # Split the data into groups based on peaks
  group1 <- subset(gene_data, peak == 1)
  group2 <- subset(gene_data, peak == 2)
  
  # Example enrichment analysis
  enrichment_analysis <- function(group) {
    return(table(group$gene))
  }
  
  enrichment1 <- enrichment_analysis(group1)
  enrichment2 <- enrichment_analysis(group2)
  
  print("Enrichment for Peak 1:")
  print(enrichment1)
  
  print("Enrichment for Peak 2:")
  print(enrichment2)
}

###For GO-MF enrichment####
for(p in 1:dim(gse_MFdf)[1]){
  print(gse_MFdf$Description[p])
  ptw_check=gse_MFdf$Description[p]
  
  gse_core_sel=gse_MFdf[gse_MFdf$Description==ptw_check,]$core_enrichment
  gse_sel_list=unlist(strsplit(gse_core_sel,"/"))
  df_combine=as.data.frame(df_combine)
  genename=df_combine[df_combine$symbol%in%gse_sel_list,]$symbol
  
  full_gene_df=as.data.frame(gene_list)
  genename_sel=full_gene_df[genename,]
  names(genename_sel)=genename
  
  # Find local maxima
  find_peaks <- function(density_data) {
    y <- density_data$y
    peaks <- which(diff(sign(diff(y))) == -2) + 1
    return(density_data$x[peaks])
  }
  
  density_data <- density(genename_sel)
  
  # Find local maxima
  peaks <- findpeaks(density_data$y, nups = 1, ndowns = 1, npeaks = 2)
  peak_positions <- density_data$x[peaks[,2]]
  print(peak_positions)
  
  # Segment data based on proximity to peaks
  assign_peak <- function(data, peaks) {
    sapply(data, function(x) {
      peak <- which.min(abs(peaks - x))
      return(peak)
    })
  }
  
  peak_assignments <- assign_peak(genename_sel, peak_positions)
  print(peak_assignments)
  
  # Create a data frame
  gene_data <- data.frame(gene = names(genename_sel), polyA_len = genename_sel, peak = peak_assignments)
  
  file_path <- outdir%&%"GO_MF/"%&%gse_MFdf$Description[p]%&%"_GOMF_PEAK.txt"
  write.table(gene_data, file=file_path, row.names=FALSE, col.names=c("GENE","PolyA_Len","PEAK_Group"), quote=FALSE,sep="\t")
  
  # Split the data into groups based on peaks
  group1 <- subset(gene_data, peak == 1)
  group2 <- subset(gene_data, peak == 2)
  
  # Example enrichment analysis
  enrichment_analysis <- function(group) {
    return(table(group$gene))
  }
  
  enrichment1 <- enrichment_analysis(group1)
  enrichment2 <- enrichment_analysis(group2)
  
  print("Enrichment for Peak 1:")
  print(enrichment1)
  
  print("Enrichment for Peak 2:")
  print(enrichment2)
}

#########################################################
#########################################################
#########################################################
###ADD sort for different subsample
virus_gene_list=gene_df$virus_mean
bac_gene_list=gene_df$bacteria_mean
full_gene_list = gene_df$mean

outdir="/gsea_new/ALL/"
gene_list=full_gene_list
names(gene_list) = gene_df$gene
gene_list = sort(gene_list, decreasing = TRUE) ##From polyA length largest to smallest
gene_list = gene_list[!duplicated(names(gene_list))]
head(gene_list)

DEG.entrez_id = mapIds(x = org.Hs.eg.db,
                       keys = names(gene_list),
                       keytype = "SYMBOL",
                       column = "ENTREZID")
gene_list_entrezid=gene_list
names(gene_list_entrezid)=DEG.entrez_id



##https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/
gse_go <- gseGO(geneList=gene_list, 
                ont ="ALL", 
                keyType = "SYMBOL", 
                nPerm = 10000, 
                minGSSize = 3, 
                maxGSSize = 800, 
                pvalueCutoff = 0.05, 
                verbose = TRUE, 
                OrgDb = org.Hs.eg.db, 
                pAdjustMethod = "BH")


gse_df=as.data.frame(gse_go)
outfile=paste(outdir,"GO_ALL.csv",sep="")
write.csv(gse_df , outfile, row.names=FALSE, quote = FALSE)

outfile1 <-paste(outdir,"/gsea_GO_ALL.pdf",sep="")
pdf(file=outfile1,width = 10,height = 10)
dotplot(gse_go, showCategory=10, split=".sign") + facet_grid(.~.sign)
dev.off()

outfile1 <-paste(outdir,"/enrich_GO_ALL.pdf",sep="")
pdf(file=outfile1,width = 10,height = 10)
ridgeplot(gse_go) + labs(x = "enrichment distribution")
dev.off()

pdf(outdir%&%"GSEA_enrichment_GO_ALL.pdf") 
for(p in 1:50){
  print(gse_go$Description[p])
  print(gseaplot(gse_go, by = "all", title = gse_go$Description[p], geneSetID = p))
  #dev.off()
}
dev.off() 

outfile0=paste(outdir,"gse_go.rds",sep="")
saveRDS(gse_go,outfile0)


gse_goBP <- gseGO(geneList=gene_list, 
                  ont ="BP", 
                  keyType = "SYMBOL", 
                  nPerm = 10000, 
                  minGSSize = 3, 
                  maxGSSize = 800, 
                  pvalueCutoff = 0.05, 
                  verbose = TRUE, 
                  OrgDb = org.Hs.eg.db, 
                  pAdjustMethod = "BH")

gse_BPdf=as.data.frame(gse_goBP)
outfile=paste(outdir,"GO_BP.csv",sep="")
write.csv(gse_BPdf , outfile, row.names=FALSE, quote = FALSE)

outfile1 <-paste(outdir,"/gsea_GO_BP.pdf",sep="")
pdf(file=outfile1,width = 10,height = 10)
dotplot(gse_goBP, showCategory=10, split=".sign") + facet_grid(.~.sign)
dev.off()

outfile1 <-paste(outdir,"/enrich_GO_BP.pdf",sep="")
pdf(file=outfile1,width = 10,height = 10)
ridgeplot(gse_goBP) + labs(x = "enrichment distribution")
dev.off()

pdf(outdir%&%"GSEA_enrichment_GO_BP.pdf") 
for(p in 1:50){
  print(gse_goBP$Description[p])
  print(gseaplot(gse_goBP, by = "all", title = gse_goBP$Description[p], geneSetID = p))
  #dev.off()
}
dev.off() 

outfile0=paste(outdir,"gse_goBP.rds",sep="")
saveRDS(gse_goBP,outfile0)


gse_goCC <- gseGO(geneList=gene_list, 
                  ont ="CC", 
                  keyType = "SYMBOL", 
                  nPerm = 10000, 
                  minGSSize = 3, 
                  maxGSSize = 800, 
                  pvalueCutoff = 0.05, 
                  verbose = TRUE, 
                  OrgDb = org.Hs.eg.db, 
                  pAdjustMethod = "BH")

gse_CCdf=as.data.frame(gse_goCC)
outfile=paste(outdir,"GO_CC.csv",sep="")
write.csv(gse_CCdf , outfile, row.names=FALSE, quote = FALSE)

outfile1 <-paste(outdir,"/gsea_GO_CC.pdf",sep="")
pdf(file=outfile1,width = 10,height = 10)
dotplot(gse_goCC, showCategory=10, split=".sign") + facet_grid(.~.sign)
dev.off()

outfile1 <-paste(outdir,"/enrich_GO_CC.pdf",sep="")
pdf(file=outfile1,width = 10,height = 10)
ridgeplot(gse_goCC) + labs(x = "enrichment distribution")
dev.off()

pdf(outdir%&%"GSEA_enrichment_GO_CC.pdf") 
for(p in 1:50){
  print(gse_goCC$Description[p])
  print(gseaplot(gse_goCC, by = "all", title = gse_goCC$Description[p], geneSetID = p))
  #dev.off()
}
dev.off() 

outfile0=paste(outdir,"gse_goCC.rds",sep="")
saveRDS(gse_goCC,outfile0)


gse_goMF <- gseGO(geneList=gene_list, 
                  ont ="MF", 
                  keyType = "SYMBOL", 
                  nPerm = 10000, 
                  minGSSize = 3, 
                  maxGSSize = 800, 
                  pvalueCutoff = 0.05, 
                  verbose = TRUE, 
                  OrgDb = org.Hs.eg.db, 
                  pAdjustMethod = "BH")

gse_MFdf=as.data.frame(gse_goMF)
outfile=paste(outdir,"GO_MF.csv",sep="")
write.csv(gse_MFdf , outfile, row.names=FALSE, quote = FALSE)

outfile1 <-paste(outdir,"/gsea_GO_MF.pdf",sep="")
pdf(file=outfile1,width = 10,height = 10)
dotplot(gse_goMF, showCategory=10, split=".sign") + facet_grid(.~.sign)
dev.off()

outfile1 <-paste(outdir,"/enrich_GO_MF.pdf",sep="")
pdf(file=outfile1,width = 10,height = 10)
ridgeplot(gse_goMF) + labs(x = "enrichment distribution")
dev.off()

pdf(outdir%&%"GSEA_enrichment_GO_MF.pdf") 
for(p in 1:50){
  print(gse_goMF$Description[p])
  print(gseaplot(gse_goMF, by = "all", title = gse_goMF$Description[p], geneSetID = p))
  #dev.off()
}
dev.off() 

outfile0=paste(outdir,"gse_goMF.rds",sep="")
saveRDS(gse_goMF,outfile0)

gse_KEGG <- gseKEGG(geneList=gene_list_entrezid, 
                    organism = "hsa", 
                    keyType = "kegg", 
                    nPerm = 10000, 
                    minGSSize = 3, 
                    maxGSSize = 800, 
                    pvalueCutoff = 0.05, 
                    verbose = TRUE, 
                    pAdjustMethod = "BH")

df_combine=cbind(names(gene_list),names(gene_list_entrezid))
colnames(df_combine)=c("symbol","ENTREZID")
outfile_=paste(outdir,"GeneNameMap.csv",sep="")
write.csv(df_combine , outfile_, row.names=FALSE, quote = FALSE)

gse_KEGGdf=as.data.frame(gse_KEGG)
outfile=paste(outdir,"KEGG.csv",sep="")
write.csv(gse_KEGGdf , outfile, row.names=FALSE, quote = FALSE)

outfile1 <-paste(outdir,"/gsea_KEGG.pdf",sep="")
pdf(file=outfile1,width = 10,height = 10)
dotplot(gse_KEGG, showCategory=10, split=".sign") + facet_grid(.~.sign)
dev.off()

outfile1 <-paste(outdir,"/enrich_KEGG.pdf",sep="")
pdf(file=outfile1,width = 10,height = 10)
ridgeplot(gse_KEGG) #+ labs(x = "enrichment distribution")
dev.off()

pdf(outdir%&%"GSEA_enrichment_KEGG.pdf") 
for(p in 1:23){
  print(gse_KEGG$Description[p])
  print(gseaplot(gse_KEGG, by = "all", title = gse_KEGG$Description[p], geneSetID =p))
  #dev.off()
}
dev.off() 

outfile0=paste(outdir,"gse_KEGG.rds",sep="")
saveRDS(gse_KEGG,outfile0)
