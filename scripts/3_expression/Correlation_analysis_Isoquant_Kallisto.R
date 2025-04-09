.libPaths(libdir)

"%&%" <- function(a,b) paste(a,b, sep = "")


#ggforce for ggsina
install=F
if(install){
  install.packages("BiocManager")
  BiocManager::install("VGAM", version = "3.8")
  BiocManager::install("ggplot2")
  BiocManager::install("gplots")
  BiocManager::install("jsonlite")
  BiocManager::install("gridExtra")
  BiocManager::install("GGally")
  BiocManager::install("binom")
  BiocManager::install("biomaRt", version = "2.4")
}
library(VGAM)
library(binom)
library(ggplot2)
#library(heatmap2)
library(gplots)
#lbirary(beeswarm)
library(jsonlite)
library(gridExtra)
#library(biomaRt)
#library(readxl)
library(writexl)
#library(GGally)
#install.packages(paste(rhome1,"ggally", sep="/"),  repos = NULL, type="source")
library(GGally)
library(ggpmisc)
library(reshape)

#rm(list = ls())
rhome1="/data"
rapidshome = paste(rhome1,"rapids_rna_seq", sep="/");
#rhome1 = "C:/Users/ljmco/bitbucket" 
source("Correlation_functions.R")

datahome= rhome1%&%"humanRNA/"

#resdir = getResultDir( "results")
resdir="results/Isoquant_Illumina_tpm"
dir.create(resdir)
betag = c("HBB" , "HBG1" , "HBG2" , "HBA1" , "HBA2" , "HBD", "HBE1","AC104389.3","AC104389.4","AC104389.5","AC104389.6")


####NOTES: Finally, I understand. We are going to exclude nanopore direct cDNA sequencing which is under HTSeq_count (cDNA); 
## We are comparing nanopore RNA sequencing with Illumina cDNA sequencing (under Illumina_Kallisto/)

#data_dir = "HtSeq-Count_files";
#data_dirK = "Kallisto";
data_dir = datahome%&%"isoquant/gene_tpm/";
data_dirK = datahome%&%"kallisto_illumina/"; ##Don't really understand the name under Illumina_Kallisto2/ Under Illumina_Kallisto only contain cDNA

exons = read.table(datahome%&%"human_exon_length_v35.txt", head=T);
mt = read.table(datahome%&%"reference_genome_set_chr_mitochondria.txt", sep="\t", head=T, as.is = T)
mt_name = mt$Approved.Symbol

gc = read.table(datahome%&%"gencode.v35.transcripts_gc.txt", head=T, comment.char='#', as.is=T)

exon_count=read.csv(datahome%&%"exon_gene_count.csv", head=T)

types = list("protein_coding", "lincRNA", c("miRNA", "snoRNA", "snRNA"))

##READ Illumina and Kallisto data from HTSeq directory.  
#This will normaliseIllumina data
file1=read.table(data_dir%&%samples_HT[1]%&%"pass_count_gene_name.txt",header=F)
exons_HT=exons[exons$Geneid%in%file1$V1,]

all_data_norm_HT = readAllData_new(data_dir, samples_HT, exons_HT) 
#all_data_unnorm_HT = readAllData(data_dir, samples, ".", exons, FALSE) ## this is not normalised for length

##Reading Kallisto results for Illumina only
##In the readAllDataKallisto_new, I set zero_thresh = length(samples) and missing_thresh = length(samples), would this be allright?
all_data_norm_K = readAllDataKallisto_new(data_dirK, samples_K, types = NULL, cnt_type="tpm",incl = TRUE , illumina_only = F, gc= gc)
.dim1<-function(f) dim(f$exons)
.dim2<-function(f) dim(f$gc)

#remove null exons(due to just one entry in category)
nullexons = which(unlist(lapply(lapply(all_data_norm_K, .dim1),is.null)))
if(length(nullexons)>0){
  all_data_norm_K = all_data_norm_K[-nullexons]
}
lapply(all_data_norm_K, .dim2)

all_data_norm_K_ = all_data_norm_K
mi = unlist(lapply(all_data_norm_K_, function(f) dim(f$exons)[2]))
all_data_norm_K_ = all_data_norm_K_[mi==max(mi)]


for(i in 1:length(all_data_norm_K)){
  all_data_norm_K_[[i]]$exons =  cbind(all_data_norm_K[[i]]$exons,all_data_norm_K[[i]]$gc$X.GC)
}
all_data_norm_K__ = all_data_norm_K_

params = fromJSON(paste(rapidshome,"biotypes.json", sep="/"), simplifyDataFrame  = F)
for(i in 1:length(params)){
  print(params[[i]])
  all_data_norm_K__  = mergeTypes1(all_data_norm_K__, params[[i]], names(params)[[i]])
}

all_data_norm_K2 = lapply(all_data_norm_K__, mergeSameName)

##JN:Q:till this point under HTSeq_count/ only contain nanopore RNA, cDNA results and under Illumina_Kallisto/ only contain illumina cDNA result
## Where is the illumina RNA results??
## However, after merging the hybrid_ does contain illumina RNA? Does illumina RNA==illumina cDNA, and what is called cDNA means nanopore cDNA??
## If so, in this case, we actually want to compare illumina cDNA with nanopore RNA????
hybrid = list()
hybrid_ = list()
hybrid3 = list()
nmes = names(all_data_norm_K2)
inds= 1:length(all_data_norm_K2)
for(i in 1:length(all_data_norm_K2))  hybrid_[[i]]  =  mergeHTKallisto_new(all_data_norm_K2[[i]], all_data_norm_HT, replace=T,datatype=names(all_data_norm_K2)[[i]]) ####################

for(i in 1:length(all_data_norm_K2))  hybrid[[i]]  =  mergeHTKallisto_new(all_data_norm_K2[[i]], all_data_norm_HT, replace=T,datatype=names(all_data_norm_K2)[[i]])

for(i in 1:length(all_data_norm_K_)) {
  print(i)
  hybrid3[[i]]  =  mergeHTKallisto_new(all_data_norm_K_[[i]], all_data_norm_HT, replace=T,datatype=names(all_data_norm_K_)[[i]])
  #if(is.null(hybrid3[[i]]$data1)) hybrid3[[i]] = NULL;
  #dim(hybrid[[k]]$data1)[1]>0)
}
names(hybrid3) = names(all_data_norm_K_)
hybrid3 = hybrid3[!unlist(lapply(hybrid3, is.null))]
names(hybrid) = nmes
names(hybrid_) = nmes
#data_exon_norm = normaliseByExonLength(all_data_norm_HT$data1,exons)

#inds = which(dimnames(result)[[1]] %in% mt_name)
#		result1 = result[inds,,drop=F]
hybrid = addExonCount(hybrid, exon_count) ##add "exonL" to exon under each datatyope
hybrid3 = addExonCount(hybrid3, exon_count)
#hybrid1 = mergeIntoOthers(hybrid, len_thresh = 10, sumthresh = 1000, sumthresh_sum = 5000)
#hybrid_1 = mergeIntoOthers(hybrid_, len_thresh = 10, sumthresh = 1000, sumthresh_sum = 5000)
#hybrid_ = hybrid_1
#hybrid = hybrid1
lens = lapply(hybrid,.getlen)
allsums = lapply(hybrid,.getsums)

#hybrid = normaliseH(hybrid, sum_all_types=T) 
#hybrid_ = normaliseH(hybrid)

#hybrid3 = normaliseH(hybrid3)
hybrid =removeGlobin(hybrid, annotate=F) ##what does the removeGlobin do?
hybrid3 =removeGlobin(hybrid3, annotate=F)

hybrid=hybrid[-length(hybrid)]
hybrid3=hybrid_[-length(hybrid3)]

h = mergeAll(hybrid) #, "S"))
##h = normaliseH_(h)

hybrid[[length(hybrid)+1]] = h
names(hybrid)[[length(hybrid)]] = "combined1" ##combined1 contains all the genes from all the datatyes 

#############
nmes = names(hybrid)
res_all_hybrid = list()
res_all_hybrid_ = list()
res_all_hybrid3 = list()


hybrid=hybrid[-3]
saveRDS(hybrid,resdir%&%"/hybrid.rds")

starti = 1
for(i in starti:length(hybrid)) {
  print(names(hybrid)[i])
  minp=1e-15   ## ??CHECK
  res_all_hybrid[[i]] = run_all_new(hybrid[[i]], data_norm=F,samples_HT, paste(resdir,"/", names(hybrid)[i], sep=""),topn = NA, sig = 1e-5, gt = F, print = T, doheatmap=T) 
  ##############
  outf2 = paste(resdir, names(hybrid)[i], "hybrid.txt", sep="/")
  write.table(ftable(hybrid3[[i]]$data_prop), file=outf2, append=F, quote=FALSE, sep="\t", row.names=F, col.names=T)
  
  outf_pearson = paste(resdir, names(hybrid)[i], "corr_pearson.txt", sep="/")
  outf_spearman = paste(resdir, names(hybrid)[i], "corr_spearman.txt", sep="/")
  corr_pearson= ftable(res_all_hybrid[[i]]$corr)
  corr_spearman= ftable(res_all_hybrid[[i]]$corr1)
  write.table(corr_pearson,file= outf_pearson, quote = FALSE, sep="\t", row.names =F, col.names =T)
  write.table(corr_spearman,file= outf_spearman, quote = FALSE, sep="\t", row.names =F, col.names =T)
}
names(res_all_hybrid) =names(hybrid)
res_all_hybrid=saveRDS(res_all_hybrid,resdir%&%"/res_all_hybrid.rds")

########################OLD code From Lachlan##################################

res_all_hybrid=readRDS(resdir%&%"/res_all_hybrid.rds")

DEHeatmaps(res_all_hybrid, resdir, nme="mRNA", thresh = 1e-5, max = 100, resname = "DE_heatmap.pdf")


enrich = calcEnrichAll2(res_all_hybrid,paste(resdir, "mt_enrichment", sep="/"), mt_name, pv_thresh = 1e-5)


useqnorm = F
##comps=names(res_all_hybrid[[index1]]) ---> [1] "CDNAvsILL" "RNAvsILL"  "CDNAvsRNA"
##only "RNAvsILL" run without error
comps1=c("RNAvsILL")
pvsl=plot_pvs(res_all_hybrid,resdir, nme="mRNA", useqnorm=useqnorm, print = T,comps=comps1,text_size=2)
##This will generate pv_corr.pdf under ./mRNA/

##what does DES means??
##comps = (names(res_all_hybrid[[index1]]$DES)) ---> comps:[1] "6201vs6679"   "GE031vsGE044" "6201vsGE031"
##Only "6201vs6679" run without error
comps_des=c("6201_vs6679_")
plot_pvs_DES(res_all_hybrid, resdir, nme="mRNA", nudge = 1e-30, useqnorm=useqnorm,comps=comps_des,text_size=5)
##This will generate pv_corr1.pdf under ./mRNA/


plotAllCorsCombined2(h, paste(resdir, "combined_",sep="/"), method="pearson" , removeZero=F,samples=samples_HT)

##No below function
##leplotAllCors(hybrid, paste(resdir, "combined",sep="/"))

lev = levels(as.factor(h$type))
lev = c("mRNA", lev[-which(lev=='mRNA')])
torm  = list( lev[-which(lev=="mRNA")], 'mRNA')
pdf(paste(resdir,"density.pdf", sep="/"), paper="a4r", width = 20,height = 20)
plotDensity(h,list(1:2, 3:4), rev(lev), rm = torm)
plotDensity(h,list(5:6, 7:8), rev(lev), rm = torm)
plotDensity(h,list(9:10, 11:12), rev(lev), rm = torm)
plotDensity(h,list(13:14, 15:16), rev(lev), rm = torm)
plotDensity(h,list(17:18, 19:20), rev(lev), rm = torm)
plotDensity(h,list(21:22, 23:24), rev(lev), rm = torm)
dev.off()

if(FALSE){
  h_ = mergeAll(hybrid_) #, "S"))
  plotAllCorsCombined(h_, paste(resdir, "combined_",sep="/"))
  plotAllCors(hybrid_, paste(resdir, "combined_",sep="/"))
  
  
  h3 = mergeAll(hybrid3) #, "S"))
  plotAllCorsCombined(h3, paste(resdir, "combined_",sep="/"))
  plotAllCors(hybrid3, paste(resdir, "combined_",sep="/"))
}


pdf(paste(resdir, "regplots.pdf",sep="/")) 
thresh1 = c(1e9,1e9)
levs = levels(hybrid$mRNA$exons$type)
#levs = "protein_coding"
levs = grep("globin", levs, v=T, inv=T)
#levs = "protein_coding"

################correlation on gene to gene level ####################

gene_df=read.table("/data/HSIAO_HOUSEKEEPING_GENES.v2023.2.Hs.grp",head=T)
HK_genelist=gene_df$HSIAO_HOUSEKEEPING_GENES

hybrid=readRDS(resdir%&%"/hybrid.rds")

hybrid_T=list()
hybrid_T_sel=list()
res_hybrid_T=list()
corr_p_list=list()
corr_s_list=list()

hybrid_T_sel2=list()
res_hybrid_T2=list()
corr_p_list2=list()
corr_s_list2=list()

i=5

hybrid_T=lapply( hybrid[[i]], function(obj) {reverse_matrix(obj)})

hybrid_T_sel=hybrid_T


#gene_list <- intersect(colnames(hybrid_T$data1),HK_genelist)
gene_list <- colnames(hybrid_T$data1)
print(length(gene_list))

hybrid_T_sel$data1 <- hybrid_T$data1[,gene_list]
hybrid_T_sel$exons<-hybrid_T$exons[,gene_list]
hybrid_T_sel$data_prop <- hybrid_T$data_prop[,gene_list]
hybrid_T_sel$data_low <- hybrid_T$data_low[,gene_list]
hybrid_T_sel$data_high <- hybrid_T$data_high[,gene_list]

res_hybrid_T = run_corr_T_new(hybrid_T_sel, data_norm=F) 
print("finish corr") 


corr_p_list=res_hybrid_T$corr_p
names(corr_p_list)=gene_list
corr_p_list=sort(corr_p_list, decreasing = TRUE, na.last = T)

corr_s_list=res_hybrid_T$corr_s
names(corr_s_list)=gene_list
corr_s_list=sort(corr_s_list, decreasing = TRUE, na.last = T)


saveRDS(hybrid_T,resdir%&%"/mRNA/hybrid_T.rds")
saveRDS(res_hybrid_T,resdir%&%"/mRNA/res_hybrid_T.rds")
saveRDS(corr_p_list,resdir%&%"/mRNA/corr_p_list.rds")
saveRDS(corr_s_list,resdir%&%"/mRNA/corr_s_list.rds")

df_pcorr=data.frame(corr_p_list)
df_scorr=data.frame(corr_s_list)

df_corr=cbind(df_pcorr,df_scorr)
colnames(df_corr)=c("pearson","spearman")
write.table(df_corr, file = resdir%&%"/mRNA/corr_genes.txt", row.names =T, col.names = T, quote = FALSE,sep="\t")


hybrid_T_sel2=list()
res_hybrid_T2=list()
corr_p_list2=list()
corr_s_list2=list()

hybrid_T=lapply( hybrid[[i]], function(obj) {reverse_matrix(obj)})

hybrid_T_sel2=hybrid_T

gene_list2 <- intersect(colnames(hybrid_T$data1),HK_genelist)
#gene_list <- colnames(hybrid_T$data1)
print(length(gene_list2))

hybrid_T_sel2$data1 <- hybrid_T$data1[,gene_list2]
hybrid_T_sel2$exons<-hybrid_T$exons[,gene_list2]
hybrid_T_sel2$data_prop <- hybrid_T$data_prop[,gene_list2]
hybrid_T_sel2$data_low <- hybrid_T$data_low[,gene_list2]
hybrid_T_sel2$data_high <- hybrid_T$data_high[,gene_list2]

res_hybrid_T2 = run_corr_T_new(hybrid_T_sel2, data_norm=F) 
print("finish corr") 


corr_p_list2=res_hybrid_T2$corr_p
names(corr_p_list2)=gene_list2
corr_p_list2=sort(corr_p_list2, decreasing = TRUE, na.last = T)

corr_s_list2=res_hybrid_T2$corr_s
names(corr_s_list2)=gene_list2
corr_s_list2=sort(corr_s_list2, decreasing = TRUE, na.last = T)

saveRDS(res_hybrid_T2,resdir%&%"/mRNA/res_hybrid_T_hg.rds")
saveRDS(corr_p_list2,resdir%&%"/mRNA/corr_p_list_hg.rds")
saveRDS(corr_s_list2,resdir%&%"/mRNA/corr_s_list_hg.rds")

df_pcorr2=data.frame(corr_p_list2)
df_scorr2=data.frame(corr_s_list2)

df_corr2=cbind(df_pcorr2,df_scorr2)
colnames(df_corr2)=c("pearson","spearman")
write.table(df_corr2, file = resdir%&%"/mRNA/corr_genes_hg.txt", row.names =T, col.names = T, quote = FALSE,sep="\t")


##########Only consider genes with average read count >10 for gene to gene correlations#############
hybrid=readRDS(resdir%&%"/hybrid.rds")

filter=100

data_all=hybrid$mRNA$data1
name_list=unlist(lapply(dimnames(data_all)[[2]], function(str)strsplit(str,"__")[[1]][2]))
name_list=dimnames(data_all)[[2]]
ill_list=name_list[grep("illumina",name_list)]
rna_list=name_list[grep("illumina",name_list,inv=T)]
data_ill=t(data_all[,ill_list])
data_rna=t(data_all[,rna_list])
ill_gene=names(which(apply(data_ill,2,mean)>filter))
rna_gene=names(which(apply(data_rna,2,mean)>filter))
sel_gene_list = intersect(ill_gene, rna_gene)
#sel_gene_list=intersect(HK_genelist,sel_gene_list)


i=5

hybrid_T=list()
hybrid_T_sel2=list()
res_hybrid_T2=list()
corr_p_list2=list()
corr_s_list2=list()

hybrid_T=lapply( hybrid[[i]], function(obj) {reverse_matrix(obj)})

hybrid_T_sel2=hybrid_T

gene_list2 <- intersect(colnames(hybrid_T$data1),sel_gene_list)
#gene_list <- colnames(hybrid_T$data1)
print(length(gene_list2))

hybrid_T_sel2$data1 <- hybrid_T$data1[,gene_list2]
hybrid_T_sel2$exons<-hybrid_T$exons[,gene_list2]
hybrid_T_sel2$data_prop <- hybrid_T$data_prop[,gene_list2]
hybrid_T_sel2$data_low <- hybrid_T$data_low[,gene_list2]
hybrid_T_sel2$data_high <- hybrid_T$data_high[,gene_list2]

res_hybrid_T2 = run_corr_T_new(hybrid_T_sel2, data_norm=F) 
print("finish corr") 


corr_p_list2=res_hybrid_T2$corr_p
names(corr_p_list2)=gene_list2
corr_p_list2=sort(corr_p_list2, decreasing = TRUE, na.last = T)

corr_s_list2=res_hybrid_T2$corr_s
names(corr_s_list2)=gene_list2
corr_s_list2=sort(corr_s_list2, decreasing = TRUE, na.last = T)

saveRDS(res_hybrid_T2,resdir%&%"/mRNA/res_hybrid_T_"%&%filter%&%".rds")
saveRDS(corr_p_list2,resdir%&%"/mRNA/corr_p_list_"%&%filter%&%".rds")
saveRDS(corr_s_list2,resdir%&%"/mRNA/corr_s_list_"%&%filter%&%".rds")

df_pcorr2=data.frame(corr_p_list2)
df_scorr2=data.frame(corr_s_list2)

df_corr2=cbind(df_pcorr2,df_scorr2)
colnames(df_corr2)=c("pearson","spearman")
write.table(df_corr2, file = resdir%&%"/mRNA/corr_genes_"%&%filter%&%".txt", row.names =T, col.names = T, quote = FALSE,sep="\t")

##########Only consider genes with all samples' read count >10 for gene to gene correlations#############
#hybrid=readRDS(resdir%&%"/hybrid.rds")

#filter=30

data_all=hybrid$mRNA$data1
name_list=unlist(lapply(dimnames(data_all)[[2]], function(str)strsplit(str,"__")[[1]][2]))
name_list=dimnames(data_all)[[2]]
ill_list=name_list[grep("illumina",name_list)]
rna_list=name_list[grep("illumina",name_list,inv=T)]
data_ill=t(data_all[,ill_list])
data_rna=t(data_all[,rna_list])

cols2keep_ill <- colSums(data_ill > filter) == nrow(data_ill)
filtered_ill <- subset(data_ill, select = cols2keep_ill)
cols2keep_rna <- colSums(data_rna > filter) == nrow(data_rna)
filtered_rna <- subset(data_rna, select = cols2keep_rna)
ill_gene=colnames(filtered_ill)
rna_gene=colnames(filtered_rna)
sel_gene_list = intersect(ill_gene, rna_gene)
#sel_gene_list=intersect(HK_genelist,sel_gene_list)

i=5

hybrid_T=list()
hybrid_T_sel2=list()
res_hybrid_T2=list()
corr_p_list2=list()
corr_s_list2=list()

hybrid_T=lapply( hybrid[[i]], function(obj) {reverse_matrix(obj)})

hybrid_T_sel2=hybrid_T

gene_list2 <- intersect(colnames(hybrid_T$data1),sel_gene_list)
#gene_list <- colnames(hybrid_T$data1)
print(length(gene_list2))

hybrid_T_sel2$data1 <- hybrid_T$data1[,gene_list2]
hybrid_T_sel2$exons<-hybrid_T$exons[,gene_list2]
hybrid_T_sel2$data_prop <- hybrid_T$data_prop[,gene_list2]
hybrid_T_sel2$data_low <- hybrid_T$data_low[,gene_list2]
hybrid_T_sel2$data_high <- hybrid_T$data_high[,gene_list2]

res_hybrid_T2 = run_corr_T_new(hybrid_T_sel2, data_norm=F) 
print("finish corr") 


corr_p_list2=res_hybrid_T2$corr_p
names(corr_p_list2)=gene_list2
corr_p_list2=sort(corr_p_list2, decreasing = TRUE, na.last = T)

corr_s_list2=res_hybrid_T2$corr_s
names(corr_s_list2)=gene_list2
corr_s_list2=sort(corr_s_list2, decreasing = TRUE, na.last = T)

df_pcorr2=data.frame(corr_p_list2)
df_scorr2=data.frame(corr_s_list2)

df_corr2=cbind(df_pcorr2,df_scorr2)
colnames(df_corr2)=c("pearson","spearman")
write.table(df_corr2, file = resdir%&%"/mRNA/corr_genes_idv_"%&%filter%&%".txt", row.names =T, col.names = T, quote = FALSE,sep="\t")



################# sel_gene correlations #####################

#hybrid=readRDS(resdir%&%"/hybrid.rds")

df_pcorr_full=read.table(paste(resdir,"/mRNA/corr_genes.txt",sep=""))
df_pcorr_na=na.omit(df_pcorr_full)

corr_cutoff_list=c(0.25,0.5,0.75,0.9)

i=5
minp=1e-15
for (pcorr_cutoff in corr_cutoff_list){
  
  pcorr_cutoff=0.25
  
  resdir_sel=paste(resdir,"/", names(hybrid)[i], "/sel_gene/pcorr",as.character(pcorr_cutoff),sep="")
  dir.create(resdir_sel)
  
  sel_gene=rownames(df_pcorr_na[df_pcorr_na$pearson>pcorr_cutoff,])
  
  file_sum <- file(resdir_sel%&%"/summary_gene.txt", "w")
  writeLines("Number of total genes: "%&%as.character(dim(df_pcorr_full)[1])%&%"; Number of non-NA: "%&%as.character(dim(df_pcorr_na)[1])%&%"; Ratio of non-NA: "%&%as.character(dim(df_pcorr_na)[1]/dim(df_pcorr_full)[1]), file_sum)
  writeLines("Number of gene with pcorr > "%&%as.character(pcorr_cutoff)%&%" is "%&%as.character(length(sel_gene))%&%"; Ratio of gene with high pcorr: "%&%as.character(length(sel_gene)/(dim(df_pcorr_na)[1])), file_sum)
  close(file_sum)
  print("Number of total genes: "%&%as.character(dim(df_pcorr_full)[1])%&%"; Number of non-NA: "%&%as.character(dim(df_pcorr_na)[1])%&%"; Ratio of non-NA: "%&%as.character(dim(df_pcorr_na)[1]/dim(df_pcorr_full)[1]))
  print("Number of gene with pcorr > "%&%as.character(pcorr_cutoff)%&%" is "%&%as.character(length(sel_gene))%&%"; Ratio of gene with high pcorr: "%&%as.character(length(sel_gene)/(dim(df_pcorr_na)[1])))
  
  full_gene_list=rownames(df_pcorr_full)
  gene_excl_list=setdiff(full_gene_list, sel_gene)
  
  print(names(hybrid)[i])
  
  res_sel_hybrid=list()
  res_sel_hybrid= run_all_new2(hybrid[[i]], data_norm=F,samples_HT, resdir_sel,topn = NA, sig = 1e-5, excl=gene_excl_list, gt = F, print = T, doheatmap=T) 
  
  outf_pearson = paste(resdir_sel, "/corr_pearson.txt", sep="/")
  outf_spearman = paste(resdir_sel, "/corr_spearman.txt", sep="/")
  corr_pearson= ftable(res_sel_hybrid$corr)
  corr_spearman= ftable(res_sel_hybrid$corr1)
  write.table(corr_pearson,file= outf_pearson, quote = FALSE, sep="\t", row.names =F, col.names =T)
  write.table(corr_spearman,file= outf_spearman, quote = FALSE, sep="\t", row.names =F, col.names =T)
  
  
  corrp_df=read.table(outf_pearson,head=T)
  rownames(corrp_df)=corrp_df[,1]
  corrp_m=as.matrix(corrp_df[,-1])
  
  type="Pearson"
  corr_m=corrp_m
  
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
  
  #gpl
  
  #gpl=corr_heatmap(corrp_m,resdir,type="Pearson")
  
  pdf(resdir_sel%&%"/Peason_Corr_heatmap.pdf")
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
  
  
  #gpl=corr_heatmap(corrsp_m,resdir,type="Spearman")
  
  pdf(resdir_sel%&%"/Spearman_Corr_heatmap.pdf")
  gpl
  dev.off()
  
  
  
  
  resdir_sel=paste(resdir,"/", names(hybrid)[i], "/sel_gene/spcorr",as.character(pcorr_cutoff),sep="")
  dir.create(resdir_sel)
  
  sel_gene=rownames(df_pcorr_na[df_pcorr_na$spearman>pcorr_cutoff,])
  
  file_sum <- file(resdir_sel%&%"/summary_gene.txt", "w")
  writeLines("Number of total genes: "%&%as.character(dim(df_pcorr_full)[1])%&%"; Number of non-NA: "%&%as.character(dim(df_pcorr_na)[1])%&%"; Ratio of non-NA: "%&%as.character(dim(df_pcorr_na)[1]/dim(df_pcorr_full)[1]), file_sum)
  writeLines("Number of gene with spcorr > "%&%as.character(pcorr_cutoff)%&%" is "%&%as.character(length(sel_gene))%&%"; Ratio of gene with high spcorr: "%&%as.character(length(sel_gene)/(dim(df_pcorr_na)[1])), file_sum)
  close(file_sum)
  print("Number of total genes: "%&%as.character(dim(df_pcorr_full)[1])%&%"; Number of non-NA: "%&%as.character(dim(df_pcorr_na)[1])%&%"; Ratio of non-NA: "%&%as.character(dim(df_pcorr_na)[1]/dim(df_pcorr_full)[1]))
  print("Number of gene with spcorr > "%&%as.character(pcorr_cutoff)%&%" is "%&%as.character(length(sel_gene))%&%"; Ratio of gene with high spcorr: "%&%as.character(length(sel_gene)/(dim(df_pcorr_na)[1])))
  
  full_gene_list=rownames(df_pcorr_full)
  gene_excl_list=setdiff(full_gene_list, sel_gene)
  
  print(names(hybrid)[i])
  
  res_sel_hybrid=list()
  res_sel_hybrid= run_all_new2(hybrid[[i]], data_norm=F,samples_HT, resdir_sel,topn = NA, sig = 1e-5, excl=gene_excl_list, gt = F, print = T, doheatmap=T) 
  
  outf_pearson = paste(resdir_sel, "/corr_pearson.txt", sep="/")
  outf_spearman = paste(resdir_sel, "/corr_spearman.txt", sep="/")
  corr_pearson= ftable(res_sel_hybrid$corr)
  corr_spearman= ftable(res_sel_hybrid$corr1)
  write.table(corr_pearson,file= outf_pearson, quote = FALSE, sep="\t", row.names =F, col.names =T)
  write.table(corr_spearman,file= outf_spearman, quote = FALSE, sep="\t", row.names =F, col.names =T)
  
  
  corrp_df=read.table(outf_pearson,head=T)
  rownames(corrp_df)=corrp_df[,1]
  corrp_m=as.matrix(corrp_df[,-1])
  
  type="Pearson"
  corr_m=corrp_m
  
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
  
  #gpl
  
  #gpl=corr_heatmap(corrp_m,resdir,type="Pearson")
  
  pdf(resdir_sel%&%"/Peason_Corr_heatmap.pdf")
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
  
  
  #gpl=corr_heatmap(corrsp_m,resdir,type="Spearman")
  
  pdf(resdir_sel%&%"/Spearman_Corr_heatmap.pdf")
  gpl
  dev.off()
  
}


