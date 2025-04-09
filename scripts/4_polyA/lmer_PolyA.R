.libPaths(libdir)

"%&%" <- function(a,b) paste(a,b, sep = "")

library(VGAM)
library(binom)
library(ggplot2)
library(gplots)
library(jsonlite)
library(gridExtra)
library(GGally)
library(ggpmisc)
library(Rsamtools)
library(tidyverse)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(stringr)
library(stringi)


source("/Correlation_functions.R")

# Path to your BAM file
bam_dir<-"/data/trans_alignment/"
out_dir<-"/data/PolyA/"

rhome1="/data/"
rapidshome = paste(rhome1,"rapids_rna_seq", sep="/");
datahome= rhome1%&%"humanRNA/From_Josh/"
gc = read.table(datahome%&%"gencode.v35.transcripts_gc.txt", head=T, comment.char='#', as.is=T)
exon_count=read.csv(datahome%&%"exon_gene_count.csv", head=T)
data_dirK = datahome%&%"kallisto_illumina/";

norm_poly=F
if (norm_poly){
  resdir="results/Nanocount_normPolyA"
}else{resdir="results/Nanocount_PolyA"}
dir.create(resdir)

metadata <- read.table("/data/sample_info.txt", header = TRUE)
rownames(metadata)=metadata$sample2

logpolyA_combined_df=NULL
for (ss in 1:length(samples_list)){
  sample=samples_list[ss]
  print(sample)
  curr_condition=metadata[metadata$sample2==sample,]$infect
  
  bam_file = bam_dir%&%grep("mapped",grep("bai",grep(sample, grep("bam", dir(bam_dir), v = T), v=T), v=T, inv = T), v=T, inv = T)
  poly_file = bam_dir%&%grep(sample, grep("stats", dir(bam_dir), v = T), v=T)
  
  polyA_df<-read.table(poly_file,head=F, comment.char='#', as.is=T)
  colnames(polyA_df)=c("qname","qwidth","polyA_len")
  
  # Open BAM file
  bamFile <- Rsamtools::BamFile(bam_file)
  seqinfor<-seqinfo(bamFile)
  #aln <- scanBam(bamFile)
  
  what=c("rname", "strand", "pos", "qwidth")
  
  # Read BAM file with 'unmapped' reads excluded
  filtered_reads <-Rsamtools::scanBam(bamFile, param=ScanBamParam(what=c("qname","rname", "strand", "pos"),flag=scanBamFlag(isUnmappedQuery=FALSE)))
  
  # Display some information about the reads
  head(filtered_reads)
  
  filtered_reads_df<-as.data.frame(filtered_reads[[1]])
  
  ##Notice: Not all mapped read has polyA length, also polyA_df also contain unmapped reads
  rpolyA_combine_df<-merge(filtered_reads_df,polyA_df,by="qname")
  
  ##From rpolyA_combine_df,one read can map to multiple transcripts, because Nanopore reads are quite long
  ##We first calculate the normalize PolyA length based on the read length
  ##Normalize the PolyA length based on the read length
  ##normalized_polyA_lengths <- polyA_ratio * mean(read_lengths) where polyA_ratio=polyA_lengths / read_lengths
  rpolyA_combine_df$norm_polyA_len=(rpolyA_combine_df$polyA_len/rpolyA_combine_df$qwidth)*mean(rpolyA_combine_df$qwidth)
  ##For each transcript, we combine all reads that have been mapped to this transcript, and take average for the polyA_length, norm_polyA_length and read_length
  #avg_rpolyA_combine_df<- aggregate(cbind(rpolyA_combine_df$polyA_len, rpolyA_combine_df$norm_polyA_len,rpolyA_combine_df$qwidth) ~ rpolyA_combine_df$rname, data = rpolyA_combine_df, FUN = mean)
  #colnames(avg_rpolyA_combine_df)<-c("transcript_name","avg_polyA_len","avg_norm_polyA_len","avg_read_len")
  #outf_poly=out_dir%&%sample%&%"_polyA_length.txt"
  #write.table(avg_rpolyA_combine_df,file=outf_poly, quote = FALSE, sep="\t", row.names =F, col.names =T)
  
  ##Instead of for each transcript, we combine all reads that have been mapped to the same genes and take average for the polyA_length, norm_polyA_length and read_length
  gene_name_df=data.frame(t(data.frame(strsplit( as.character(rpolyA_combine_df$rname),"\\|"))))
  colnames(gene_name_df)=c("transcript_ID","gene_ID","OTTUMG","OTTUMT","gene_I","gene_name","length","data_type")
  rpolyA_combine_df=cbind(rpolyA_combine_df,gene_name_df$gene_name,gene_name_df$data_type)
  colnames(rpolyA_combine_df)<-c("qname","rname","strand","pos","qwidth","polyA_len","norm_polyA_len","gene_name","data_type")
  ##Here we need to first combine different data_type to subcategory
  params = fromJSON(paste(rapidshome,"biotypes.json", sep="/"), simplifyDataFrame  = F)
  
  reversed_param_map <- list()
  # Iterate through the original list and reverse the key-value pairs
  for (category in names(params)) {
    for (rna_type in params[[category]]) {
      reversed_param_map[[rna_type]] <- category
    }
  }
  
  print(reversed_param_map)
  
  rpolyA_combine_df$data_param=reversed_param_map[rpolyA_combine_df$data_type]
  colnames(rpolyA_combine_df)<-c("qname","rname","strand","pos","qwidth","polyA_len","norm_polyA_len","gene_name","data_type","data_params")
  
  ##For each gene, we combine all reads that have been mapped to this transcript, and take average for the polyA_length, norm_polyA_length and read_length
  ##However, a same genes may belongs to different data_params groups, therefore, we need to separate them
  ##For example, AAAS belongs to mRNA and others
  params_groups <- unlist(unique(rpolyA_combine_df$data_params))
  sub_grp=params_groups[1]
  rpolyA_subset_df<-subset(rpolyA_combine_df, data_params == sub_grp)
  rpolyA_subset_df$log_polyA_len=log(rpolyA_subset_df$polyA_len)
  
  rpolyA_sample_df=rpolyA_subset_df[, c("qname","polyA_len","log_polyA_len","gene_name")]
  rpolyA_sample_df$replicate=sample
  rpolyA_sample_df$condition=curr_condition
  
  logpolyA_combined_df=rbind(logpolyA_combined_df,rpolyA_sample_df)
}

write.table(logpolyA_combined_df,file=resdir%&%"/mRNA/logpolyA_read_based.txt", quote = FALSE, sep="\t", row.names =F, col.names =T)

input=read.table("/data/Nano_PolyA_Data.txt",header=T)
target_gene_list=input$Gene

## https://github.com/jchang97/mixedmodel/blob/main/tailfindr_mixed_model.R
results1<-NULL
results_r<-NULL
for (curr_gene in target_gene_list) {
  ##consider fixed effect:rand_age_weeks + rand_age + dem_gender+dem_weight+dem_ethnicity
  print(curr_gene)
   
  gene_logpolyA_df=logpolyA_combined_df[logpolyA_combined_df$gene_name==curr_gene,]
  fm <- lmer(log_polyA_len~ condition + (1 | replicate), data = gene_logpolyA_df)
  
  sum_df=as.data.frame(summary(fm)$coefficients)
  
  ## Anova-like table of random-effect terms using likelihood ratio tests:
  raov <-ranova(fm)
  # Extract p-value from the LRT for the random effect
  p_value_random_effect <- raov$`Pr(>Chisq)`[2:length(raov$`Pr(>Chisq)`)]
  # Extract p-values for fixed effects
  fixed_effects_p_values1 <- sum_df$`Pr(>|t|)`
  
  # Split the adjusted p-values back into fixed and random
  adjusted_fixed_effects_p_values1 <- p.adjust(fixed_effects_p_values1, method = "BH")
  adjusted_random_effect_p_value <- p.adjust(p_value_random_effect, method = "BH")
  
  raov=raov[-1,]
  sum_df$adjusted_P=adjusted_fixed_effects_p_values1 
  raov$adjusted_P=adjusted_random_effect_p_value
  sum_df$Gene=curr_gene
  sum_df=sum_df[-1,]
  raov$Gene=curr_gene
  
  results1<-rbind(results1,sum_df)
  results_r <- rbind(results_r, raov)
}

## when backtransform the effect-size, for positive value apply exp(log-value), for negative value apply (-1)*( 1/exp(log-value))
log_eff=results1$Estimate
backlog_eff <- ifelse(log_eff >= 0, exp(log_eff), -1 * (1 / exp(log_eff)))
results1$backlog_eff=backlog_eff 

write.csv(results1, resdir%&%"/mRNA/DES_lmer/fixed_term_summary.csv", quote = FALSE, row.names = TRUE)
write.csv(results_r, resdir%&%"/mRNA/DES_lmer/random_term_summary.csv", quote = FALSE, row.names = TRUE)

logFC <- results1$Estimate
p_value <- results1$adjusted_P
gene_names <- results1$Gene

# Create a dataframe for plotting
plot_data <- data.frame(
  gene = gene_names,
  logFC = logFC,
  p_value = p_value,
  negLogP = -log10(p_value),
  significance = case_when(
    p_value < 0.05 & abs(logFC) > 0.5 ~ "Significant",
    p_value >= 0.05 & abs(logFC) > 0.5 ~ "Marginally Significant",
    TRUE ~ "Not Significant"
  )
)

# Create the volcano plot using ggplot2
volcano_plot <- ggplot(plot_data, aes(x = logFC, y = negLogP, color = significance)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Significant" = "red2", "Marginally Significant" = "forestgreen", "Not Significant" = "grey30")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
  geom_text(data = subset(plot_data, significance == "Significant"), aes(label = gene), vjust = -1, hjust = 0.5,color="black",size = 2) +
  theme_minimal() +
  theme(
    legend.position = "none",  # Remove the legend
    axis.line = element_line(color = "black"),  # Add solid black lines for axes
    axis.title = element_text(size = 14),  # Increase size of axis titles
    axis.text = element_text(size = 12),  # Increase size of axis labels
    plot.title = element_text(size = 16, face = "bold"),  # Increase size of the plot title
    legend.text = element_text(size = 12)  
  ) +
  labs(
    x = "LogFold Change",
    y = "-log10(adjusted P-value)",
  )

# Display the plot
print(volcano_plot)

outdir=resdir%&%"/mRNA/DES_lmer/"
pdf(outdir%&%"VolcanoPlot_DP.pdf",width=5,height=4)
print(volcano_plot)
dev.off()
########## draw raincloud ############
sig_gene_list=results1[results1$adjusted_P<0.05,]$Gene
for (tmp_gene in sig_gene_list) {
  print(tmp_gene)
  tmp_df=logpolyA_combined_df[logpolyA_combined_df$gene_name==tmp_gene,]
  poly_df=tmp_df[, c("polyA_len","condition")]
  colnames(poly_df)=c("y_axis","x_axis")
  
  df_1x1 <- data_1x1( 
    array_1 = poly_df[poly_df$x_axis=="Bacteria",]$y_axis, #first set of values
    array_2 = poly_df[poly_df$x_axis=="Virus",]$y_axis, #second set of values
    jit_distance = .09,
    jit_seed = 321) 
  
  
  raincloud_1_h <- raincloud_1x1(
    data = df_1x1, 
    colors = (c('darkorange','darkgreen')), 
    fills = (c('darkorange','darkgreen')), 
    size = 0.1, 
    alpha = .5, 
    ort = 'h') +
    
    scale_x_continuous(breaks=c(1,2), labels=c("Bacteria", "Virus"), limits=c(0, 3)) +
    xlab("Conditions") + 
    ylab("PolyA Length (nt)") +
    labs(title = tmp_gene) +
    theme_classic(base_size = 12)+
    theme(plot.title = element_text(size = 15, face = "bold"),axis.text.x = element_text(angle = 45, vjust = 1, 
                                                                                         size = 10, hjust = 1),axis.text.y = element_text(size = 10))
  
  raincloud_1_h
  
  outdir=resdir%&%"/mRNA/DES_lmer/"
  pdf(outdir%&%tmp_gene%&%"_raincloudplots.pdf",width=4,height=3)
  print(raincloud_1_h)
  dev.off()
}
