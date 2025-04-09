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

#write.table(logpolyA_combined_df,file=resdir%&%"/mRNA/logpolyA_read_based.txt", quote = FALSE, sep="\t", row.names =F, col.names =T)
logpolyA_combined_df = read.csv("/data/logpolyA_read_based.txt", header=T, sep="\t")
input=read.table("/data/Nano_PolyA_Data.txt",header=T)
target_gene_list=input$Gene

gene_count <- logpolyA_combined_df %>% count(gene_name)
p <- ggplot(gene_count_no_HBG, aes(y=n)) + 
  geom_boxplot() + scale_y_continuous(trans='log10')
p
gene_count_no_HBG <- tail(gene_count[order(gene_count$n, decreasing=TRUE), ], -3)

target_gene_list
gene_target <- filter(gene_count,
                       gene_name %in% target_gene_list)
gene_target

p <- ggplot(gene_target, aes(y=n)) + 
  geom_boxplot() +  scale_y_continuous(trans='log10')
p
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

results1

#plot
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

#plot_data[plot_data$significance == "Significant",]$gene
new_gene_list = plot_data[plot_data$significance == "Significant",]$gene


#redo with limited gene list
results1<-NULL
results_r<-NULL
for (curr_gene in new_gene_list) {
  ##consider fixed effect:rand_age_weeks + rand_age + dem_gender+dem_weight+dem_ethnicity
  print(curr_gene)
  
  
  gene_logpolyA_df=logpolyA_combined_df[logpolyA_combined_df$gene_name==curr_gene,]
  if(curr_gene == "RNF10"){
    print(gene_logpolyA_df %>% count(condition))
  }else{
    next
  }
  
  sample_length <- length(gene_logpolyA_df)
  for(n in 1:100){
    #sub_gene_df <- sample(gene_logpolyA_df[,c("log_polyA_len","condition", "replicate")], sample_length, replace=TRUE)
    #sub_gene_df <- gene_logpolyA_df[sample(nrow(gene_logpolyA_df), size=sample_length, replace=TRUE), ]
    ins = sample.int(nrow(gene_logpolyA_df), size=nrow(gene_logpolyA_df))
    sub_gene_df <- gene_logpolyA_df[inds,]
    while(1){
      if(nrow(unique(sub_gene_df$condition))==1){
        ins = sample.int(nrow(gene_logpolyA_df), size=nrow(gene_logpolyA_df))
        sub_gene_df <- gene_logpolyA_df[inds,]
      }
      else{
        break
      }
    }
    fm <- lmer(log_polyA_len~ condition + (1 | replicate), data = sub_gene_df)
    
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
    sum_df$Run_ID = n
    sum_df=sum_df[-1,]
    raov$Gene=curr_gene
    
    results1<-rbind(results1,sum_df)
    results_r <- rbind(results_r, raov)
    
  }
  
}

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

#outdir=resdir%&%"/mRNA/DES_lmer/"
#pdf(outdir%&%"VolcanoPlot_DP.pdf",width=5,height=4)
#print(volcano_plot)
#dev.off()

