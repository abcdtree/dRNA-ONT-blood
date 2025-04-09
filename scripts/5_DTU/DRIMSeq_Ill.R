"%&%" <- function(a,b) paste(a,b, sep = "")

pasilla_metadata <- read.table("metadata_DRIMSeq.txt",
                               header = TRUE, as.is = TRUE)
## Load counts
pasilla_counts <- read.table("OUT.transcript_DRIMSeq_tpm.tsv",
                             header = TRUE, as.is = TRUE,check.names=FALSE)

library(DRIMSeq)
pasilla_samples <- data.frame(sample_id = pasilla_metadata$Sample,
                              group = pasilla_metadata$condition)
levels(pasilla_samples$group)
## NULL
d <- dmDSdata(counts = pasilla_counts, samples = pasilla_samples)
d

head(counts(d), 3)

head(samples(d), 3)

pdf(file="plotData.pdf",width = 10,height = 6)
plotData(d)
dev.off()

gene_id_subset <- readLines("gene_id_subset_DRIMSeq_All.txt")
d <- d[names(d) %in% gene_id_subset, ]
d
## An object of class dmDSdata
## with 42 genes and 7 samples
## * data accessors: counts(), samples()

# Check what is the minimal number of replicates per condition
table(samples(d)$group)
##
## CTL KD
## 4 3
d <- dmFilter(d, min_samps_gene_expr = 12, min_samps_feature_expr = 4,
              min_gene_expr = 10, min_feature_expr = 10)

## Create the design matrix
design_full <- model.matrix(~ group, data = samples(d))
design_full
## (Intercept) groupKD
## 1 1 0
## 2 1 0
## 3 1 0
## 4 1 1
## 5 1 1
## 6 1 1
## 7 1 0
## attr(,"assign")
## [1] 0 1
## attr(,"contrasts")
## attr(,"contrasts")$group
## [1] "contr.treatment"
## To make the analysis reproducible
set.seed(123)
## Calculate precision
d <- dmPrecision(d, design = design_full)
## ! Using a subset of 0.1 genes to estimate common precision !
## ! Using common_precision = 44.3452 as prec_init !
## ! Using 0 as a shrinkage factor !
d
## An object of class dmDSprecision
## with 26 genes and 7 samples
## * data accessors: counts(), samples()
## design()
## mean_expression(), common_precision(), genewise_precision()
head(mean_expression(d), 3)
## gene_id mean_expression
## 1 FBgn0000256 2622.286
## 2 FBgn0020309 13217.714
## 3 FBgn0259735 11992.903
common_precision(d)
## [1] 44.34515
head(genewise_precision(d))
## gene_id genewise_precision
## 1 FBgn0000256 93.388660
## 2 FBgn0020309 6.196615
## 3 FBgn0259735 72.027232
## 4 FBgn0032785 11.380309
## 5 FBgn0040297 13.085074
## 6 FBgn0032979 103.009092

pdf(file="plotPrecision.pdf",width = 10,height = 6)
plotPrecision(d)
dev.off()

library(ggplot2)
ggp <- plotPrecision(d)
ggp + geom_point(size = 4)

pdf(file="plotPrecision2.pdf",width = 10,height = 6)
ggp
dev.off()

d <- dmFit(d, design = design_full, verbose = 1)
## * Fitting the DM model..
## Using the one way approach.
## Took 0.1268 seconds.
## * Fitting the BB model..
## Using the one way approach.
## Took 0.0758 seconds.
d
## An object of class dmDSfit
## with 26 genes and 7 samples
## * data accessors: counts(), samples()
## design()
## mean_expression(), common_precision(), genewise_precision()
## proportions(), coefficients()
## Get fitted proportions
head(proportions(d))
## gene_id feature_id GSM461176 GSM461177 GSM461178 GSM461179
## 1 FBgn0000256 FBtr0290077 0.357375706 0.357375706 0.357375706 0.076553229
## 2 FBgn0000256 FBtr0290078 0.043422163 0.043422163 0.043422163 0.270030880
## 3 FBgn0000256 FBtr0290082 0.006049622 0.006049622 0.006049622 0.001885836
## 4 FBgn0000256 FBtr0077511 0.286678646 0.286678646 0.286678646 0.222608942
## 5 FBgn0000256 FBtr0290081 0.016852677 0.016852677 0.016852677 0.036981836
## 6 FBgn0000256 FBtr0077513 0.280207783 0.280207783 0.280207783 0.378704927
## GSM461180 GSM461181 GSM461182
## 1 0.076553229 0.076553229 0.357375706
## 2 0.270030880 0.270030880 0.043422163
## 3 0.001885836 0.001885836 0.006049622
## 4 0.222608942 0.222608942 0.286678646
## 5 0.036981836 0.036981836 0.016852677
## 6 0.378704927 0.378704927 0.280207783
## Get the DM regression coefficients (gene-level)
head(coefficients(d))
## gene_id feature_id X.Intercept. groupKD
## 1 FBgn0000256 FBtr0290077 3.6366530 -1.88148245
## 2 FBgn0000256 FBtr0290078 1.5288353 1.48688523
## 3 FBgn0000256 FBtr0290082 -0.4421389 -1.50630555
## 4 FBgn0000256 FBtr0077511 3.4162273 -0.59362640
## 5 FBgn0000256 FBtr0290081 0.5823749 0.44523627
## 6 FBgn0000256 FBtr0077513 3.3933968 -0.03945519
## Get the BB regression coefficients (feature-level)
head(coefficients(d), level = "feature")
## gene_id feature_id X.Intercept. groupKD
## 1 FBgn0000256 FBtr0290077 3.6366530 -1.88148245
## 2 FBgn0000256 FBtr0290078 1.5288353 1.48688523
## 3 FBgn0000256 FBtr0290082 -0.4421389 -1.50630555
## 4 FBgn0000256 FBtr0077511 3.4162273 -0.59362640
## 5 FBgn0000256 FBtr0290081 0.5823749 0.44523627
## 6 FBgn0000256 FBtr0077513 3.3933968 -0.039455


d <- dmTest(d, coef = "groupVirus", verbose = 1)
## * Fitting the DM model..
## Using the one way approach.
## Took 0.0789 seconds.
## * Calculating likelihood ratio statistics..
## Took 6e-04 seconds.
## * Fitting the BB model..
## Using the one way approach.
## Took 0.0313 seconds.
## * Calculating likelihood ratio statistics..
## Took 3e-04 seconds.
design(d)
## (Intercept)
## 1 1
## 2 1
## 3 1
## 4 1
## 5 1
## 6 1
## 7 1
head(results(d), 3)

## gene_id lr df pvalue adj_pvalue
## 1 FBgn0000256 146.1896729 6 4.942434e-29 1.285033e-27
## 2 FBgn0020309 17.9269288 4 1.275344e-03 4.467073e-03
## 3 FBgn0259735 0.9374648 2 6.257950e-01 7.747938e-01

design_null <- model.matrix(~ 1, data = samples(d))
design_null
## (Intercept)
## 1 1
## 2 1
## 3 1
## 4 1
## 5 1
## 6 1
## 7 1
## attr(,"assign")
## [1] 0
d <- dmTest(d, design = design_null)
head(results(d), 3)
## gene_id lr df pvalue adj_pvalue
## 1 FBgn0000256 146.1896729 6 4.942434e-29 1.285033e-27
## 2 FBgn0020309 17.9269288 4 1.275344e-03 4.467073e-03
## 3 FBgn0259735 0.9374648 2 6.257950e-01 7.747938e-01

contrast <- c(0, 1)
d <- dmTest(d, contrast = contrast)
design(d)
## [,1]
## 1 -1
## 2 -1
## 3 -1
## 4 -1
## 5 -1
## 6 -1
## 7 -1
head(results(d), 3)
## gene_id lr df pvalue adj_pvalue
## 1 FBgn0000256 146.1896729 6 4.942434e-29 1.285033e-27
## 2 FBgn0020309 17.9269288 4 1.275344e-03 4.467073e-03
## 3 FBgn0259735 0.9374648 2 6.257950e-01 7.747938e-01
head(results(d, level = "feature"), 3)
## gene_id feature_id lr df pvalue adj_pvalue
## 1 FBgn0000256 FBtr0290077 87.619933 1 7.932120e-21 4.521309e-19
## 2 FBgn0000256 FBtr0290078 72.329719 1 1.820846e-17 6.919213e-16
## 3 FBgn0000256 FBtr0290082 1.341238 1 2.468158e-01 5.542627e-01

plotPValues(d)

pdf(file="plotPValues.pdf",width = 5,height = 3)
plotPValues(d)
dev.off()

plotPValues(d, level = "feature")
pdf(file="plotPValues_featire.pdf",width = 5,height = 3)
plotPValues(d, level = "feature")
dev.off()

res <- results(d)
res <- res[order(res$pvalue, decreasing = FALSE), ]

saveRDS(d,"d.rds")
write.csv(res, "res.csv", row.names=FALSE, quote = FALSE)

res_sig=res[res$adj_pvalue<0.05,]
pdf("groupplot.pdf",width = 10,height = 6) 
for(p in 1:dim(res_sig)[1]){
  top_gene_id <- res$gene_id[p]
  print(plotProportions(d, gene_id = top_gene_id, group_variable = "group"))
}
dev.off()

pdf("lineplot.pdf",width = 10,height = 6) 
for(p in 1:dim(res_sig)[1]){
  top_gene_id <- res$gene_id[p]
  print(plotProportions(d, gene_id = top_gene_id, group_variable = "group",plot_type = "lineplot"))
}
dev.off()

pdf("ribbonplot.pdf",width = 10,height = 6) 
for(p in 1:dim(res_sig)[1]){
  top_gene_id <- res$gene_id[p]
  print(plotProportions(d, gene_id = top_gene_id, group_variable = "group",plot_type = "ribbonplot"))
}
dev.off()

## Coordinate system already present. Adding new coordinate system, which will
## replace the existing one.

library(stageR)

pScreen<-results(d)$pvalue
names(pScreen)<-results(d)$gene_id

## Assign transcript-level pvalues to the confirmation stage
pConfirmation<-matrix(results(d,level="feature")$pvalue,ncol=1)
rownames(pConfirmation)<-results(d,level="feature")$feature_id
## Create the gene-transcript mapping
tx2gene<-results(d,level="feature")[,c("feature_id","gene_id")]
## Create the stageRTx object and perform the stage-wise analysis
stageRObj<-stageRTx(pScreen= pScreen,pConfirmation= pConfirmation,
                    pScreenAdjusted=FALSE,tx2gene= tx2gene)
stageRObj<-stageWiseAdjustment(object= stageRObj,method="dtu",
                               alpha=0.05,allowNA=TRUE)
getSignificantGenes(stageRObj)
getSignificantTx(stageRObj)

#padj<-getAdjustedPValues(stageRObj,order=TRUE,onlySignificantGenes=FALSE)
#head(padj)
