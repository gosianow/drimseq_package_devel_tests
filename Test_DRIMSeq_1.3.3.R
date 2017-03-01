
library(DRIMSeq)

library(devtools)
load_all()


### EXAMPLES ###
# -----------------------------------------------------------------------------
# Create dmDSdata object 
# -----------------------------------------------------------------------------
## Get kallisto transcript counts from the 'PasillaTranscriptExpr' package

library(PasillaTranscriptExpr)

data_dir  <- system.file("extdata", package = "PasillaTranscriptExpr")

## Load metadata
pasilla_metadata <- read.table(file.path(data_dir, "metadata.txt"), 
	header = TRUE, as.is = TRUE)

## Load counts
pasilla_counts <- read.table(file.path(data_dir, "counts.txt"), 
	header = TRUE, as.is = TRUE)

## Create a pasilla_samples data frame
pasilla_samples <- data.frame(sample_id = pasilla_metadata$SampleName, 
  group = pasilla_metadata$condition)
levels(pasilla_samples$group)

## Create a dmDSdata object
d <- dmDSdata(counts = pasilla_counts, samples = pasilla_samples)

## Use a subset of genes, which is defined in the following file
gene_id_subset <- readLines(file.path(data_dir, "gene_id_subset.txt"))

d <- d[names(d) %in% gene_id_subset, ]


# -----------------------------------------------------------------------------
# Differential transcript usage analysis - simple two group comparison 
# -----------------------------------------------------------------------------

## Filtering
## Check what is the minimal number of replicates per condition 
table(samples(d)$group)

d <- dmFilter(d, min_samps_gene_expr = 7, min_samps_feature_expr = 3,
  min_gene_expr = 10, min_feature_expr = 10)

plotData(d)

## Create the design matrix
design_full <- model.matrix(~ group, data = samples(d))

## To make the analysis reproducible
set.seed(123)
## Calculate precision
d <- dmPrecision(d, design = design_full, one_way = FALSE)

set.seed(123)
d1 <- dmPrecision(d, design = design_full)

plotPrecision(d)

head(mean_expression(d))

common_precision(d)
genewise_precision(d)

common_precision(d1)
genewise_precision(d1)


## Fit full model proportions
d <- dmFit(d, design = design_full)

## Get fitted proportions
head(proportions(d))
## Get the DM regression coefficients (gene-level) 
head(coefficients(d))
## Get the BB regression coefficients (feature-level) 
head(coefficients(d), level = "feature")

## Fit null model proportions and perform the LR test to detect DTU
d <- dmTest(d, coef = "groupKD")
design(d)
head(results(d))


## Plot the gene-level p-values
plotPValues(d)

## Get the gene-level results
head(results(d))

## Plot feature proportions for a top DTU gene
res <- results(d)
res <- res[order(res$pvalue, decreasing = FALSE), ]

top_gene_id <- res$gene_id[1]


plotProportions(d, gene_id = top_gene_id, group_variable = "group")

plotProportions(d, gene_id = top_gene_id, group_variable = "group", 
  plot_type = "lineplot")

plotProportions(d, gene_id = top_gene_id, group_variable = "group", 
  plot_type = "ribbonplot")

plotProportions(d, gene_id = top_gene_id, group_variable = "group", 
  plot_type = "boxplot1")

plotProportions(d, gene_id = top_gene_id, group_variable = "group", 
  plot_type = "boxplot2")



plotProportions(d, gene_id = top_gene_id, group_variable = "sample_id")

plotProportions(d, gene_id = top_gene_id, group_variable = "sample_id", 
  plot_type = "lineplot")

plotProportions(d, gene_id = top_gene_id, group_variable = "sample_id", 
  plot_type = "ribbonplot")

plotProportions(d, gene_id = top_gene_id, group_variable = "sample_id", 
  plot_type = "boxplot1")

plotProportions(d, gene_id = top_gene_id, group_variable = "sample_id", 
  plot_type = "boxplot2")

# ------------------------------------------

library(stageR)

pScreen <- results(d)$pvalue
names(pScreen) <- results(d)$gene_id

pConfirmation <- matrix(results(d, level = "feature")$pvalue, ncol = 1)
rownames(pConfirmation) <- results(d, level = "feature")$feature_id

tx2gene <- results(d, level = "feature")[, c("feature_id", "gene_id")]

stageRObj <- stageRTx(pScreen = pScreen, pConfirmation = pConfirmation, 
	pScreenAdjusted = FALSE, tx2gene = tx2gene)

stageRObj <- stageWiseAdjustment(object = stageRObj, method = "dtu", 
	alpha = 0.05)

getSignificantGenes(stageRObj)

getSignificantTx(stageRObj)

padj <- getAdjustedPValues(stageRObj, order = TRUE, 
	onlySignificantGenes = FALSE)

head(padj)



# ------------------------------------------

pasilla_samples2 <- data.frame(sample_id = pasilla_metadata$SampleName, 
  group = pasilla_metadata$condition, 
  library_layout = pasilla_metadata$LibraryLayout)

d2 <- dmDSdata(counts = pasilla_counts, samples = pasilla_samples2)

## Subsetting to a vignette runnable size
d2 <- d2[names(d2) %in% gene_id_subset, ]

## Filtering
d2 <- dmFilter(d2, min_samps_gene_expr = 7, min_samps_feature_expr = 3,
  min_gene_expr = 10, min_feature_expr = 10)

## Create the design matrix
design_full2 <- model.matrix(~ group + library_layout, data = samples(d2))
design_full2

## To make the analysis reproducible
set.seed(123)

## Calculate precision
d2 <- dmPrecision(d2, design = design_full2)

common_precision(d2)
genewise_precision(d2)

## Fit proportions
d2 <- dmFit(d2, design = design_full2, verbose = 1)

## Test for DTU
d2 <- dmTest(d2, coef = "groupKD", verbose = 1)
design(d2)
head(results(d2), 3)

## Plot p-value distribution
plotPValues(d2)

## Plot the top significant gene
res2 <- results(d2)
res2 <- res2[order(res2$pvalue, decreasing = FALSE), ]
top_gene_id2 <- res2$gene_id[1]


pdf("prop.pdf")
plotProportions(d2, gene_id = top_gene_id2, group_variable = "group")
dev.off()

x <- d
gene_id <- top_gene_id
group_variable  <- "group"
plot_type = "barplot"
order = TRUE
plot_fit = TRUE
plot_main = TRUE
group_colors = NULL
feature_colors = NULL






# -----------------------------------------------------------------------------
# Create dmSQTLdata object
# -----------------------------------------------------------------------------
# Use subsets of data defined in the GeuvadisTranscriptExpr package

library(GeuvadisTranscriptExpr)

geuv_counts <- GeuvadisTranscriptExpr::counts
geuv_genotypes <- GeuvadisTranscriptExpr::genotypes
geuv_gene_ranges <- GeuvadisTranscriptExpr::gene_ranges
geuv_snp_ranges <- GeuvadisTranscriptExpr::snp_ranges

colnames(geuv_counts)[c(1,2)] <- c("feature_id", "gene_id")
colnames(geuv_genotypes)[4] <- "snp_id"
geuv_samples <- data.frame(sample_id = colnames(geuv_counts)[-c(1,2)])

d <- dmSQTLdata(counts = geuv_counts, gene_ranges = geuv_gene_ranges,
  genotypes = geuv_genotypes, snp_ranges = geuv_snp_ranges, 
  samples = geuv_samples, window = 5e3)

# -----------------------------------------------------------------------------
# sQTL analysis - simple group comparison
# -----------------------------------------------------------------------------

## Filtering
d <- dmFilter(d, min_samps_gene_expr = 70, min_samps_feature_expr = 5,
  minor_allele_freq = 5, min_gene_expr = 10, min_feature_expr = 10)

plotData(d)

## To make the analysis reproducible
set.seed(123)
## Calculate precision
d <- dmPrecision(d)

plotPrecision(d)

## Fit full model proportions
d <- dmFit(d)

## Fit null model proportions, perform the LR test to detect tuQTLs 
## and use the permutation approach to adjust the p-values
d <- dmTest(d)

## Plot the gene-level p-values
plotPValues(d)

## Get the gene-level results
head(results(d))

res <- results(d)
res <- res[order(res$pvalue, decreasing = FALSE), ]

top_gene_id <- res$gene_id[1]
top_snp_id <- res$snp_id[1]


library(devtools)
load_all()


pdf("props.pdf")
plotProportions(d, gene_id = top_gene_id, snp_id = top_snp_id)

plotProportions(d, gene_id = top_gene_id, snp_id = top_snp_id,
  plot_type = "boxplot2", order = FALSE)

dev.off()














