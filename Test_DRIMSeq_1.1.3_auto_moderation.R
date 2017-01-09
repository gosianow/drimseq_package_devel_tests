

rwd <- "/home/gosia/multinomial_project/package_devel/Test_DRIMSeq_1_1_3_auto_moderation"
dir.create(rwd, recursive = TRUE)
setwd(rwd)

library(DRIMSeq)

#############################
### Create dmDSdata object
#############################
### Get kallisto transcript counts from 'PasillaTranscriptExpr' package

library(PasillaTranscriptExpr)

data_dir  <- system.file("extdata", package = "PasillaTranscriptExpr")

metadata <- read.table(file.path(data_dir, "metadata.txt"),
 header = TRUE, as.is = TRUE)
metadata

counts <- read.table(file.path(data_dir, "counts.txt"),
 header = TRUE, as.is = TRUE)
head(counts)

# Create a dmDSdata object
d <- dmDSdata(counts = counts[, metadata$SampleName],
 gene_id = counts$gene_id, feature_id = counts$feature_id,
 sample_id = metadata$SampleName, group = metadata$condition)

plotData(d, out_dir = "./")

# # Use a subset of genes, which is defined in the following file
# gene_id_subset <- readLines(file.path(data_dir, "gene_id_subset.txt"))
# d <- d[names(d) %in% gene_id_subset, ]
# 
# plotData(d)


###################################
### Differential splicing analysis
###################################
# If possible, use BPPARAM = BiocParallel::MulticoreParam() with more workers

### Filtering
# Check what is the minimal number of replicates per condition
table(samples(d)$group)

d <- dmFilter(d, min_samps_gene_expr = 7, min_samps_feature_expr = 3,
 min_samps_feature_prop = 1, min_feature_prop = 0.05)
plotData(d, out_dir = "./")

### Calculate dispersion
d <- dmDispersion(d, mean_expression = TRUE, common_dispersion = FALSE, genewise_dispersion = TRUE, BPPARAM = BiocParallel::MulticoreParam(workers = 10))
plotDispersion(d, out_dir = "./")


d <- dmDispersion(d, disp_moderation = "none", mean_expression = TRUE, common_dispersion = FALSE, genewise_dispersion = TRUE, BPPARAM = BiocParallel::MulticoreParam(workers = 10))
plotDispersion(d, out_dir = "./none_")



head(mean_expression(d))
common_dispersion(d)
head(genewise_dispersion(d))

### Fit full model proportions
d <- dmFit(d, BPPARAM = BiocParallel::SerialParam())

head(proportions(d))
head(statistics(d))

### Fit null model proportions and test for DS
d <- dmTest(d, BPPARAM = BiocParallel::SerialParam())
plotTest(d)

head(proportions(d))
head(statistics(d))
head(results(d))

### Plot feature proportions for top DS gene
res <- results(d)
res <- res[order(res$pvalue, decreasing = FALSE), ]

gene_id <- res$gene_id[1]

plotFit(d, gene_id = gene_id)
plotFit(d, gene_id = gene_id, plot_type = "lineplot")
plotFit(d, gene_id = gene_id, plot_type = "ribbonplot")
















