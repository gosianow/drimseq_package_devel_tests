rwd <- "/home/gosia/multinomial_project/package_devel/Test_DRIMSeq_1_1_3_sQTL_permutations"
dir.create(rwd, recursive = TRUE)
setwd(rwd)

library(DRIMSeq)

# Use subsets of data defined in GeuvadisTranscriptExpr package
library(GeuvadisTranscriptExpr)

counts <- GeuvadisTranscriptExpr::counts
genotypes <- GeuvadisTranscriptExpr::genotypes
gene_ranges <- GeuvadisTranscriptExpr::gene_ranges
snp_ranges <- GeuvadisTranscriptExpr::snp_ranges

# Make sure that samples in counts and genotypes are in the same order
sample_id <- colnames(counts[, -(1:2)])

d <- dmSQTLdataFromRanges(counts = counts[, sample_id],
   gene_id = counts$Gene_Symbol, feature_id = counts$TargetID,
   gene_ranges = gene_ranges, genotypes = genotypes[, sample_id],
   snp_id = genotypes$snpId, snp_ranges = snp_ranges, sample_id = sample_id,
   window = 5e3, BPPARAM = BiocParallel::SerialParam())

plotData(d, out_dir = "./")



#############################
### sQTL analysis
#############################
# If possible, use BPPARAM = BiocParallel::MulticoreParam() with more workers


### Filtering
d <- dmFilter(d, min_samps_gene_expr = 70, min_samps_feature_expr = 5,
   min_samps_feature_prop = 0, minor_allele_freq = 5,
   BPPARAM = BiocParallel::MulticoreParam(workers = 10))
plotData(d, out_dir = "./")


### Calculate dispersion
d <- dmDispersion(d, mean_expression = TRUE, common_dispersion = FALSE, genewise_dispersion = TRUE, BPPARAM = BiocParallel::MulticoreParam(workers = 10))
plotDispersion(d, out_dir = "./")

### Fit full model proportions
d <- dmFit(d, BPPARAM = BiocParallel::MulticoreParam(workers = 10))



### Fit null model proportions and test for sQTLs
d <- dmTest(d, BPPARAM = BiocParallel::MulticoreParam(workers = 10))
plotTest(d, out_dir = "./")

head(results(d))



### Plot feature proportions for top sQTL
res <- results(d)
res <- res[order(res$pvalue, decreasing = FALSE), ]

gene_id <- res$gene_id[1]
snp_id <- res$snp_id[1]

plotFit(d, gene_id, snp_id, out_dir = "./gg")
plotFit(d, gene_id, snp_id, plot_type = "boxplot2", order = FALSE, out_dir = "./gg")
plotFit(d, gene_id, snp_id, plot_type = "ribbonplot", out_dir = "./gg")








