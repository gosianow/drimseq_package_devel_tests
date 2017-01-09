

library(devtools)
library(PasillaTranscriptExpr)
library(DRIMSeq)


##############################################################################
# Test arguments
##############################################################################


rwd <- "/home/gosia/multinomial_project/package_devel/Test_DRIMSeq_1_3_2_beta_binomial"
package_dir <- "/home/gosia/R/package_devel/DRIMSeq"



##############################################################################
# Read in the arguments
##############################################################################

rm(list = ls())

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

cat(paste0(args, collapse = "\n"), fill = TRUE)

##############################################################################


dir.create(rwd, recursive = TRUE, showWarnings = FALSE)
setwd(rwd)


load_all(package_dir)

########################################################
# Examples
########################################################


#############################
### Create dmDSdata object
#############################
### Get kallisto transcript counts from 'PasillaTranscriptExpr' package


data_dir  <- system.file("extdata", package = "PasillaTranscriptExpr")

# Load metadata
pasilla_metadata <- read.table(file.path(data_dir, "metadata.txt"), header = TRUE, 
  as.is = TRUE)

# Load counts
pasilla_counts <- read.table(file.path(data_dir, "counts.txt"), header = TRUE, 
  as.is = TRUE)

# Create a samples data frame
pasilla_samples <- data.frame(sample_id = pasilla_metadata$SampleName, 
  group = pasilla_metadata$condition)

# Create a dmDSdata object
d <- dmDSdata(counts = pasilla_counts, samples = pasilla_samples)
d

plotData(d, out_dir = "./")

# Use a subset of genes, which is defined in the following file
gene_id_subset <- readLines(file.path(data_dir, "gene_id_subset.txt"))
d <- d[names(d) %in% gene_id_subset, ]
d

plotData(d, out_dir = "./")

data_dmDSdata <- d


###################################
### Differential splicing analysis
###################################
# If possible, use BPPARAM = BiocParallel::MulticoreParam() with more workers

d <- data_dmDSdata

head(counts(d))
samples(d)
head(names(d))
length(d)
d[1:20, ]
d[1:20, 1:3]

### Filtering
# Check what is the minimal number of replicates per condition 
table(samples(d)$group)
d <- dmFilter(d, min_samps_gene_expr = 7, min_samps_feature_expr = 3, 
  min_samps_feature_prop = 0)

plotData(d, "./")

### Calculate dispersion
d <- dmDispersion(d, verbose = 1, BPPARAM = BiocParallel::MulticoreParam(workers = 2))

plotDispersion(d, out_dir = "./")

head(mean_expression(d))
common_dispersion(d)
head(genewise_dispersion(d))





### Fit full model proportions
d <- dmFit(d, verbose = 0, BPPARAM = BiocParallel::SerialParam())




head(proportions(d))


head(statistics(d))




load_all(package_dir)



### Fit null model proportions and test for DS
d <- dmTest(d, verbose = 1, BPPARAM = BiocParallel::SerialParam())



plotTest(d, out_dir = "./")


head(proportions(d))

head(statistics(d))

head(results(d))





### Plot feature proportions for top DS gene
res <- results(d)
res <- res[order(res$pvalue, decreasing = FALSE), ]

gene_id <- res$gene_id[1]

plotFit(d, gene_id = gene_id, out_dir = "./")
plotFit(d, gene_id = gene_id, plot_type = "lineplot", out_dir = "./")
plotFit(d, gene_id = gene_id, plot_type = "ribbonplot", out_dir = "./")



pvalue_gene = d@results_gene 
pvalue_feature = d@results_feature
FDR = 0.05
verbose = FALSE


table <- dmDS_two_stage_test(pvalue_gene, pvalue_feature, FDR = 0.05, 
  verbose = FALSE)


table <- dmTwoStageTest(d)

































