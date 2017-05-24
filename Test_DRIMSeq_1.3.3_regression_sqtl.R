

library(devtools)
library(GeuvadisTranscriptExpr)
library(DRIMSeq)


##############################################################################
# Test arguments
##############################################################################


# rwd <- "/home/gosia/multinomial_project/package_devel/Test_DRIMSeq_1_3_3_regression"
# package_dir <- "/home/gosia/R/package_devel/DRIMSeq"


rwd <- "/Users/gosia/Dropbox/UZH/package_devel/Test_DRIMSeq_1_3_3_regression"
package_dir <- "/Users/gosia/Dropbox/UZH/package_devel/DRIMSeq"


##############################################################################
# Read in the arguments
##############################################################################

# rm(list = ls())
# 
# args <- (commandArgs(trailingOnly = TRUE))
# for (i in 1:length(args)) {
#   eval(parse(text = args[[i]]))
# }
# 
# cat(paste0(args, collapse = "\n"), fill = TRUE)

##############################################################################


dir.create(rwd, recursive = TRUE, showWarnings = FALSE)
setwd(rwd)


load_all(package_dir)

########################################################
# Examples
########################################################

#############################
### Create dmSQTLdata object
#############################

# Use subsets of data defined in GeuvadisTranscriptExpr package
library(GeuvadisTranscriptExpr)

counts <- GeuvadisTranscriptExpr::counts
genotypes <- GeuvadisTranscriptExpr::genotypes
gene_ranges <- GeuvadisTranscriptExpr::gene_ranges
snp_ranges <- GeuvadisTranscriptExpr::snp_ranges

colnames(counts)[c(1,2)] <- c("feature_id", "gene_id")
colnames(genotypes)[4] <- "snp_id"
samples <- data.frame(sample_id = colnames(counts)[-c(1,2)])

d <- dmSQTLdata(counts = counts, gene_ranges = gene_ranges,
  genotypes = genotypes, snp_ranges = snp_ranges, samples = samples,
  window = 5e3, BPPARAM = BiocParallel::SerialParam())

plotData(d)


data_dmSQTLdata <- d

#############################
### sQTL analysis
#############################
# If possible, use BPPARAM = BiocParallel::MulticoreParam() with more workers

d <- data_dmSQTLdata

head(names(d))
length(d)
d[1:10, ]
d[1:10, 1:10]

### Filtering
d <- dmFilter(d, min_samps_gene_expr = 70, min_samps_feature_expr = 5,
  min_samps_feature_prop = 0, minor_allele_freq = 5,
  BPPARAM = BiocParallel::SerialParam())

plotData(d)


### Calculate dispersion
d <- dmDispersion(d, BPPARAM = BiocParallel::SerialParam())

plotDispersion(d)

### Fit full model proportions
d <- dmFit(d, BPPARAM = BiocParallel::SerialParam())



### Fit null model proportions and test for sQTLs
d <- dmTest(d, BPPARAM = BiocParallel::SerialParam())



plotPValues(d)



head(results(d))





# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------

y <- x@counts[[1]]

genotype <- x@genotypes[[1]][1, ] 

df <- data.frame(y = 1, genotype = factor(genotype))

model.matrix(y ~ genotype, data = df)





genotype <- x@genotypes[[1]][11, ] 

df <- data.frame(y = 1, genotype = factor(genotype, levels = c("0", "1", "2")))

design <- model.matrix(y ~ genotype, data = df)




groups <- edgeR::designAsFactor(design)

nlevels(groups) == ncol(design)





genotype <- x@genotypes[[1]][11, ] 

df <- data.frame(y = 1, genotype = genotype + 1)

design <- model.matrix( ~ genotype, data = df)




groups <- edgeR::designAsFactor(design)

nlevels(groups) == ncol(design)


# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------

load_all(package_dir)

g <- 1
counts <- x@counts
genotypes <- x@genotypes
disp <- 10
prop_mode <- "constrOptim"
prop_tol <- 1e-12
verbose <- TRUE


coef_mode = "optim"
coef_tol = 1e-12
BPPARAM = BiocParallel::SerialParam()


x <- d
permutation_mode = "all_genes"
one_way = TRUE
prop_mode = "constrOptim"
prop_tol = 1e-12
coef_mode = "optim"
coef_tol = 1e-12
verbose = 0
BPPARAM = BiocParallel::SerialParam()

max_nr_perm_cycles = 10
max_nr_min_nr_sign_pval = 1e3

max_nr_perm = 1e6
max_nr_sign_pval = 1e2











