

library(devtools)
library(PasillaTranscriptExpr)
library(DRIMSeq)
library(edgeR)

##############################################################################
# Test arguments
##############################################################################


rwd <- "/home/gosia/multinomial_project/package_devel/Test_DRIMSeq_1_3_3_regression"
package_dir <- "/home/gosia/R/package_devel/DRIMSeq"


# rwd <- "/Users/gosia/Dropbox/UZH/package_devel/Test_DRIMSeq_1_3_3_regression"
# package_dir <- "/Users/gosia/Dropbox/UZH/package_devel/DRIMSeq"


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

plotData(d)

# Use a subset of genes, which is defined in the following file
gene_id_subset <- readLines(file.path(data_dir, "gene_id_subset.txt"))
d <- d[names(d) %in% gene_id_subset, ]
d

plotData(d)

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


plotData(d)


# ---------------------------------------------------

counts = d@counts

design = as.matrix(model.matrix( ~ group, data = d@samples))
design


dispersion = 10
one_way = FALSE
prop_mode = "constrOptim"
prop_tol = 1e-12
coef_mode = "optim"
coef_tol = 1e-12
verbose = FALSE
BPPARAM = BiocParallel::SerialParam()




# ---------------------------------------------------

load_all(package_dir)

y <- counts[[1]]
coef_mode = "optim"
coef_tol = 1e-12
disp = 10


dm_fitRegression(y, design, 
  disp, coef_mode = "optim", coef_tol = 1e-12)
  

groups <- edgeR::designAsFactor(design)
groups <- factor(groups, labels = paste0("gr", levels(groups)))
ngroups <- nlevels(groups)
lgroups <- levels(groups)
igroups <- lapply(lgroups, function(gr){which(groups == gr)})
names(igroups) <- lgroups


dm_fitManyGroups(y, ngroups, lgroups, igroups, 
  disp, prop_mode = "constrOptim", prop_tol = 1e-12)


# ---------------------------------------------------

dmDS_estimateCommonDispersion(counts, design, disp_adjust = TRUE, 
  disp_interval = c(0, 1e+5), disp_tol = 1e-01, prop_mode = "constrOptim", 
  prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::SerialParam())
  

# ---------------------------------------------------


load_all(package_dir)

fit1 <- dmDS_fit(counts, design, dispersion, one_way = FALSE, 
  prop_mode = "constrOptim", prop_tol = 1e-12, 
  coef_mode = "optim", coef_tol = 1e-12,
  verbose = FALSE, 
  BPPARAM = BiocParallel::SerialParam())


y <- counts[[1]]
disp = 10
fit <- fit1$fit[[1]]


bb_fitRegression(y = y, design = design, disp = disp, fit = fit)


f1 <- bbDS_fit(counts = counts[1, ], fit = fit1$fit[1, ], design = design, dispersion = disp,
  one_way = TRUE, verbose = TRUE, BPPARAM = BiocParallel::SerialParam())


f2 <- bbDS_fit(counts = counts[1, ], fit = fit1$fit[1, ], design = design, dispersion = disp,
  one_way = FALSE, verbose = TRUE, BPPARAM = BiocParallel::SerialParam())


round(c(f1$coef[[1]]), 6) == round(c(f2$coef[[1]]), 6)

round(f1$lik, 6) == round(f2$lik, 6)

f1$lik == f2$lik


# ---------------------------------------------------

b <- b_init
x <- design

m_Hessian_regG(b, x, y)




# ---------------------------------------------------
# Run the DM pipeline
# ---------------------------------------------------

load_all(package_dir)

design = as.matrix(model.matrix( ~ group, data = samples(d)))
design


set.seed(1234)

### Calculate dispersion
d <- dmDispersion(d, design = design, common_dispersion = FALSE, one_way = FALSE, verbose = 1, 
  coef_mode = "optim",
  BPPARAM = BiocParallel::MulticoreParam(workers = 2))


plotDispersion(d)


head(mean_expression(d))
common_dispersion(d)
head(genewise_dispersion(d))


# ---------------------------------------------------


load_all(package_dir)

### Fit full model proportions
d <- dmFit(d, design = design, one_way = FALSE, verbose = 1, BPPARAM = BiocParallel::SerialParam())


head(proportions(d))
head(coefficients(d))

head(coefficients(d, level = "feature"))


x = d
gene_id = "FBgn0000256"
group_var = "group"
plot_type = "barplot" 
order = TRUE
plot_fit = TRUE
plot_main = TRUE


plotProportions(d, gene_id = gene_id, group_var = "group", group_colors = c("red", "blue"))


plotProportions(d, gene_id = gene_id, group_var = "group", plot_type = "lineplot")

plotProportions(d, gene_id = gene_id, group_var = "group", plot_type = "ribbonplot")

plotProportions(d, gene_id = gene_id, group_var = "group", plot_type = "boxplot1")

plotProportions(d, gene_id = gene_id, group_var = "group", plot_type = "boxplot2")


# ---------------------------------------------------


load_all(package_dir)

### Fit null model proportions and test for DS
d <- dmTest(d, coef = 2, one_way = FALSE, verbose = 1, BPPARAM = BiocParallel::SerialParam())


d <- dmTest(d, contrast = c(0, -1), verbose = 1, BPPARAM = BiocParallel::SerialParam())


x <- d
FDR = 0.05
verbose = 0
BPPARAM = BiocParallel::SerialParam()
pvalue_gene = x@results_gene[, c("gene_id", "pvalue"), drop = FALSE]
pvalue_feature = x@results_feature[, c("gene_id", "feature_id", "pvalue"), drop = FALSE]


dmTwoStageTest(d)




plotPValues(d)


head(proportions(d))
head(coefficients(d))
head(results(d))



































