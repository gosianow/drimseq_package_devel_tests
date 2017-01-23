

library(devtools)
library(PasillaTranscriptExpr)
library(DRIMSeq)
library(edgeR)

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


# ---------------------------------------------------

load_all(package_dir)


counts = d@counts

design = as.matrix(model.matrix( ~ group, data = d@samples))
design

dispersion = 10
prop_mode = "constrOptim"
prop_tol = 1e-12
verbose = FALSE
BPPARAM = BiocParallel::SerialParam()


y <- counts[[1]]
coef_mode = "optim"
coef_tol = 1e-12
disp = 10


dm_fitRegression(y, design, 
  disp, coef_mode = "optim", coef_tol = 1e-12, verbose = FALSE)
  

groups <- edgeR::designAsFactor(design)
groups <- factor(groups, labels = paste0("gr", levels(groups)))
ngroups <- nlevels(groups)
lgroups <- levels(groups)
igroups <- lapply(lgroups, function(gr){which(groups == gr)})
names(igroups) <- lgroups


dm_fitManyGroups(y, ngroups, lgroups, igroups, 
  disp, prop_mode = "constrOptim", prop_tol = 1e-12, verbose = FALSE)


dmDS_fit(counts, design, dispersion,
  prop_mode = "constrOptim", prop_tol = 1e-12, verbose = FALSE, 
  BPPARAM = BiocParallel::SerialParam())



dmDS_estimateCommonDispersion(counts, design, disp_adjust = TRUE, 
  disp_interval = c(0, 1e+5), disp_tol = 1e-01, prop_mode = "constrOptim", 
  prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::SerialParam())
  
# ---------------------------------------------------

b <- b_init
x <- design

m_Hessian_regG(b, x, y)





# ---------------------------------------------------

load_all(package_dir)


### Calculate dispersion
d <- dmDispersion(d, design = design, verbose = 1, 
  BPPARAM = BiocParallel::MulticoreParam(workers = 2))


plotDispersion(d, out_dir = "./")

head(mean_expression(d))
common_dispersion(d)
head(genewise_dispersion(d))


# ---------------------------------------------------


load_all(package_dir)

### Fit full model proportions
d <- dmFit(d, design = design, verbose = 1, BPPARAM = BiocParallel::SerialParam())


head(proportions(d))
head(statistics(d))





# ---------------------------------------------------


load_all(package_dir)

### Fit null model proportions and test for DS
d <- dmTest(d, coef = 2, verbose = 1, BPPARAM = BiocParallel::SerialParam())


d <- dmTest(d, contrast = c(0, -1), verbose = 1, BPPARAM = BiocParallel::SerialParam())


dmTwoStageTest(d)




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




































