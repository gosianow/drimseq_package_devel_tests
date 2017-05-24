### Load libraries
library('DRIMSeq')
library('PasillaTranscriptExpr')

### Lets used the pasilla data
data_dir  <- system.file("extdata", package = "PasillaTranscriptExpr")
pasilla_metadata <- read.table(file.path(data_dir, "metadata.txt"), header = TRUE, as.is = TRUE)
pasilla_counts <- read.table(file.path(data_dir, "counts.txt"), header = TRUE, as.is = TRUE)
gene_id_subset <- readLines(file.path(data_dir, "gene_id_subset.txt"))



### pasilla_samples design matrix with 3 groups
# basicly splitting ctrl into two

pasilla_samples <- data.frame(sample_id = pasilla_metadata$SampleName, group = pasilla_metadata$condition)

pasilla_samples <- data.frame(sample_id = pasilla_metadata$SampleName, group = c('Nr1','Nr1','Nr2','Nr3','Nr3','Nr3','Nr2'))

levels(pasilla_samples$group)


### Create and subset dmData
d <- dmDSdata(counts = pasilla_counts, samples = pasilla_samples)
d <- d[names(d) %in% gene_id_subset, ]
### Filtering
d <- dmFilter(d, min_samps_gene_expr = 7, min_samps_feature_expr = 3, min_gene_expr = 10, min_feature_expr = 10)



## Create the design matrix without intercept due to no clear 'ground state'

design_full2 <- model.matrix(~ -1 + group, data = pasilla_samples)
design_full2

### Calculate precision and fit
d2 <- dmPrecision(d, design = design_full2)

genewise_precision(d2)

d2 <- dmFit(d2, design = design_full2, bb_model = TRUE)

head(proportions(d2))



contrast1 <- c(-1, 1, 0)

t21 <- dmTest(d2, contrast = contrast1, verbose = TRUE)

results(t21)



contrast2 <- c(-1, 0, 1)

t22 <- dmTest(d2, contrast = contrast2, verbose = TRUE)

results(t22)



contrast3 <- c(0, -1, 1)

t23 <- dmTest(d2, contrast = contrast3, verbose = TRUE)

results(t23)




## Create the design matrix with intercept

design_full1 <- model.matrix(~ group, data = pasilla_samples)





