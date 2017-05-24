### Set up R session
library('DRIMSeq')
library('PasillaTranscriptExpr')
data_dir  <- system.file("extdata", package = "PasillaTranscriptExpr")
pasilla_metadata <- read.table(file.path(data_dir, "metadata.txt"), header = TRUE, as.is = TRUE)
pasilla_counts <- read.table(file.path(data_dir, "counts.txt"), header = TRUE, as.is = TRUE)

pasilla_samples <- data.frame(sample_id = pasilla_metadata$SampleName, group = pasilla_metadata$condition)
levels(pasilla_samples$group)

d <- dmDSdata(counts = pasilla_counts, samples = pasilla_samples)

gene_id_subset <- readLines(file.path(data_dir, "gene_id_subset.txt"))
d <- d[names(d) %in% gene_id_subset, ]

design_full <- model.matrix(~ group, data = samples(d))
design_full


### With filtering
d1 <- dmFilter(d, min_samps_gene_expr = 7, min_samps_feature_expr = 3, min_gene_expr = 10, min_feature_expr = 10)

d1 <- dmPrecision(d1, design = design_full)

d1 <- dmFit(d1, design = design_full, bb_model=TRUE)



### Using one way approach

d2 <- dmFilter(d)

d2 <- dmPrecision(d2, design = design_full)

genewise_precision(d2)

d2 <- dmFit(d2, design = design_full, bb_model=TRUE)

## Get fitted proportions

p2 <- proportions(d2)

c2 <- counts(d2)

c2[c2$gene_id %in% p2[is.na(p2[, 3]), "gene_id"], ]



### Using regression

d3 <- dmFilter(d)

d3 <- dmPrecision(d3, design = design_full, one_way = FALSE)

genewise_precision(d3)

d3 <- dmFit(d3, design = design_full, bb_model = TRUE, one_way = FALSE)


## Get fitted proportions

p3 <- proportions(d3)

c3 <- counts(d3)

c3[c3$gene_id %in% p3[is.na(p3[, 3]), "gene_id"], ]




### Using one way approach + some filtering (but with min_samps_feature_expr = 1:3, there is still one gene with NAs)

d2 <- dmFilter(d, min_samps_feature_expr = 4, min_feature_expr = 1)

d2 <- dmPrecision(d2, design = design_full)

genewise_precision(d2)

d2 <- dmFit(d2, design = design_full, bb_model=TRUE)

## Get fitted proportions

p2 <- proportions(d2)

c2 <- counts(d2)

c2[c2$gene_id %in% p2[is.na(p2[, 3]), "gene_id"], ]




### Using regression + some filtering (but with min_samps_feature_expr = 1:2, there is still one gene with NAs)

d3 <- dmFilter(d, min_samps_feature_expr = 3, min_feature_expr = 1)

d3 <- dmPrecision(d3, design = design_full, one_way = FALSE)

genewise_precision(d3)

d3 <- dmFit(d3, design = design_full, bb_model = TRUE, one_way = FALSE)


## Get fitted proportions

p3 <- proportions(d3)

c3 <- counts(d3)

c3[c3$gene_id %in% p3[is.na(p3[, 3]), "gene_id"], ]








