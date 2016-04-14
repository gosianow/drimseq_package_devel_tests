

# R32dev
# Created 7 Apr 2016


rwd <- "/home/gosia/multinomial_project/package_devel/Test_DRIMSeq_0.3.3_sQTL_permutations"
dir.create(rwd, recursive = TRUE)

setwd(rwd)



library(DRIMSeq)
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)



### Load sQTL data 

data_dir  <- system.file("extdata", package = "DRIMSeq")
# gene_ranges with names!
gene_ranges <- import(paste0(data_dir, "/genes_subset.bed"))
names(gene_ranges) <- mcols(gene_ranges)$name
counts <- read.table(paste0(data_dir, "/TrQuantCount_CEU_subset.tsv"),
  header = TRUE, sep = "\t", as.is = TRUE)
genotypes <- read.table(paste0(data_dir, "/genotypes_CEU_subset.tsv"),
  header = TRUE, sep = "\t", as.is = TRUE)
# snp_ranges with names!
snp_ranges <- GRanges(Rle(genotypes$chr), IRanges(genotypes$start,
  genotypes$end))
names(snp_ranges) <- genotypes$snpId
## Check if samples in count and genotypes are in the same order
all(colnames(counts[, -(1:2)]) == colnames(genotypes[, -(1:4)]))
## [1] TRUE
sample_id <- colnames(counts[, -(1:2)])

d <- dmSQTLdataFromRanges(counts = counts[, -(1:2)], gene_id = counts$geneId,
  feature_id = counts$trId, gene_ranges = gene_ranges,
  genotypes = genotypes[, -(1:4)], snp_id = genotypes$snpId,
  snp_ranges = snp_ranges, sample_id = sample_id, window = 5e3,
  BPPARAM = BiocParallel::MulticoreParam(workers = 1))
d



d <- dmFilter(d, min_samps_gene_expr = 70, min_samps_feature_expr = 5,
  min_samps_feature_prop = 5, minor_allele_freq = 5,
  BPPARAM = BiocParallel::MulticoreParam(workers = 1))



d1 <- dmDispersion(d, common_dispersion = FALSE, speed = TRUE, BPPARAM = BiocParallel::MulticoreParam(workers = 5))

plotDispersion(d1, out_dir = "speed_TRUE_")

ggp1 <- plotDispersion(d1)



###########################################################################
### Comapre dispersion when speed = TRUE and speed = FALSE
###########################################################################

d2 <- dmDispersion(d, common_dispersion = FALSE, speed = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = 5))

plotDispersion(d2, out_dir = "speed_FALSE_")

ggp2 <- plotDispersion(d2)



### Plot both on one figure

ggp3 <- ggp2 +
  geom_point(data = ggp1$data, aes(x = mean_expression, y = dispersion), color = "black", size = 0.5)

pdf("speed_both_dispersion_vs_mean.pdf")
print(ggp3)
dev.off()



### Expectation that there are NAs for the cases where null theta is higher than full theta
### It is not the case.

disp1 <- unlist(d1@genewise_dispersion)
disp2 <- unlist(d2@genewise_dispersion)


table(disp1 <= disp2)
table(disp1 <= disp2, rowSums(is.na(d1@genotypes@unlistData)) == 0)





###########################################################################
### Implement permutations based on all genes
###########################################################################


library(devtools)
load_all("/home/gosia/R/package_devel/DRIMSeq")


d <- d1

d <- dmFit(d, BPPARAM = BiocParallel::MulticoreParam(workers = 5))


d <- dmTest(d, permutations = "all_genes", verbose = TRUE, BPPARAM = BiocParallel::MulticoreParam(workers = 10))

d <- dmTest(d, permutations = "per_gene", verbose = TRUE, BPPARAM = BiocParallel::MulticoreParam(workers = 10))


plotTest(d)




x <- d
test = "lr"
prop_mode = "constrOptimG"
prop_tol = 1e-12
verbose = TRUE 
BPPARAM = BiocParallel::MulticoreParam(workers = 10)


fit_null <- dmSQTL_fitOneModel(counts = x@counts, genotypes = x@genotypes, dispersion = slot(x, x@dispersion), model = "null", prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)


### Number of samples used for the analysis
n <- unlist(lapply(1:length(x@counts), function(g){
  sum(!is.na(x@counts[[g]][1, ]))
}))


results <- dmSQTL_test(fit_full = x@fit_full, fit_null = fit_null, test = test, n = n, verbose = verbose, BPPARAM = BPPARAM)




# library(gtools)
# gtools::permutations(4, 4)


max_nr_perm_cycles <- 10
max_nr_min_nr_sign_pval <- 1e3


dmSQTL_permutations_all_genes <- function(x, fit_null, results, max_nr_perm_cycles = 10, max_nr_min_nr_sign_pval = 1e3, prop_mode, prop_tol, n, test, verbose, BPPARAM){
  
  fit_full <- x@fit_full
  
  nr_perm_tot <- 0
  nr_perm_cycles <- 0
  min_nr_sign_pval <- 0
  n <- ncol(x@counts)
  
  
  pval <- results$pvalue
  nas <- is.na(pval)
  pval <- pval[!nas]
  pval <- factor(pval)
  sum_sign_pval <- rep(0, length(pval))
  
  # ds_genes <- results$adj_pvalue < 0.1
  
  while(nr_perm_cycles < max_nr_perm_cycles && min_nr_sign_pval < max_nr_min_nr_sign_pval){
    
    permutation <- sample(n, n)
    
    ### Permute counts for all genes
    counts <- x@counts[, permutation, drop = FALSE]
    
    fit_full_perm <- dmSQTL_fitOneModel(counts = counts, genotypes = x@genotypes, dispersion = slot(x, x@dispersion), model = "full", prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
    
    fit_null_perm <- dmSQTL_fitOneModel(counts = counts, genotypes = x@genotypes, dispersion = slot(x, x@dispersion), model = "null", prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
    
    results_perm <- dmSQTL_test(fit_full = fit_full_perm, fit_null = fit_null_perm, test = test, n = n, verbose = verbose, BPPARAM = BPPARAM)
    
    nr_perm <- nrow(results_perm)
    nr_perm_tot <- nr_perm_tot + nr_perm
    nr_perm_cycles <- nr_perm_cycles + 1
    
    
    ### Count how many pval_permuted is lower than pval from the model
    pval_perm <- results_perm$pvalue
    nas_perm <- is.na(pval_perm)
    pval_perm <- pval_perm[!nas_perm]
    pval_perm_cut <- cut(pval_perm, c(-1, levels(pval), 2), right=FALSE)
    pval_perm_sum <- table(pval_perm_cut)
    pval_perm_cumsum <- cumsum(pval_perm_sum)[-length(pval_perm_sum)]
    names(pval_perm_cumsum) <- levels(pval)
    sum_sign_pval <- sum_sign_pval + pval_perm_cumsum[pval]
    
    pval_adj <- (sum_sign_pval + 1) / (nr_perm_tot + 1)
    
    min_nr_sign_pval <- min(sum_sign_pval)
    
    
  }
  
  
  pval_out <- rep(NA, nrow(results))
  pval_out[!nas] <- pval_adj
  
  return(pval_out)
  
}





pval_permutations <- dmSQTL_permutations_all_genes(x, fit_null, results, max_nr_perm_cycles = 3, max_nr_min_nr_sign_pval = 1e3, prop_mode, prop_tol, verbose, n, test, BPPARAM)



ggp <- dm_plotPvalues(pvalues = results$pvalue)

pdf(paste0("hist_pvalues_dm.pdf"))
print(ggp)
dev.off()



ggp <- dm_plotPvalues(pvalues = pval_permutations)

pdf(paste0("hist_pvalues_permutations_all_genes.pdf"))
print(ggp)
dev.off()




###########################################################################
### Implement permutations foe each gene individually 
###########################################################################


library(devtools)
load_all("/home/gosia/R/package_devel/DRIMSeq")


d <- d1

d <- dmFit(d, BPPARAM = BiocParallel::MulticoreParam(workers = 5))


# d <- dmTest(d, BPPARAM = BiocParallel::MulticoreParam(workers = 5))


x <- d
test = "lr"
prop_mode = "constrOptimG"
prop_tol = 1e-12
verbose = 0 
BPPARAM = BiocParallel::MulticoreParam(workers = 5)


fit_null <- dmSQTL_fitOneModel(counts = x@counts, genotypes = x@genotypes, dispersion = slot(x, x@dispersion), model = "null", prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)


### Number of samples used for the analysis
n <- unlist(lapply(1:length(x@counts), function(g){
  sum(!is.na(x@counts[[g]][1, ]))
}))


results <- dmSQTL_test(fit_full = x@fit_full, fit_null = fit_null, test = test, n = n, verbose = verbose, return_list = TRUE, BPPARAM = BPPARAM)



max_nr_perm = 100
max_nr_sign_pval = 10
verbose = TRUE

dmSQTL_permutations_per_gene <- function(x, fit_null, results, max_nr_perm = 1e6, max_nr_sign_pval = 1e2, prop_mode, prop_tol, n, test, verbose = TRUE, BPPARAM){
  
  fit_full <- x@fit_full
  n <- ncol(x@counts)
  
  pval <- lapply(results, function(x){
    pval_tmp <- x$pvalue
    pval_tmp <- pval_tmp[!is.na(pval_tmp)]
    pval_tmp <- factor(pval_tmp)
    return(pval_tmp)
  })
  
  results_width <- unlist(lapply(pval, length))
  nas <- results_width == 0
  genes2permute <- which(!nas)
  
  sum_sign_pval <- vector("list", length(pval))
  sum_sign_pval[!nas] <- split(rep(0, sum(results_width)), factor(rep(1:length(results_width), times = results_width)))
  
  nr_perm_tot <- rep(0, length(x@counts))
  nr_perm_tot[nas] <- NA
  
  min_nr_sign_pval <- rep(0, length(x@counts))
  min_nr_sign_pval[nas] <- NA
  
  
  while(length(genes2permute) > 0){
    
    if(verbose)
    message(paste0(length(genes2permute), " genes left for permutation..\n"))
    
    permutation <- sample(n, n)
    
    ### Permute counts for all genes that need additional permutations
    counts <- x@counts[genes2permute, permutation]
    genotypes <- x@genotypes[genes2permute, ]
    dispersion <- slot(x, x@dispersion)
    if(is.list(dispersion))
      dispersion <- dispersion[genes2permute]
    
    fit_full_perm <- dmSQTL_fitOneModel(counts = counts, genotypes = genotypes, dispersion = dispersion, model = "full", prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
    
    fit_null_perm <- dmSQTL_fitOneModel(counts = counts, genotypes = genotypes, dispersion = dispersion, model = "null", prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
    
    results_perm <- dmSQTL_test(fit_full = fit_full_perm, fit_null = fit_null_perm, test = test, n = n, return_list = TRUE, verbose = verbose, BPPARAM = BPPARAM)
    
    
    ### Count how many pval_permuted is lower than pval from the model
    update_nr_sign_pval <- lapply(1:length(results_perm), function(i, results_perm, pval, genes2permute){
      # i = 1
      
      pval_perm <- results_perm[[i]]$pvalue
      pval_perm <- pval_perm[!is.na(pval_perm)]
      
      pval_perm_cut <- cut(pval_perm, c(-1, levels(pval[[genes2permute[i]]]), 2), right=FALSE)
      pval_perm_sum <- table(pval_perm_cut)
      
      pval_perm_cumsum <- cumsum(pval_perm_sum)[-length(pval_perm_sum)]
      names(pval_perm_cumsum) <- levels(pval[[genes2permute[i]]])
      nr_sign_pval <- pval_perm_cumsum[pval[[genes2permute[i]]]]
      
      return(nr_sign_pval)
      
    }, results_perm = results_perm, pval = pval, genes2permute = genes2permute)
    
    
    ### Update values in sum_sign_pval
    for(i in 1:length(update_nr_sign_pval)){
      sum_sign_pval[[genes2permute[i]]] <- sum_sign_pval[[genes2permute[i]]] + update_nr_sign_pval[[i]]
    }
    
    
    nr_perm <- unlist(lapply(results_perm, nrow))
    nr_perm_tot[genes2permute] <- nr_perm_tot[genes2permute] + nr_perm
    
    min_nr_sign_pval[genes2permute] <- unlist(lapply(sum_sign_pval[genes2permute], min))
    
    ### Update genes2permute
    genes2permute <- which(nr_perm_tot < max_nr_perm & min_nr_sign_pval < max_nr_sign_pval)
    
  }
  
  
  ### Calculate permutation adjusted p-values
  pval_adj <- lapply(1:length(results), function(i, results, sum_sign_pval, nr_perm_tot){
    
    pval_tmp <- results[[i]]$pvalue
    nas <- is.na(pval_tmp)
    
    if(sum(!nas) == 0)
      return(pval_tmp)
    
    pval_tmp[!nas] <- (sum_sign_pval[[i]] + 1) / (nr_perm_tot[i] + 1)
    
    return(pval_tmp)
    
    }, results = results, sum_sign_pval = sum_sign_pval, nr_perm_tot = nr_perm_tot)
  
  
  pval_out <- unlist(pval_adj)
  
  return(pval_out)

}


pval_permutations <- dmSQTL_permutations_per_gene(x, fit_null, results, max_nr_perm = 1e6, max_nr_sign_pval = 1e2, prop_mode, prop_tol, n, test, verbose = TRUE, BPPARAM)














