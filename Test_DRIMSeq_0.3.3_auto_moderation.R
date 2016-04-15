# R32dev2
# Created 14 Apr 2016



##############################################################################

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(args)

# rwd='/home/gosia/multinomial_project/package_devel/Test_DRIMSeq_0.3.3_auto_moderation'
# count_method='kallisto'


##############################################################################


dir.create(rwd, recursive = TRUE)

setwd(rwd)


##############################################################################
# Use Pasilla data
##############################################################################

library(BiocParallel)
library(DRIMSeq)
library(limma)
library(reshape2)
library(ggplot2)

##############################################################################

rwd='/home/Shared/data/seq/brooks_pasilla/'
workers=15

model='model_full'

method_out <- "drimseq_0_3_3"

if(workers > 1){
  BPPARAM <- MulticoreParam(workers = workers)
}else{
  BPPARAM <- SerialParam()
}

out_dir <- "./"

##########################################################################
# load metadata
##########################################################################

metadata <- read.table(paste0(rwd, "3_metadata/metadata.xls"), stringsAsFactors = FALSE, sep="\t", header=TRUE) 

metadata_org <- metadata


##########################################################################
# DRIMSeq analysis - prepare counts
##########################################################################


count_dir <- paste0(rwd, "2_counts/", count_method, "/")


### load counts
counts_list <- lapply(1:length(metadata_org$sampleName), function(i){
  # i = 1
  cts <- read.table(paste0(count_dir, metadata_org$sampleName[i], ".txt"), header = FALSE, as.is = TRUE)
  colnames(cts) <- c("group_id", metadata_org$sampleName[i])  
  return(cts)
})

counts <- Reduce(function(...) merge(..., by = "group_id", all=TRUE, sort = FALSE), counts_list)
counts <- counts[!grepl(pattern = "_", counts$group_id),]


### Prepare data
group_split <- strsplit2(counts[,1], ":")
counts <- counts[, -1]
### order the samples like in metadata!!!
counts <- counts[, metadata_org$sampleName]


metadata <- metadata_org
all(colnames(counts) == metadata$sampleName)

d <- dmDSdata(counts = counts, gene_id = group_split[, 1], feature_id = group_split[, 2], sample_id = metadata$sampleName, group = metadata$condition)

### Filtering
table(samples(d)$group)
d <- dmFilter(d, min_samps_gene_expr = 0, min_samps_feature_expr = 0, min_samps_feature_prop = 0, min_gene_expr = 0, min_feature_expr = 0, min_feature_prop = 0, max_features = Inf)



##########################################################################
# DRIMSeq analysis - estimate dispersion
##########################################################################


### DRIMSeq pipelines : genewise_dispersion

common_disp <- as.numeric(read.table(paste0(rwd, method_out, "/",  model, "/", count_method, "/", "common_dispersion.txt")))
common_disp



# genewise dispersion
d <- dmDispersion(d, common_dispersion = FALSE, genewise_dispersion = TRUE, disp_mode = "grid", disp_init = common_disp, disp_moderation = "none", disp_prior_df = 1, verbose = TRUE, BPPARAM = BPPARAM)


common_dispersion(d) <- common_disp

plotDispersion(d, out_dir = paste0(count_method, "_moderation_none_"))


######################################
# Implement automatic moderation
######################################

library(devtools)
load_all("/home/gosia/R/package_devel/DRIMSeq")


x <- d

mean_expression = TRUE
common_dispersion = FALSE
genewise_dispersion = TRUE
disp_adjust = TRUE
disp_mode = "grid"
disp_interval = c(0, 1e+05)
disp_tol = 1e-08
disp_init = common_disp
disp_init_weirMoM = TRUE
disp_grid_length = 21
disp_grid_range = c(-10, 10)
disp_moderation = "trended"
disp_prior_df = 1
disp_span = 0.3
prop_mode = "constrOptimG"
prop_tol = 1e-12
verbose = TRUE

###### Inside the dmDispersion function

mean_expression <- dm_estimateMeanExpression(counts = x@counts, verbose = verbose, BPPARAM = BPPARAM)


# genewise_dispersion <- dmDS_estimateTagwiseDispersion(counts = x@counts, samples = x@samples, mean_expression = mean_expression, disp_adjust = disp_adjust, disp_mode = disp_mode, disp_interval = disp_interval, disp_tol = disp_tol, disp_init = disp_init, disp_init_weirMoM = disp_init_weirMoM, disp_grid_length = disp_grid_length, disp_grid_range = disp_grid_range, disp_moderation = disp_moderation, disp_prior_df = disp_prior_df, disp_span = disp_span, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)


counts = x@counts
samples = x@samples



###### Inside the dmDS_estimateTagwiseDispersion function


inds <- 1:length(counts)

group <- samples$group
ngroups <- nlevels(group)
lgroups <- levels(group)
nlibs <- length(group)

igroups <- lapply(lgroups, function(gr){which(group == gr)})
names(igroups) <- lgroups



###### Switch to grid

splinePts <- seq(from = disp_grid_range[1], to = disp_grid_range[2], length = disp_grid_length)



disp_grid_length <- length(splinePts)


splineDisp <- disp_init * 2^splinePts

### calculate the likelihood for each gene at the spline dispersion points
seq_disp_grid_length <- seq(disp_grid_length)

loglikL <- BiocParallel::bplapply(inds, function(g){
  # g = 1237
  # print(g)
  
  ll <- numeric(disp_grid_length)
  
  for(i in seq_disp_grid_length){
    # i = 1
    
    out <- dm_profileLikTagwise(gamma0 = splineDisp[i], y = counts[[g]], ngroups = ngroups, lgroups = lgroups, igroups = igroups, disp_adjust = disp_adjust, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose)
    
    if(is.na(out)){
      ll <- rep(NA, disp_grid_length)
      break
    }
    
    ll[i] <- out
    
  }
  
  return(ll)
  
}, BPPARAM = BPPARAM)



loglik <- do.call(rbind, loglikL)

not_nas <- complete.cases(loglik)        

loglik <- loglik[not_nas, , drop = FALSE]


loglik_orig <- loglik





### Check where the grid is maximized 
grid_max <- apply(loglik, 1, which.max)
lik_max <- apply(loglik, 1, max)

range <- t(apply(loglik, 1, range))
lik_span <- range[, 2] - range[, 1]



### In the calculation of moderation, do not take into account genes that have dispersion on the top and bottom boundry of the grid (skipp 4 last grid points and 1 first grid point)
not_boundry <- grid_max < (disp_grid_length - 3) & grid_max > 1


boundry_last <- grid_max >= disp_grid_length - 3
table(boundry_last)



############################################################################
### Auto moderation for adjustment to common dispersion

prefix_mod <- paste0(count_method, "_common_")

## Calculate common likelihood
if(sum(not_boundry) == length(not_boundry)){
  moderation <- colMeans(loglik)
}else{
  moderation <- colMeans(loglik[not_boundry, , drop = FALSE])
}

moderation_span <- max(moderation) - min(moderation)




# Plot common likelihood and boundry likelihoods

df_liks_boundry <- data.frame(t(loglik[boundry_last, , drop = FALSE]))

df_liks_boundry <- data.frame(splinePts = splinePts, moderation = moderation, df_liks_boundry)

df_liks_boundrym <- melt(df_liks_boundry, id.vars = "splinePts", variable.name = "gene", value.name = "loglik")
df_liks_boundrym$moderation <- grepl("moderation", df_liks_boundrym$gene)

ggp <- ggplot(df_liks_boundrym, aes(x = splinePts, y = -log10(-loglik), group = gene, color = moderation, size = moderation)) +
  geom_line() +
  scale_size_manual(values = c(0.5, 2)) +
  coord_cartesian(ylim = range(-log10(-lik_max)))

pdf(paste0(prefix_mod, "liks_boundry_and_moderation.pdf"))
print(ggp)
dev.off()



# Plot lik max versus mean expression 

df_lik_max_mean <- data.frame(loglik_max = lik_max, mean_expression = mean_expression[not_nas], boundry = boundry_last)

ggp <- ggplot(df_lik_max_mean, aes(x = log10(mean_expression), y = -log10(-loglik_max), color = boundry)) +
  geom_point(alpha = 1, size = 0.5)

pdf(paste0(prefix_mod, "liks_max_versus_mean_expression.pdf"))
print(ggp)
dev.off()



# Plot boxplot of likelihood span for boundry and not boundry genes

df_lik_span <- data.frame(loglik_span = lik_span, boundry = boundry_last)

ggp <- ggplot(df_lik_span, aes(x = boundry, y = log10(loglik_span), fill = boundry)) +
  geom_boxplot(outlier.size = 0) +
  geom_point(position = position_jitter(width = 0.5), alpha = 0.2) +
  geom_abline(intercept = log10(moderation_span), slope = 0, linetype = 2) 

pdf(paste0(prefix_mod, "liks_span_log.pdf"))
print(ggp)
dev.off()


ggp <- ggplot(df_lik_span, aes(x = boundry, y = loglik_span, fill = boundry)) +
  geom_boxplot(outlier.size = 0) +
  geom_point(position = position_jitter(width = 0.5), alpha = 0.2) +
  geom_abline(intercept = moderation_span, slope = 0, linetype = 2) +
  coord_cartesian(ylim = c(0, 2*moderation_span))

pdf(paste0(prefix_mod, "liks_span.pdf"))
print(ggp)
dev.off()


# Plot lik span versus mean expression 

df_lik_span_mean <- data.frame(loglik_span = lik_span, mean_expression = mean_expression[not_nas], boundry = boundry_last)

ggp <- ggplot(df_lik_span_mean, aes(x = log10(mean_expression), y = log10(loglik_span), color = boundry)) +
  geom_point(alpha = 1, size = 0.5) +
  geom_abline(intercept = log10(moderation_span), slope = 0, linetype = 2)

pdf(paste0(prefix_mod, "liks_span_versus_mean_expression_log.pdf"))
print(ggp)
dev.off()


ggp <- ggplot(df_lik_span_mean, aes(x = log10(mean_expression), y = loglik_span, color = boundry)) +
  geom_point(alpha = 1, size = 0.5) +
  geom_abline(intercept = moderation_span, slope = 0, linetype = 2) +
  coord_cartesian(ylim = c(0, 2*moderation_span))

pdf(paste0(prefix_mod, "liks_span_versus_mean_expression.pdf"))
print(ggp)
dev.off()



### Calculate the ratio between moderation lik span and lik span of boundry genes

loglik_span_boundry <- df_lik_span[boundry_last, "loglik_span"]

span_ratio <- moderation_span / loglik_span_boundry

priorN <- quantile(1/span_ratio, 0.5)

message(paste0("! Using ", round(priorN, 2), " as a shrinkage factor !"))

loglik <- sweep(loglik, 2, priorN * moderation, FUN = "+")

grid_max_new <- apply(loglik, 1, which.max)

grid_max_changed <- grid_max_new != grid_max


### Plot priorN versus mean expression

df_priorN <- data.frame(priorN = 1/span_ratio, mean_expression = df_lik_span_mean[boundry_last, "mean_expression"])

ggp <- ggplot(df_priorN, aes(x = log10(mean_expression), y = priorN)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = priorN, slope = 0, linetype = 2) +
  coord_cartesian(xlim = log10(range(df_lik_span_mean$mean_expression))) +
  ggtitle(paste0("priorN = ", round(priorN, 4)))

pdf(paste0(prefix_mod, "priorn_versus_mean_expression.pdf"))
print(ggp)
dev.off()


pdf(paste0(prefix_mod, "priorn_versus_mean_expression_smoothscatter.pdf"))
smoothScatter(x = log10(df_priorN$mean_expression), y = df_priorN$priorN, xlim = log10(range(df_lik_span_mean$mean_expression)), main = paste0("priorN = ", round(priorN, 4)))
abline(a = priorN, b = 0, lty = 2, lwd = 2)
dev.off()



### Plot lik span versus mean expression + mark genes that changed their lik maximum

df_lik_span_mean <- data.frame(loglik_span = lik_span, mean_expression = mean_expression[not_nas], boundry = boundry_last, grid_max_changed = grid_max_changed)

ggp <- ggplot(df_lik_span_mean, aes(x = log10(mean_expression), y = log10(loglik_span), color = boundry)) +
  geom_point(alpha = 1, size = 0.5) +
  geom_abline(intercept = log10(moderation_span), slope = 0, linetype = 2) +
  facet_wrap(~grid_max_changed)

pdf(paste0(prefix_mod, "liks_span_versus_mean_expression_grid_max_changed_log.pdf"), 14, 7)
print(ggp)
dev.off()


ggp <- ggplot(df_lik_span_mean, aes(x = log10(mean_expression), y = loglik_span, color = boundry)) +
  geom_point(alpha = 1, size = 0.5) +
  geom_abline(intercept = moderation_span, slope = 0, linetype = 2) +
  coord_cartesian(ylim = c(0, 2*moderation_span)) +
  facet_wrap(~grid_max_changed)

pdf(paste0(prefix_mod, "liks_span_versus_mean_expression_grid_max_changed.pdf"), 14, 7)
print(ggp)
dev.off()







############################################################################
### Auto moderation for adjustment to trended dispersion

loglik <- loglik_orig


prefix_mod <- paste0(count_method, "_trended_")

mean_expression <- mean_expression[not_nas]

if(sum(not_boundry) == length(not_boundry)){
  
  o <- order(mean_expression)
  oo <- order(o)
  width <- floor(disp_span * nrow(loglik))
  
  moderation <- edgeR::movingAverageByCol(loglik[o,], width = width)[oo,]
  
}else{
  
  ### Use non boundry genes for calculating the moderation
  mean_expression_not_boundry <- mean_expression[not_boundry]
  loglik_not_boundry <- loglik[not_boundry, , drop = FALSE]
  
  o <- order(mean_expression_not_boundry)
  oo <- order(o)
  
  width <- floor(disp_span * nrow(loglik_not_boundry))
  
  moderation_not_boundry <- edgeR::movingAverageByCol(loglik_not_boundry[o, , drop = FALSE], width = width)[oo, , drop = FALSE]
  
  ### Fill in moderation values for the boundy genes
  moderation <- matrix(NA, nrow = nrow(loglik), ncol = ncol(loglik))
  
  moderation[not_boundry, ] <- moderation_not_boundry
  
  o <- order(mean_expression)
  oo <- order(o)
  
  moderation <- moderation[o, , drop = FALSE]
  not_boundry <- not_boundry[o]
  
  ### Last value in not_boundry must be TRUE
  if(not_boundry[length(not_boundry)] == FALSE){
    
    last_true <- max(which(not_boundry))
    moderation[length(not_boundry), ] <- moderation[last_true, ]
    
    not_boundry[length(not_boundry)] <- TRUE
    
  }
  
  not_boundry_diff <- diff(not_boundry, lag = 1)
  
  not_boundry_cumsum <- cumsum(not_boundry)
  
  ### Values used for filling in the boundry NAs - switch from FALSE to TRUE
  replacement_indx <- which(not_boundry_diff == 1) + 1
  
  replaced_indx <- which(!not_boundry)
  
  replaced_freq <- as.numeric(table(not_boundry_cumsum[replaced_indx]))
  
  moderation_boundry  <- moderation[rep(replacement_indx, times = replaced_freq), , drop = FALSE]
  
  moderation[!not_boundry, ] <- moderation_boundry
  
  moderation <- moderation[oo, , drop = FALSE]
  
}


moderation_span <- apply(moderation, 1, function(x){max(x) - min(x)})


# Plot lik span versus mean expression 

df_lik_span_mean <- data.frame(loglik_span = lik_span, mean_expression = mean_expression, boundry = boundry_last)

df_moderation_span <- data.frame(moderation_span = moderation_span, mean_expression = mean_expression, boundry = boundry_last)

ggp <- ggplot(df_lik_span_mean, aes(x = log10(mean_expression), y = log10(loglik_span), color = boundry)) +
  geom_point(alpha = 1, size = 0.5) +
  geom_point(data = df_moderation_span, aes(x = log10(mean_expression), y = log10(moderation_span)), color = "darkgrey")

pdf(paste0(prefix_mod, "liks_span_versus_mean_expression_log.pdf"))
print(ggp)
dev.off()


ggp <- ggplot(df_lik_span_mean, aes(x = log10(mean_expression), y = loglik_span, color = boundry)) +
  geom_point(alpha = 1, size = 0.5) +
  geom_point(data = df_moderation_span, aes(x = log10(mean_expression), y = moderation_span), color = "darkgrey") +
  coord_cartesian(ylim = c(0, 2*max(moderation_span)))

pdf(paste0(prefix_mod, "liks_span_versus_mean_expression.pdf"))
print(ggp)
dev.off()



### Calculate the ratio between moderation lik span and lik span of boundry genes

loglik_span_boundry <- df_lik_span_mean[boundry_last, "loglik_span"]
moderation_span_boundry <- df_moderation_span[boundry_last, "moderation_span"]


span_ratio <- moderation_span_boundry / loglik_span_boundry

priorN <- quantile(1/span_ratio, 0.5)

message(paste0("! Using ", round(priorN, 2), " as a shrinkage factor !"))

loglik <- loglik + priorN * moderation 

grid_max_new <- apply(loglik, 1, which.max)

grid_max_changed <- grid_max_new != grid_max




### Plot priorN versus mean expression

df_priorN <- data.frame(priorN = 1/span_ratio, mean_expression = df_lik_span_mean[boundry_last, "mean_expression"])

ggp <- ggplot(df_priorN, aes(x = log10(mean_expression), y = priorN)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = priorN, slope = 0, linetype = 2) +
  coord_cartesian(xlim = log10(range(df_lik_span_mean$mean_expression))) +
  ggtitle(paste0("priorN = ", round(priorN, 4)))

pdf(paste0(prefix_mod, "priorn_versus_mean_expression.pdf"))
print(ggp)
dev.off()


pdf(paste0(prefix_mod, "priorn_versus_mean_expression_smoothscatter.pdf"))
smoothScatter(x = log10(df_priorN$mean_expression), y = df_priorN$priorN, xlim = log10(range(df_lik_span_mean$mean_expression)), main = paste0("priorN = ", round(priorN, 4)))
abline(a = priorN, b = 0, lty = 2, lwd = 2)
dev.off()


### Plot lik span versus mean expression + mark genes that changed their lik maximum

df_lik_span_mean <- data.frame(loglik_span = lik_span, mean_expression = mean_expression, boundry = boundry_last, grid_max_changed = grid_max_changed)

ggp <- ggplot(df_lik_span_mean, aes(x = log10(mean_expression), y = log10(loglik_span), color = boundry)) +
  geom_point(alpha = 1, size = 0.5) +
  geom_point(data = df_moderation_span, aes(x = log10(mean_expression), y = log10(moderation_span)), color = "darkgrey") +
  facet_wrap(~grid_max_changed)

pdf(paste0(prefix_mod, "liks_span_versus_mean_expression_grid_max_changed_log.pdf"), 14, 7)
print(ggp)
dev.off()


ggp <- ggplot(df_lik_span_mean, aes(x = log10(mean_expression), y = loglik_span, color = boundry)) +
  geom_point(alpha = 1, size = 0.5) +
  geom_point(data = df_moderation_span, aes(x = log10(mean_expression), y = moderation_span), color = "darkgrey") +
  coord_cartesian(ylim = c(0, 2*max(moderation_span))) +
  facet_wrap(~grid_max_changed)

pdf(paste0(prefix_mod, "liks_span_versus_mean_expression_grid_max_changed.pdf"), 14, 7)
print(ggp)
dev.off()




##########################################################################
# DRIMSeq analysis - new dispersion
##########################################################################


# moderation to common dispersion
d <- dmDispersion(d, common_dispersion = FALSE, genewise_dispersion = TRUE, disp_mode = "grid", disp_init = common_disp, disp_moderation = "common", disp_prior_df = 0.01, verbose = TRUE, BPPARAM = BPPARAM)


common_dispersion(d) <- common_disp

plotDispersion(d, out_dir = paste0(count_method, "_moderation_common_new_"))


# moderation to trended dispersion
d <- dmDispersion(d, common_dispersion = FALSE, genewise_dispersion = TRUE, disp_mode = "grid", disp_init = common_disp, disp_moderation = "trended", disp_prior_df = 0.01, verbose = TRUE, BPPARAM = BPPARAM)


common_dispersion(d) <- common_disp

plotDispersion(d, out_dir = paste0(count_method, "_moderation_trended_new_"))


















