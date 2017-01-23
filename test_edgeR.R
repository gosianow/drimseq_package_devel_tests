nlibs <- 9
ngenes <- 100
dispersion.true <- 0.1

# Make first gene respond to covariate x
x <- 0:(nlibs-1)
design <- model.matrix(~x)
beta.true <- cbind(Beta1=2,Beta2=c(2,rep(0,ngenes-1)))
mu.true <- 2^(beta.true %*% t(design))

# Generate count data
y <- rnbinom(ngenes*nlibs,mu=mu.true,size=1/dispersion.true)
y <- matrix(y,ngenes,nlibs)
colnames(y) <- paste0("x", 1:nlibs)
rownames(y) <- paste("gene",1:ngenes,sep=".")

d <- DGEList(y)

# Normalize
d <- calcNormFactors(d)


fit <- glmFit(d, design, dispersion=dispersion.true)

design0 <- design[, -ncol(design), drop = FALSE]

fit0 <- glmFit(d, design0, dispersion=dispersion.true)

designAsFactor(design)
designAsFactor(design0)



new_x <- factor(c(rep("A", nlibs/3), rep("B", nlibs/3), rep("C", nlibs/3)))
new_design <- model.matrix(~new_x)

# Fit the NB GLMs
fit <- glmFit(d, new_design, dispersion=dispersion.true)


new_design0 <- new_design[, -ncol(new_design), drop = FALSE]

fit0 <- glmFit(d, new_design0, dispersion=dispersion.true)

designAsFactor(new_design)
designAsFactor(new_design0)




new_x <- factor(c(rep("A", nlibs/3), rep("B", nlibs/3), rep("C", nlibs/3)))
new_design <- model.matrix(~ 0 + new_x)

# Fit the NB GLMs
fit <- glmFit(d, new_design, dispersion=dispersion.true)
fit$method

new_design0 <- new_design[, -ncol(new_design), drop = FALSE]

fit0 <- glmFit(d, new_design0, dispersion=dispersion.true)
fit0$method

designAsFactor(new_design)
designAsFactor(new_design0)






# Likelihood ratio tests for trend
results <- glmLRT(fit, coef=2)
topTags(results)

# Estimate the dispersion (may be unreliable with so few genes)
d <- estimateGLMCommonDisp(d, design, verbose=TRUE)
		 
		 
		 
		 