

# ------------------------------------------------
# Test Hessian from nnet
# ------------------------------------------------

require(foreign)
require(nnet)
require(ggplot2)
require(reshape2)

ml <- read.dta("http://www.ats.ucla.edu/stat/data/hsbdemo.dta")

with(ml, table(ses, prog))

ml$prog2 <- relevel(ml$prog, ref = "academic")


test <- multinom(prog2 ~ ses + write, data = ml, Hess = TRUE)

test$Hessian

nnet:::multinomHess

test$fitted
coef(test)


object <- test
Z <- model.matrix(prog2 ~ ses + write, data = ml)


multinomHess <- function (object, Z = model.matrix(object)){
  
  probs <- object$fitted
  coefs <- coef(object)
  
  if (is.vector(coefs)) {
    coefs <- t(as.matrix(coefs))
    probs <- cbind(1 - probs, probs)
  }
  
  coefdim <- dim(coefs)
  p <- coefdim[2L]
  k <- coefdim[1L]
  ncoefs <- k * p
  kpees <- rep(p, k)
  n <- dim(Z)[1L]
  
  # Hessian 
  info <- matrix(0, ncoefs, ncoefs) 
  
  Names <- dimnames(coefs)
  
  if (is.null(Names[[1L]])){
    Names <- Names[[2L]]
  }else{ Names <- as.vector(outer(Names[[2L]], Names[[1L]], function(name2,
    name1) paste(name1, name2, sep = ":")))
  }
  
  dimnames(info) <- list(Names, Names)
  
  # ni
  row.totals <- object$weights
  
  
  for (i in seq_len(n)) {
    # i = 1
    
    Zi <- Z[i, ]
    Zi
    
    xbar <- Zi * rep(probs[i, -1, drop = FALSE], kpees)
    xbar
    
    # They use the first feature as a reference
    for (j in seq_len(k + 1)) {  # 
      # j = 3
      
      x <- matrix(0, p, k + 1L)
      x
      x[, j] <- Zi
      x
      x <- x[, -1, drop = FALSE]
      x
      
      x <- x - xbar
      x
      
      dim(x) <- c(1, ncoefs)
      x
      
      info_tmp <- (row.totals[i] * probs[i, j] * crossprod(x))
      info_tmp
      
      info <- info + info_tmp
      info
      
      
    }
    
  }
  
  info
  
  
}




# ------------------------------------------------
# Compare Hessian from nnet and mlogit and mnlogit
# ------------------------------------------------


library(mlogit)


## Cameron and Trivedi's Microeconometrics p.493 There are two
## alternative specific variables : price and catch one individual
## specific variable (income) and four fishing mode : beach, pier, boat,
## charter

data("Fishing", package = "mlogit")

Fish <- mlogit.data(Fishing, varying = c(2:9), shape = "wide", choice = "mode")

## a pure "multinomial model"

fit1 <- mlogit(mode ~ 0 | income, data = Fish)

summary(fit1)

fit1$hessian # Negative values on diagonal




library(nnet)

## which can also be estimated using multinom (package nnet)

fit2 <- multinom(mode ~ income, data = Fishing, Hess = TRUE)

summary(fit2)

fit2$Hessian # Positive values on diagonal




library(mnlogit)

fit3 <- mnlogit(formula(mode ~ 1 | income), Fish)

summary(fit3)

fit3$hessian # Positive values on diagonal



# ------------------------------------------------
# Compare Hessian from nnet and mlogit and mnlogit
# ------------------------------------------------

library(mlogit)
library(nnet)
library(mnlogit)


ml <- read.dta("http://www.ats.ucla.edu/stat/data/hsbdemo.dta")

ml$prog2 <- relevel(ml$prog, ref = "academic")


mlwide <- mlogit.data(ml, shape = "wide", choice = "prog2")




fit1 <- mlogit(prog2 ~ 0 | ses + write, data = mlwide)

summary(fit1)

coef(fit1)

fit1$hessian # Negative values on diagonal




## which can also be estimated using multinom (package nnet)

fit2 <- multinom(prog2 ~ ses + write, data = ml, Hess = TRUE)

summary(fit2)

fit2$Hessian # Positive values on diagonal




fit3 <- mnlogit(formula(prog2 ~ 1 | ses + write), mlwide)

summary(fit3)

fit3$hessian # Positive values on diagonal




### Test my function 

prop <- fit2$fitted

x <- model.matrix(fit2)

m <- fit2$weights[,1]

p <- ncol(x)
q <- ncol(prop)


H <- matrix(0, p*(q-1), p*(q-1))

rownames(H) <- colnames(H) <- paste0(rep(colnames(prop)[-1], each = p), ":", rep(colnames(x), q-1))

jp_index <- matrix(1:p, nrow = q-1, ncol = p, byrow = TRUE) + p * 0:(q-2)

for(jp in 2:q){
  
  for(jpp in 2:q){
    # jp = 1; jpp = 1
    
    W <- m * prop[, jp] * (prop[, jpp] - as.numeric(jp == jpp))
    
    h <- t(x) %*% (W * x)
    
    H[jp_index[jp-1, ], jp_index[jpp-1, ]] <- h
    
  }
  
}

H

round(fit2$Hessian, 6) == - round(H, 6)


# Using the symetry of H

H <- matrix(0, p*(q-1), p*(q-1))

rownames(H) <- colnames(H) <- paste0(rep(colnames(prop)[-1], each = p), ":", rep(colnames(x), q-1))

jp_index <- matrix(1:p, nrow = q-1, ncol = p, byrow = TRUE) + p * 0:(q-2)

for(jp in 2:q){
  
  for(jpp in 2:jp){
    # jp = 2; jpp = 2
    
    W <- m * prop[, jp] * (prop[, jpp] - as.numeric(jp == jpp))
    
    h <- t(x) %*% (W * x)
    
    H[jp_index[jp-1, ], jp_index[jpp-1, ]] <- h
    
    if(!jp == jpp){
      H[jp_index[jpp-1, ], jp_index[jp-1, ]] <- t(h)
    }
    
  }
  
}

H

round(fit2$Hessian, 6) == - round(H, 6)





























