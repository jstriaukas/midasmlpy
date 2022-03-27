rm(list=ls())
#setwd("~/Documents/GitHub/midasmlpy/midasmlpy/data_files")

library(midasml)

set.seed(123)
x = matrix(rnorm(100 * 50), 100, 50)
beta = c(5,4,3,2,1,rep(0, times = 45))
y = x%*%beta + rnorm(100)
gamma = 0.5
gindex = sort(rep(1:10,times=5))
nlambda = 100L
method = "single"
nf = 0
lambda.factor = 1e-04
lambda = NULL
pf = rep(1, 50)
dfmax = 50 + 1
pmax = min(dfmax * 1.2, 50)
standardize = FALSE
intercept = TRUE
eps = 1e-08 
maxit = 1000000L
peps = 1e-08

fit <- sglfit(x, y, gamma = gamma, nlambda = nlambda, method = method, nf = nf, lambda.factor = lambda.factor, 
       lambda = lambda, pf = pf, gindex = gindex, 
       dfmax = dfmax, pmax = pmax, standardize = standardize, 
       intercept = intercept, eps = eps, maxit = maxit, peps = peps)


# compare:
fit$b0
fit$a0
fit$beta
fit$df
fit$dim
fit$lambda
fit$nf
fit$npasses
fit$jerr
fit$dimx




