###############################################
# Example

#Load packages
library(MASS)
library(Matrix)
library(Rcpp)
library(RcppArmadillo)

p = 30                  
n = 150                     # sample size for each model.

rho = 0.0                   # misspecification ratio
sparsity <- 0.08            # sparsity level
ua <- 0.4
ub <- 0.6

K = 4       # number of models/networks 


Omega <- Gen_Precison_Mats(K, p, rho, ua, ub, sparsity)
Sigma <- vector("list", K)
x <- vector("list", K)

for (k in 1:K){
  Sigma[[k]] <- solve(Omega[[k]])
  x[[k]] <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma[[k]])
  x[[k]] <- scale(x[[k]], center = T, scale = T)
}

matlist <- c()

for(k in 1:K){
  matlist <- c(matlist, combn(0:(K-1), k, simplify = FALSE))
}

fit <- BJNS(x, matlist)





