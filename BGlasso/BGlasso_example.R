rm(list=ls())
library('BayesianGLasso') # Must install this package
source('Sims/ROC_calcs.R')
library(huge)

# Example of how to use the Bayesian graphical lasso by Wang et al. (2012) 
n=200
p=100
set.seed(123)
data.sf = huge::huge.generator(n=n, d=p,graph = 'scale-free') 
theta.true = data.sf$omega # The precision matrix
theta.true[which(theta.true<10e-5,arr.ind=T)]=0
x.sf = mvtnorm::rmvnorm(n, sigma=solve(theta.true)) # Sample data
x.sf.scaled = scale(x.sf) # Scale columns/variables.
data.sf$sparsity # True sparsity: 0.02

# Using the standard prior choices as proposed by Wang et al. (2012) in their package
res.bglasso = blockGLasso(x.sf.scaled,iterations = 1000, burnIn = 100, lambdaPriora = 1, lambdaPriorb = 1/10)
theta.est.bglasso = apply(simplify2array(res.bglasso$Omegas), c(1,2), mean) # Posterior mean of MCMC samples

# Threshold the same way as GHS, varying the threshold for the credible interval to use. 

