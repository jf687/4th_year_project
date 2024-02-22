rm(list=ls())
source('libraries.R')
source('Sims/ROC_calcs.R')


# Example of how to use the Bayesian graphical lasso by Wang et al. (2012)
bglasso_ROC_sim <- function(n=200,p=100,plot = F){
  set.seed(123)
  data.sf = huge::huge.generator(n=n, d=p,graph = 'scale-free') 
  theta.true = data.sf$omega # The precision matrix
  theta.true[which(theta.true<10e-5,arr.ind=T)]=0
  y.sf = mvtnorm::rmvnorm(n, sigma=solve(theta.true)) # Sample data
  y.sf.scaled = scale(y.sf) # Scale columns/variables.
  data.sf$sparsity # True sparsity: 0.02
  
  # Using the standard prior choices as proposed by Wang et al. (2012) in their package
  res.bglasso = blockGLasso(y.sf.scaled,iterations = 1000, burnIn = 100, lambdaPriora = 1, lambdaPriorb = 1/10)
  theta.est.bglasso = apply(simplify2array(res.bglasso$Omegas), c(1,2), mean) # Posterior mean of MCMC samples
  
  # Threshold the same way as GHS, varying the threshold for the credible interval to use. 
  
  
  bool_up <- upper.tri(theta.est.bglasso)
  
  labels <- theta.true[bool_up]>0
  predictions <- theta.est.bglasso[bool_up]
  
  
  pred <- ROCR::prediction(predictions, labels)
  perf <- ROCR::performance(pred, measure = "tpr", x.measure = "fpr")
  
  auc_value <- ROCR::performance(pred, "auc")@y.values[[1]]
  prec_value <- ROCR::performance(pred, "prec")@y.values[[1]]
  
  if(plot){
    plot(perf, main = '', xlim = c(0, 1), ylim = c(0, 1), lwd = 2, colorize = TRUE)
  }
  
  return(list(n = n, p = p, auc_value = auc_value, prec_value = prec_value))
}

