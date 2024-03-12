rm(list=ls())
source('GHS/GHS.R')
source('Sims/ROC_calcs.R')
library(huge)


# Example of how to use GHS (graphical horseshoe)
GHS_ROC_sim <- function(n=50, p = 100){
  set.seed(123)
  data.sf = huge::huge.generator(n=n, d=p,graph = 'scale-free') 
  theta.true = data.sf$omega # The precision matrix
  theta.true[which(theta.true<10e-5,arr.ind=T)]=0
  y.sf = mvtnorm::rmvnorm(n, sigma=solve(theta.true)) # Sample data
  y.sf.scaled = scale(y.sf) # Scale columns/variables.
  data.sf$sparsity # True sparsity: 0.02
  
  ghs.res = GHS(t(y.sf.scaled)%*%y.sf.scaled,n,burnin=100,nmc=1000)
  theta.est.ghs = apply(ghs.res$thetas.sampled, c(1,2), mean) # Posterior mean of MCMC samples

  # To threshold, Li et al. (2017) use the symmetric central 50% posterior credible intervals for variable selection. 
  # That is, if the 50% posterior credible interval of an off-diagonal element of does not contain zero, that element is considered a discovery, and vice versa. 
  theta.ghs.credible.lower = apply(ghs.res$thetas.sampled, c(1,2), FUN = function(s) quantile(s, 0.25))
  theta.ghs.credible.upper = apply(ghs.res$thetas.sampled, c(1,2), FUN = function(s) quantile(s, 0.75))
  theta.ghs.selected = theta.est.ghs
  theta.ghs.selected[(theta.ghs.credible.lower < 0 & theta.ghs.credible.upper > 0)] = 0
  sparsity = sparsity(theta.ghs.selected)
  precision = precision(theta.true!=0, theta.ghs.selected!=0)
  recall = recall(theta.true!=0, theta.ghs.selected!=0)

  # To vary the threshold for ROC curves etc, you can increase/decrease the width of the credible intervals used to determine nonzero effects. 
  # I.e., replace 0.25 and 0.75 by x and 1.0-x for any 0<x<0.5.
  quant.thresh = seq(0,0.5,length.out=50)
  theta.path.ghs = list()
  TPR.ghs = rep(0,length(quant.thresh))
  FPR.ghs = rep(0,length(quant.thresh))
  
  pb <- progress_bar$new(format = "[:bar] :percent ETA: :eta", total = length(quant.thresh))
  
  for(x in 1:length(quant.thresh)){
    
    pb$tick()
    theta.ghs.credible.lower = apply(ghs.res$thetas.sampled, c(1,2), FUN = function(s) quantile(s, quant.thresh[x]))
    theta.ghs.credible.upper = apply(ghs.res$thetas.sampled, c(1,2), FUN = function(s) quantile(s, 1-quant.thresh[x]))
    theta.ghs.selected = theta.est.ghs
    theta.ghs.selected[(theta.ghs.credible.lower < 0 & theta.ghs.credible.upper > 0)] = 0
    diag(theta.ghs.selected) = 1
    theta.path.ghs[[x]] = theta.ghs.selected
    TPR.ghs[x] = TPR(theta.true!=0, theta.ghs.selected!=0)
    FPR.ghs[x] = FPR(theta.true!=0, theta.ghs.selected!=0)
  }
  
  bool_up <- upper.tri(theta.true)
  
  labels <- theta.true[bool_up]>0
  predictions <- theta.est.ghs[bool_up]
  
  
  pred <- ROCR::prediction(predictions, labels)
  perf <- ROCR::performance(pred, measure = "tpr", x.measure = "fpr")
  
  thresh <- return_closest_threshold(pred@cutoffs[[1]])

  auc_value <- ROCR::performance(pred, "auc")@y.values[[1]]
  prec_value <- ROCR::performance(pred, "prec")@y.values[[1]][thresh]
  rec_value <- ROCR::performance(pred, "rec")@y.values[[1]][thresh]
  
  
  plot(perf, main = '', xlim = c(0, 1), ylim = c(0, 1), lwd = 2, colorize = TRUE)
  
  # Plot ROC curve manually
  
  plot(FPR.ghs,TPR.ghs,type='l', col='red')
  
  x = c(0,FPR.ghs)
  y = c(0,TPR.ghs)
  AUC=flux::auc(x, y)

  return(list(n = n, p = p, AUC = AUC,auc_value = auc_value, precision = precision, recall = recall, sparsity = sparsity))

}




