
source('GHS/GHS.R')
source('Sims/ROC_calcs.R')
library(huge)


# Example of how to use GHS (graphical horseshoe)
GHS_ROC_sim <- function(data.sf, y.sf, x = 0.25){
  if(x >0.5){
    x <- 0.25
  }
  data.sf = huge::huge.generator(n=n, d=p,graph = 'scale-free') 
  theta.true = data.sf$omega # The precision matrix
  theta.true[which(theta.true<10e-5,arr.ind=T)]=0
  
  n <- dim(data.sf$data)[[1]]
  p <- dim(data.sf$data)[[2]]
  
  #data.sf$sparsity # True sparsity: 0.02
  
  ghs.res = GHS(t(y.sf)%*%y.sf,n,burnin=100,nmc=1000)
  theta.est.ghs = apply(ghs.res$thetas.sampled, c(1,2), mean) # Posterior mean of MCMC samples
  
  lambda = mean(ghs.res$lambdas.sampled)

  # To threshold, Li et al. (2017) use the symmetric central 50% posterior credible intervals for variable selection. 
  # That is, if the 50% posterior credible interval of an off-diagonal element of does not contain zero, that element is considered a discovery, and vice versa. 
  theta.ghs.credible.lower = apply(ghs.res$thetas.sampled, c(1,2), FUN = function(s) quantile(s, 0.25))
  theta.ghs.credible.upper = apply(ghs.res$thetas.sampled, c(1,2), FUN = function(s) quantile(s, 0.75))
  theta.ghs.selected = theta.est.ghs
  theta.ghs.selected[(theta.ghs.credible.lower < 0 & theta.ghs.credible.upper > 0)] = 0
  
  
  diag(theta.ghs.selected) <- 0
  diag(theta.true) <- 0
  
  browser()
  
  sparsity = sparsity(theta.ghs.selected)
  precision = precision(theta.true!= 0, theta.ghs.selected!= 0)
  recall = recall(theta.true!= 0, theta.ghs.selected!= 0)
  
  # To vary the threshold for ROC curves etc, you can increase/decrease the width of the credible intervals used to determine nonzero effects. 
  # I.e., replace 0.25 and 0.75 by x and 1.0-x for any 0<x<0.5.
  quant.thresh = seq(x,1-x,length.out=50)
  theta.path.ghs = list()
  TPR.ghs = rep(0,length(quant.thresh))
  FPR.ghs = rep(0,length(quant.thresh))
  
  pb <- progress_bar$new(format = "[:bar] :percent ETA: :eta", total = length(quant.thresh))
  
  for(l in 1:length(quant.thresh)){
    
    pb$tick()
    theta.ghs.credible.lower = apply(ghs.res$thetas.sampled, c(1,2), FUN = function(s) quantile(s, quant.thresh[l]))
    theta.ghs.credible.upper = apply(ghs.res$thetas.sampled, c(1,2), FUN = function(s) quantile(s, 1-quant.thresh[l]))
    theta.ghs.selected = theta.est.ghs
    theta.ghs.selected[(theta.ghs.credible.lower < 0 & theta.ghs.credible.upper > 0)] = 0
    diag(theta.ghs.selected) = 1
    theta.path.ghs[[l]] = theta.ghs.selected
    TPR.ghs[l] = TPR(theta.true!=0, theta.ghs.selected!=0)
    FPR.ghs[l] = FPR(theta.true!=0, theta.ghs.selected!=0)
  }
  
  x = c(0,FPR.ghs)
  y = c(0,TPR.ghs)
  AUC=flux::auc(c(0, FPR.ghs), c(0,TPR.ghs))

  return(list(n = n, p = p, AUC = AUC, precision = precision, lambda = lambda,
              recall = recall, sparsity = sparsity, FPR.ghs = FPR.ghs, TPR.ghs = TPR.ghs))

}

GHS_average <- function(ghs.sf, N){
  
  res = list()
  seed <- 0
  
  n <- dim(ghs.sf$data)[[1]]
  p <- dim(ghs.sf$data)[[2]]
  
  AUC.list = c()
  REC.list = c()
  PREC.list = c()
  SPARS.list = c()
  lambda.list = c()
  
  for(i in 1:N){
    print(i)
     
    set.seed(seed)
    y = mvtnorm::rmvnorm(n, rep(0, p), ghs.sf$sigma)
    y = scale(y)
    
    
    op <- GHS_ROC_sim(ghs.sf, y)
    
    AUC.list = c(AUC.list, op$auc_value)
    REC.list = c(REC.list, op$recall)
    PREC.list = c(PREC.list, op$precision)
    SPARS.list = c(SPARS.list, op$sparsity)
    lambda.list = c(lambda.list, op$lambda)
    seed <- seed + 1
  }
  
  res$AUC = list(mean = mean(AUC.list), SE = sd(AUC.list)/sqrt(N))
  res$REC = list(mean = mean(REC.list), SE = sd(REC.list)/sqrt(N))
  res$PREC = list(mean = mean(PREC.list), SE = sd(PREC.list)/sqrt(N))
  res$SPARS = list(mean = mean(SPARS.list), SE = sd(SPARS.list)/sqrt(N))
  res$AUC = list(mean = mean(AUC.list), SE = sd(AUC.list)/sqrt(N))
  res$lambda = list(mean = mean(lambda.list), SE = sd(AUC.list)/sqrt(N))
  
  res$n = n
  res$p = p
  
  return(res)
  
}

plot.GHS <- function(op){
  plot(op$FPR.ghs,op$TPR.ghs,type='l', col='red', xlab = "FPR", ylab = "TPR", main = "ROC curve for GHS method")
}




