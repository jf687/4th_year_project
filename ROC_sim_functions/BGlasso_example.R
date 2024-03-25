
source('~/4th_year_project/aux_sim_funcs/ROC_calcs.R')


#data.sf = huge::huge.generator(n=n, d=p,graph = 'scale-free')
# Example of how to use the Bayesian graphical lasso by Wang et al. (2012)
#y.sf = mvtnorm::rmvnorm(n, sigma=solve(theta.true)) # Sample data
#y = scale(y.sf) # Scale columns/variables.


BGLasso_ROC_sim <- function(data.sf, y){

  theta.true = data.sf$omega # The precision matrix
  theta.true[which(theta.true<10e-5,arr.ind=T)]=0
  
  n <- dim(data.sf$data)[[1]]
  p <- dim(data.sf$data)[[2]]

  iterations = 1000
  burnIn = 100
  
  res.bglasso <- BayesianGLasso::blockGLasso(y, iterations = 1000, burnIn = 100, lambdaPriora = 1, lambdaPriorb = 1/10)
  res.bglasso$Omegas = simplify2array(res.bglasso$Omegas)
  theta.est.bglasso = apply(simplify2array(res.bglasso$Omegas), c(1,2), mean) # Posterior mean of MCMC samples

  lambda <- mean(res.bglasso$lambdas[burnIn:burnIn + iterations])
  
  # To threshold, Li et al. (2017) use the symmetric central 50% posterior credible intervals for variable selection. 
  # That is, if the 50% posterior credible interval of an off-diagonal element of does not contain zero, that element is considered a discovery, and vice versa. 
  theta.bglasso.credible.lower = apply(res.bglasso$Omegas, c(1,2), FUN = function(s) quantile(s, 0.25))
  theta.bglasso.credible.upper = apply(res.bglasso$Omegas, c(1,2), FUN = function(s) quantile(s, 0.75))
  theta.bglasso.selected = theta.est.bglasso
  theta.bglasso.selected[(theta.bglasso.credible.lower < 0 & theta.bglasso.credible.upper > 0)] = 0
  
  precision = precision(theta.true!=0, theta.bglasso.selected!=0)
  recall = recall(theta.true!=0, theta.bglasso.selected!=0)

  diag(theta.bglasso.selected) <- 0
  sparsity = sparsity(theta.bglasso.selected)
  
  # To vary the threshold for ROC curves etc, you can increase/decrease the width of the credible intervals used to determine nonzero effects. 
  # I.e., replace 0.25 and 0.75 by x and 1.0-x for any 0<x<0.5.
  quant.thresh = seq(0.25,0.75,length.out=50)
  theta.path.BGLasso = list()
  TPR.BGLasso = rep(0,length(quant.thresh))
  FPR.BGLasso = rep(0,length(quant.thresh))
  
  pb <- progress_bar$new(format = "[:bar] :percent ETA: :eta", total = length(quant.thresh))
  
  for(x in 1:length(quant.thresh)){
    
    pb$tick()
    theta.bglasso.credible.lower = apply(res.bglasso$Omegas, c(1,2), FUN = function(s) quantile(s, quant.thresh[x]))
    theta.bglasso.credible.upper = apply(res.bglasso$Omegas, c(1,2), FUN = function(s) quantile(s, 1-quant.thresh[x]))
    theta.bglasso.selected = theta.est.bglasso
    theta.bglasso.selected[(theta.bglasso.credible.lower < 0 & theta.bglasso.credible.upper > 0)] = 0
    theta.path.BGLasso[[x]] = theta.bglasso.selected
    TPR.BGLasso[x] = TPR(theta.true!=0, theta.bglasso.selected!=0)
    FPR.BGLasso[x] = FPR(theta.true!=0, theta.bglasso.selected!=0)
  }
  
  x = c(0,FPR.BGLasso)
  y = c(0,TPR.BGLasso)
  AUC=flux::auc(x, y)
  
  return(list(n = n, p = p, AUC = AUC, precision = precision, recall = recall, 
              sparsity = sparsity, lambda = lambda, FPR.BGLasso = FPR.BGLasso, TPR.BGLasso = TPR.BGLasso))
  
}

BGLasso_average <- function(data.sf, N){
  
  res = list()
  seed <- 0
  
  n <- dim(data.sf$data)[[1]]
  p <- dim(data.sf$data)[[2]]
  
  AUC.list = c()
  REC.list = c()
  PREC.list = c()
  SPARS.list = c()
  lambda.list = c()
  
  for(i in 1:N){
    
    set.seed(seed)
    y = mvtnorm::rmvnorm(n, rep(0, p), data.sf$sigma)
    y = scale(y)
    
    op <- BGLasso_ROC_sim(data.sf, y)
    
    AUC.list = c(AUC.list, op$AUC)
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

plot.BGLasso <- function(op){
  plot(op$FPR.BGLasso,op$TPR.BGLasso,type='l', col='red', xlab = "FPR", ylab = "TPR", main = "ROC curve for BGLasso method")
}

