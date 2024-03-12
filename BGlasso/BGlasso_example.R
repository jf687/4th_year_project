rm(list=ls())
source('libraries.R')
source('Sims/ROC_calcs.R')

#data.sf = huge::huge.generator(n=n, d=p,graph = 'scale-free')
# Example of how to use the Bayesian graphical lasso by Wang et al. (2012)
BGlasso_ROC_sim <- function(data.sf, n=200,p=100,plot = F){
  set.seed(123)
  theta.true = data.sf$omega # The precision matrix
  theta.true[which(theta.true<10e-5,arr.ind=T)]=0
  y.sf = mvtnorm::rmvnorm(n, sigma=solve(theta.true)) # Sample data
  y.sf.scaled = scale(y.sf) # Scale columns/variables.
  data.sf$sparsity # True sparsity: 0.02
  
  # Using the standard prior choices as proposed by Wang et al. (2012) in their package
  res.bglasso = BayesianGLasso::blockGLasso(y.sf.scaled,iterations = 1000, burnIn = 100, lambdaPriora = 1, lambdaPriorb = 1/10)
  theta.est.bglasso = apply(simplify2array(res.bglasso$Omegas), c(1,2), mean) # Posterior mean of MCMC samples
  
  # Threshold the same way as GHS, varying the threshold for the credible interval to use. 
  
  
  bool_up <- upper.tri(theta.est.bglasso)
  
  labels <- theta.true[bool_up]>0
  predictions <- theta.est.bglasso[bool_up]
  
  
  pred <- ROCR::prediction(predictions, labels)
  perf <- ROCR::performance(pred, measure = "tpr", x.measure = "fpr")

  
  thresh <- return_closest_threshold(pred@cutoffs[[1]])
  auc_value <- ROCR::performance(pred, "auc")@y.values[[1]]
  prec_value <- ROCR::performance(pred, "prec")@y.values[[1]][thresh]
  rec_value <- ROCR::performance(pred, "rec")@y.values[[1]][thresh]
  
  if(plot){
    plot(perf, main = '', xlim = c(0, 1), ylim = c(0, 1), lwd = 2, colorize = TRUE)
  }
  
  return(list(n = n, p = p, auc_value = auc_value, prec_value = prec_value, rec_value = rec_value, sparsity = data.sf$sparsity))
}

multiple_BGlasso <- function(N, n_, p_){
  AUC <- 0
  REC <- 0
  PREC <- 0
  SPARS <- 0
  
  data.sf = huge::huge.generator(n=n_, d=p_, graph = 'scale-free')
  for(i in 1:N){
    op <- BGlasso_ROC_sim(data.sf, n=n_,p=p_,plot = F)
    AUC <- AUC + op$auc_value
    REC <- REC + op$rec_value
    PREC <- PREC + op$prec_value
    SPARS <- SPARS + op$sparsity
  }
  res$n = n_
  res$p = p_
  res$AUC = AUC/N
  res$REC = REC/N
  res$PREC = PREC/N
  res$SPARS =  SPARS/N
  
  return(res)
  
}

BGlasso_run <- function(ns, ps){
  pb <- progress_bar$new(format = "[:bar] :percent ETA: :eta", total = length(ns)*length(ps))
  
  n.list = list()
  p.list = list()
  auc.list = list()
  rec.list = list()
  prec.list = list()
  sparsity.list = list()
  
  for(n in ns){
    for(p in ps){
      pb$tick()
      BGlasso <- multiple_BGlasso(20,n,p) #BGlasso_ROC_sim(n = n, p = p)
      n.list <- c(n.list, n)
      p.list <- c(p.list, p)
      auc.list <- c(auc.list, BGlasso$AUC)
      rec.list <- c(rec.list, BGlasso$REC)
      prec.list <- c(rec.list, BGlasso$PREC)
      sparsity.list <- c(sparsity.list, BGlasso$SPARS)
      
      
    }
  }
  
  out$ns = n.list
  out$ps = p.list
  out$aucs = auc.list
  out$recs = rec.list
  out$precs = prec.list
  out$sparsities = sparsity.list
  return(out)
}
