source('~/4th_year_project/Sims/ROC_calcs.R')



y = mvtnorm::rmvnorm(n, rep(0, p), sigma.true)
y = scale(y)


GLasso_ROC_sim <- function(GLasso.sf, y){
  
  adj_mat <- sf$theta
  omega.true <- sf$omega
  sigma.true <- sf$sigma
  
  res=list()
  
  n <- dim(GLasso.sf$data)[[1]]
  p <- dim(GLasso.sf$data)[[2]]
  
  
  
  n.lambda = 20 # Varying the penalty parameter to get different FDRs
  
#  res$precisions.glasso = matrix(0,1,n.lambda)
#  res$recalls.glasso = matrix(0,1,n.lambda)
#  res$TPR.glasso = matrix(0,1,n.lambda)
#  res$FPR.glasso = matrix(0,1,n.lambda)
#  
  
  res.glasso = huge(y, method = 'glasso', nlambda=n.lambda, verbose = F)
  
  
  sim.obj = huge.select(res.glasso, criterion = 'stars', stars.thresh = 0.05)
  res.glasso.omega.opt = sim.obj$opt.icov

  res$lambda = sim.obj$opt.lambda
  res$precision = precision(abs(res.glasso.omega.opt)>1e-5, abs(omega.true)>1e-5)
  res$recall = recall(abs(res.glasso.omega.opt)>1e-5, abs(omega.true)>1e-5)
  res$TPR = TPR(abs(res.glasso.omega.opt)>1e-5, abs(omega.true)>1e-5)
  res$FPR = FPR(abs(res.glasso.omega.opt)>1e-5, abs(omega.true)>1e-5)
  res$sparsity = sparsity(abs(res.glasso.omega.opt)>1e-5)
  
  bool_up <- upper.tri(res.glasso.omega.opt)
  labels <- adj_mat[bool_up]
  predictions <- res.glasso.omega.opt[bool_up]
  pred <- ROCR::prediction(predictions, labels)
  res$auc_value <- ROCR::performance(pred, "auc")@y.values[[1]]
  
  return(res)
}

GLasso_average <- function(GLasso.sf, N){
  
  res = list()
  seed <- 0
  
  
  n <- dim(GLasso.sf$data)[[1]]
  p <- dim(GLasso.sf$data)[[2]]
  
  AUC.list = c()
  REC.list = c()
  PREC.list = c()
  SPARS.list = c()
  lambda.list = c()
  
  for(i in 1:N){
    print(i)
    
    set.seed(seed)
    y = mvtnorm::rmvnorm(n, rep(0, p), GLasso.sf$sigma)
    y = scale(y)
    
    op <- GLasso_ROC_sim(GLasso.sf, y)
    
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
  res$lambda = list(mean = mean(lambda.list), SE = sd(lambda.list)/sqrt(N))
  
  res$n = n
  res$p = p
  
  return(res)
  
}



plot_ROC = function(sim.obj,  cutoff=NULL){
  df.plot = data.frame(TPR=sim.obj$mean.TPR.glasso, FPR=sim.obj$mean.FPR.glasso)[-1,]
  
  # If TPR = 1 is achieved, can extrapolate
  if(max(df.plot$TPR==1, na.rm=T)){
    df.plot.extra = data.frame(TPR=1, FPR=1)
    df.plot = rbind(df.plot, df.plot.extra)
  }
  if(is.null(cutoff)){
    cutoff=max(df.plot$FPR, na.rm=T)
  }

  ggplot(df.plot, aes(x=FPR, y=TPR))+geom_point(colour='red')+geom_line(colour='red')+theme_bw()+ylim(0,1)+xlim(0,cutoff)+geom_abline(slope=1, linetype='dashed',color='grey')
}




