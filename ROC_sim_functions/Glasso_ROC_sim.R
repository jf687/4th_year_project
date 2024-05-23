source('~/4th_year_project/aux_sim_funcs/ROC_calcs.R')

na.rm = T

GLasso_ROC_sim <- function(GLasso.sf, y){
  
  adj_mat <- sf$theta
  omega.true <- sf$omega
  sigma.true <- sf$sigma
  
  res=list()
  
  n <- dim(GLasso.sf$data)[[1]]
  p <- dim(GLasso.sf$data)[[2]]
  
  
  
  n.lambda = 50 # Varying the penalty parameter to get different FDRs
  
  start.time <- Sys.time()
  res.glasso = huge(y, method = 'glasso', nlambda=n.lambda, verbose = F)
  end.time <- Sys.time()
  
  sim.obj = huge.select(res.glasso, criterion = 'stars', stars.thresh = 0.05)
  res.glasso.omega.opt = sim.obj$opt.icov
  all.FPRs = list()
  for(i in length(res.glasso$icov)){
    TPRs[[i]] = TPR(theta.true!=0, res.glasso$icov!=0)
    FPRs[[i]] = FPR(theta.true!=0, res.glasso$icov!=0)
  }

  res$lambda = sim.obj$opt.lambda
  res$precision = precision(abs(omega.true)>1e-15, abs(res.glasso.omega.opt)>1e-15)
  res$recall = recall(abs(omega.true)>1e-15, abs(res.glasso.omega.opt)>1e-15)
  res$TPR = TPR(abs(omega.true)>1e-15, abs(res.glasso.omega.opt)>1e-15)
  res$FPR = FPR(abs(omega.true)>1e-15, abs(res.glasso.omega.opt)>1e-15)
  res$sparsity = sparsity(abs(res.glasso.omega.opt)>1e-15)
  res$time = end.time - start.time
  
  bool_up <- upper.tri(res.glasso.omega.opt)
  labels <- adj_mat[bool_up]
  predictions <- res.glasso.omega.opt[bool_up]
  
  pred <- ROCR::prediction(predictions, labels)
  res$perf <- ROCR::performance(pred, measure = "tpr", x.measure = "fpr")
  res$auc_value <- ROCR::performance(pred, "auc")@y.values[[1]]
  return(res)
}

GLasso_average <- function(GLasso.sf, N, plot = FALSE){
  
  res = list()
  seed <- 0
  
  
  n <- dim(GLasso.sf$data)[[1]]
  p <- dim(GLasso.sf$data)[[2]]
  
  AUC.list = c()
  REC.list = c()
  PREC.list = c()
  SPARS.list = c()
  lambda.list = c()
  time.list = c()
  
  for(i in 1:N){
    print(i)
    
    set.seed(seed)
    y = mvtnorm::rmvnorm(n, rep(0, p), GLasso.sf$sigma)
    y = scale(y)
    
    op <- GLasso_ROC_sim(GLasso.sf, y)
    
    if(plot & i == 1 ){
      print("trying to plot GLasso ROC")
      plot_ROC(op)
    }
    
    AUC.list = c(AUC.list, op$auc_value)
    REC.list = c(REC.list, op$recall)
    PREC.list = c(PREC.list, op$precision)
    SPARS.list = c(SPARS.list, op$sparsity)
    lambda.list = c(lambda.list, op$lambda)
    time.list = c(time.list, op$time)
    
    seed <- seed + 1
  }
  
  res$AUC = list(mean = mean(AUC.list), SE = sd(AUC.list)/sqrt(N))
  res$REC = list(mean = mean(REC.list), SE = sd(REC.list)/sqrt(N))
  res$PREC = list(mean = mean(PREC.list), SE = sd(PREC.list)/sqrt(N))
  res$SPARS = list(mean = mean(SPARS.list), SE = sd(SPARS.list)/sqrt(N))
  res$AUC = list(mean = mean(AUC.list), SE = sd(AUC.list)/sqrt(N))
  res$lambda = list(mean = mean(lambda.list), SE = sd(lambda.list)/sqrt(N))
  res$time = list(mean = mean(time.list), SE = sd(time.list)/sqrt(N))
  
  res$n = n
  res$p = p
  
  return(res)
  
}



plot_ROC = function(sim.obj,  cutoff=NULL){
  # sim.obj is an output from the GLasso_ROC_sim function
  df.plot = data.frame(TPR=sim.obj$mean.TPR.glasso, FPR=sim.obj$mean.FPR.glasso)[-1,]
  browser()
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

plot_GLasso_ROC = function(res){
  # res needs to be result from glasso_roc_sim
  auc_value = res$auc_value
  perf = res$perf
  plot(perf, main = '', xlim = c(0, 1), ylim = c(0, 1), lwd = 2, colorize = TRUE)
  title(main = runif(1))
  #title(main = paste('ROC plot using GLasso method with data size n = 200, p = 100.'), sub = bquote(v[0] = v0 v[1] = list_hyper$v1)))
  abline(a= 0, b = 1, col = "gray80")
}

test_GLasso = function(){
  set.seed(999)
  sf = huge::huge.generator(n=200, d=100,graph = 'scale-free', prob = 0.02)
  y = mvtnorm::rmvnorm(n, rep(0, p), sf$sigma)
  y = scale(y)
  opGL <- GLasso_ROC_sim(sf,y)
  plot_ROC(opGL)
  plot_GLasso_ROC(opGL)
}


