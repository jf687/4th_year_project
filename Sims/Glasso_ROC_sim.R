source('~/4th_year_project/Sims/ROC_calcs.R')

generate_frequentist_preds_labels <- function(sf, n = 200, p = 100){
  
  adj_mat <- sf$theta
  omega.true <- sf$omega
  sigma.true <- sf$sigma
  
  y = mvtnorm::rmvnorm(n, rep(0, p), sigma.true)
  y = scale(y)
  
  res=list()
  
  
  n.lambda = 20 # Varying the penalty parameter to get different FDRs
  
  res$precisions.glasso = matrix(0,1,n.lambda)
  res$recalls.glasso = matrix(0,1,n.lambda)
  res$TPR.glasso = matrix(0,1,n.lambda)
  res$FPR.glasso = matrix(0,1,n.lambda)
  
  
  res.glasso = huge(y, method = 'glasso', nlambda=n.lambda, verbose = F)
  
  
  sim.obj = huge.select(res.glasso, criterion = 'ebic', ebic.gamma = 0) # stars.thresh = 0.05)
  
  
  res.glasso.omega.opt = sim.obj$opt.icov # A list of precision matrices
  
  ## threshold e-10
  ## find SE
  
  
  
  
  res$precisions.glasso = unlist(lapply(res.glasso.omegas, FUN = function(s) precision(abs(s)>1e-5, abs(omega.true)>1e-5)))
  res$recalls.glasso = unlist(lapply(res.glasso.omegas, FUN = function(s) recall(abs(s)>1e-5, abs(omega.true)>1e-5)))
  res$TPR.glasso[1,] = unlist(lapply(res.glasso.omegas, FUN = function(s) TPR(abs(s)>1e-5, abs(omega.true)>1e-5)))
  res$FPR.glasso[1,] = unlist(lapply(res.glasso.omegas, FUN = function(s) FPR(abs(s)>1e-5, abs(omega.true)>1e-5)))
  res$SPARSITY.glasso = unlist(lapply(res.glasso.omegas, FUN = function(s) sparsity(s)))
  
  
  # Average over all N simulations (i.e., find the average performance for each lambda/PPI)
  
  res$mean.precisions.glasso = mean(res$precisions.glasso)
  res$mean.recalls.glasso = mean(res$recalls.glasso)
  res$mean.TPR.glasso = colMeans(res$TPR.glasso)
  res$mean.FPR.glasso = colMeans(res$FPR.glasso) 

  res$mean.SPARSITY.glasso = mean(res$SPARSITY.glasso)
  
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




freq_run <- function(ns = list(50,100,200,400), ps = list(25,50,100,200,300,400,500), plot= F, save.data = F){
  pb <- progress_bar$new(format = "[:bar] :percent ETA: :eta", total = length(ns)*length(ps))
  
  auc.glasso.list = list()
  prec.glasso.list = list()
  rec.glasso.list = list()
  n.list = list()
  p.list = list()
  sparsity.glasso.list = list()

  sf = huge::huge.generator(n=n, d=p,graph = 'scale-free', prob = 0.02)
  
  adj_mat <- sf$theta
  omega.true <- sf$omega
  sigma.true <- sf$sigma
  
  y = mvtnorm::rmvnorm(n, rep(0, p), sigma.true)
  y = scale(y)
  
  for(n in ns){
    for(p in ps){
      pb$tick()
      res <- generate_frequentist_preds_labels(y, n = 200, p = 100, scale.data = T, include_all = T)
      if(plot){
        plot_ROC(res, include_all = T)
      }
      auc.glasso <- get_AUC(res, method = 'Glasso')
      auc.glasso <- auc.glasso$complete_AUC
      prec.glasso <- mean(res$mean.precisions.glasso)
      rec.glasso <- mean(res$mean.recalls.glasso)
      sparsity.glasso <- mean(res$mean.SPARSITY.glasso)

      n.list <- c(n.list, n)
      p.list <- c(p.list, p)
      auc.glasso.list <- c(auc.glasso.list, auc.glasso)
      prec.glasso.list <- c(prec.glasso.list, prec.glasso)
      rec.glasso.list <- c(rec.glasso.list, rec.glasso)
      sparsity.glasso.list <- c(sparsity.glasso.list, sparsity.glasso)
    }
    
  }
  print('Run complete')
  op <- list(n = n.list, p = p.list, auc = auc.glasso.list, prec = prec.glasso.list, rec =  rec.glasso.list, sparsity = sparsity.glasso.list)
  if(save.data){
  save(op, file='op.rda')
  }
  return(op)
}

colMax <- function (colData) {
  apply(colData, MARGIN=c(2), max)
}


GLasso_multi_rep <- function(N = 50, n = 200, p = 100){
  
  auc.glasso.list = list()
  prec.glasso.list = list()
  rec.glasso.list = list()
  sparsity.glasso.list = list()
  
  
  sf = huge::huge.generator(n=n, d=p,graph = 'scale-free', prob = 0.02)
  
  adj_mat <- sf$theta
  omega.true <- sf$omega
  sigma.true <- sf$sigma
  
  y = mvtnorm::rmvnorm(n, rep(0, p), sigma.true)
  y = scale(y)
  
  for(l in 1:N){
    
    print(l)
    
    res <- generate_frequentist_preds_labels(sf, n, p)
    
  
  
    auc.glasso <- get_AUC(res, method = 'Glasso')
    auc.glasso <- auc.glasso$complete_AUC
    prec.glasso <- mean(res$mean.precisions.glasso)
    rec.glasso <- mean(res$mean.recalls.glasso)
    sparsity.glasso <- mean(res$mean.SPARSITY.glasso)
    
  
    auc.glasso.list <- c(auc.glasso.list, auc.glasso)
    prec.glasso.list <- c(prec.glasso.list, prec.glasso)
    rec.glasso.list <- c(rec.glasso.list, rec.glasso)
    sparsity.glasso.list <- c(sparsity.glasso.list, sparsity.glasso)
    
    if(l == 1){
      browser()
      cutoff <- max(res$mean.FPR.glasso)
      df.plot = data.frame(TPR=res$mean.TPR.glasso, FPR=res$mean.FPR.glasso)[-1,]
      ggplot(df.plot, aes(x=FPR, y=TPR))+geom_point(colour='red')+geom_line(colour='red')+theme_bw()+ylim(0,1)+xlim(0,cutoff)+geom_abline(slope=1, linetype='dashed',color='grey')
    }
  }

  return(list(n = n, p = p, AUC = mean(auc.glasso.list[[1]]), precision = mean(prec.glasso.list[[1]]), recall = mean(rec.glasso.list[[1]]), sparsity = mean(sparsity.glasso.list[[1]])))
                            
                            
}