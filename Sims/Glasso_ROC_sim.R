source('~/4th_year_project/Sims/ROC_calcs.R')

generate_frequentist_preds_labels <- function(n = 200, p = 100, scale.data = T, include_all = T){
  
  
  res=list()
  res$method = method
  
  pb <- progress_bar$new(format = "[:bar] :percent ETA: :eta", total = N)
  n.lambda = 20 # Varying the penalty parameter to get different FDRs
  
  res$precisions.glasso = matrix(0,N,n.lambda)
  res$recalls.glasso = matrix(0,N,n.lambda)
  res$TPR.glasso = matrix(0,N,n.lambda)
  res$FPR.glasso = matrix(0,N,n.lambda)
  
  if(include_all){
  res$precisions.ct = matrix(0,N,n.lambda)
  res$recalls.ct = matrix(0,N,n.lambda)
  res$TPR.ct = matrix(0,N,n.lambda)
  res$FPR.ct = matrix(0,N,n.lambda)

  res$precisions.mb = matrix(0,N,n.lambda)
  res$recalls.mb = matrix(0,N,n.lambda)
  res$TPR.mb = matrix(0,N,n.lambda)
  res$FPR.mb = matrix(0,N,n.lambda)
  
  res$precisions.tiger = matrix(0,N,n.lambda)
  res$recalls.tiger = matrix(0,N,n.lambda)
  res$TPR.tiger = matrix(0,N,n.lambda)
  res$FPR.tiger = matrix(0,N,n.lambda)
  }
  
  
  
  
  for(i in 1:N){
    
    adj_mat <- adj_gen(0.02,p)
    while(all(adj_mat ==0)){
      adj_mat <- adj_gen(0.02,p)
    }
    omega.true <- prec_from_adj(adj_mat,p)
    sigma.true <- solve(omega.true)
    
    y = mvtnorm::rmvnorm(n, rep(0, p), sigma.true)
    pb$tick()
    if(scale.data){
      y = scale(y)
    }
    
    res.glasso = huge(y, method = method, nlambda=n.lambda, verbose = F)
    res.glasso.omegas = res.glasso$icov # A list of precision matrices
    res$precisions.glasso[i,] = unlist(lapply(res.glasso.omegas, FUN = function(s) precision(abs(s)>1e-5, abs(omega.true)>1e-5)))
    res$recalls.glasso[i,] = unlist(lapply(res.glasso.omegas, FUN = function(s) recall(abs(s)>1e-5, abs(omega.true)>1e-5)))
    res$TPR.glasso[i,] = unlist(lapply(res.glasso.omegas, FUN = function(s) TPR(abs(s)>1e-5, abs(omega.true)>1e-5)))
    res$FPR.glasso[i,] = unlist(lapply(res.glasso.omegas, FUN = function(s) FPR(abs(s)>1e-5, abs(omega.true)>1e-5)))
    
    if(include_all){
      
      #res.ct = huge(y, method = 'ct', nlambda = n.lambda, verbose = F, scr = T)
      #res.ct.omegas = res.ct$icov # A list of precision matrices
      #browser()
      #res$precisions.ct[i,] = unlist(lapply(res.ct.omegas, FUN = function(s) precision(abs(s)>1e-5, abs(omega.true)>1e-5)))
      #res$recalls.ct[i,] = unlist(lapply(res.ct.omegas, FUN = function(s) recall(abs(s)>1e-5, abs(omega.true)>1e-5)))
      #res$TPR.ct[i,] = unlist(lapply(res.ct.omegas, FUN = function(s) TPR(abs(s)>1e-5, abs(omega.true)>1e-5)))
      #res$FPR.ct[i,] = unlist(lapply(res.ct.omegas, FUN = function(s) FPR(abs(s)>1e-5, abs(omega.true)>1e-5)))
      
      #res.mb = huge(y, method = 'mb', nlambda = n.lambda, verbose = F, scr = T)
      #res.mb.omegas = res.mb$icov # A list of precision matrices
      #res$precisions.mb[i,] = unlist(lapply(res.mb.omegas, FUN = function(s) precision(abs(s)>1e-5, abs(omega.true)>1e-5)))
      #res$recalls.mb[i,] = unlist(lapply(res.mb.omegas, FUN = function(s) recall(abs(s)>1e-5, abs(omega.true)>1e-5)))
      #res$TPR.mb[i,] = unlist(lapply(res.mb.omegas, FUN = function(s) TPR(abs(s)>1e-5, abs(omega.true)>1e-5)))
      #res$FPR.mb[i,] = unlist(lapply(res.mb.omegas, FUN = function(s) FPR(abs(s)>1e-5, abs(omega.true)>1e-5)))
      
      
      res.tiger = huge(y, method = 'tiger', nlambda = n.lambda, verbose = F)
      res.tiger.omegas = res.tiger$icov # A list of precision matrices
      res$precisions.tiger[i,] = unlist(lapply(res.tiger.omegas, FUN = function(s) precision(abs(s)>1e-5, abs(omega.true)>1e-5)))
      res$recalls.tiger[i,] = unlist(lapply(res.tiger.omegas, FUN = function(s) recall(abs(s)>1e-5, abs(omega.true)>1e-5)))
      res$TPR.tiger[i,] = unlist(lapply(res.tiger.omegas, FUN = function(s) TPR(abs(s)>1e-5, abs(omega.true)>1e-5)))
      res$FPR.tiger[i,] = unlist(lapply(res.tiger.omegas, FUN = function(s) FPR(abs(s)>1e-5, abs(omega.true)>1e-5)))
    }
    
  }
  
  
  # Average over all N simulations (i.e., find the average performance for each lambda/PPI)
  # Glasso
  res$mean.precisions.glasso = colMeans(res$precisions.glasso)
  res$mean.recalls.glasso = colMeans(res$recalls.glasso)
  res$mean.TPR.glasso = colMeans(res$TPR.glasso)
  res$mean.FPR.glasso = colMeans(res$FPR.glasso) 
  
  if(include_all){
  #res$mean.precisions.ct = colMeans(res$precisions.ct)
  #res$mean.recalls.ct = colMeans(res$recalls.ct)
  #res$mean.TPR.ct = colMeans(res$TPR.ct)
  #res$mean.FPR.ct = colMeans(res$FPR.ct) 
  
  #res$mean.precisions.mb = colMeans(res$precisions.mb)
  #res$mean.recalls.mb = colMeans(res$recalls.mb)
  #res$mean.TPR.mb = colMeans(res$TPR.mb)
  #res$mean.FPR.mb = colMeans(res$FPR.mb) 
  
  res$mean.precisions.tiger = colMeans(res$precisions.tiger)
  res$mean.recalls.tiger = colMeans(res$recalls.tiger)
  res$mean.TPR.tiger = colMeans(res$TPR.tiger)
  res$mean.FPR.tiger = colMeans(res$FPR.tiger) 
  }
  return(res)
}

plot_ROC = function(sim.obj, include_all = T,  cutoff=NULL){
  df.plot = data.frame(TPR=sim.obj$mean.TPR.glasso, FPR=sim.obj$mean.FPR.glasso)[-1,]
  
  if(include_all){
    #df.plot2 = data.frame(TPR=sim.obj$mean.TPR.ct, FPR=sim.obj$mean.FPR.ct)
    #df.plot3 = data.frame(TPR=sim.obj$mean.TPR.mb, FPR=sim.obj$mean.FPR.mb)
    df.plot2 = data.frame(TPR=sim.obj$mean.TPR.tiger, FPR=sim.obj$mean.FPR.tiger)[-1,]
    
    df.plot=rbind(df.plot, df.plot2)
    df.plot$method = c(rep('Glasso',length(sim.obj$mean.TPR.glasso)-1), rep('tiger',length(sim.obj$mean.TPR.tiger)-1))
    browser()
    }
  # If TPR = 1 is achieved, can extrapolate
  if(! include_all & max(df.plot$TPR==1, na.rm=T)){
    df.plot.extra = data.frame(TPR=1, FPR=1)
    df.plot = rbind(df.plot, df.plot.extra)
  }
  browser()
  if(include_all){
    if(max(df.plot$TPR==1, na.rm = T)){# & df.plot$method=='Glasso', na.rm=T)){
      df.plot.extra = data.frame(TPR=1, FPR=1, method='Glasso')
      df.plot = rbind(df.plot, df.plot.extra)
    #}
    #if(max(df.plot$TPR==1 & df.plot$method=='tiger', na.rm=T)){
      df.plot.extra <- data.frame(TPR=1, FPR=1, method='tiger')
      df.plot <-  rbind(df.plot, df.plot.extra)
    }
  }
  if(is.null(cutoff)){
    cutoff=max(df.plot$FPR, na.rm=T)
  }
  if(include_all){
    ggplot(df.plot, aes(x=FPR, y=TPR, colour=method))+geom_point()+geom_line()+theme_bw()+ylim(0,1)+xlim(0,1)+geom_abline(slope=1, linetype='dashed',color='grey')
  }
  else{
    ggplot(df.plot, aes(x=FPR, y=TPR))+geom_point(colour='red')+geom_line(colour='red')+theme_bw()+ylim(0,1)+xlim(0,cutoff)+geom_abline(slope=1, linetype='dashed',color='grey')
  }
}