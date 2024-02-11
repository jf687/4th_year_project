## Code for generating cut-off ROC and precision-recall curves

source('libraries.R')
source('simulations.R')
source('ROC_calcs.R')

# Function for plotting ROC curves

perform_ROC_simulation = function(omega.true, n,list_hyper, list_init, N=100, include.glasso=T, include.ssl=T, scale.data = T){
  adj_mat <- adj_gen(0.02,p)
  while(all(adj_mat ==0)){
    print('¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬')
    adj_mat <- adj_gen(0.02,p)
  }
  omega.true <- prec_from_adj(adj_mat,p)
  sigma.true <- solve(omega.true)
  
  res=list()
  p=nrow(omega.true)
  sigma.true = solve(omega.true)
  
  if(include.glasso){
    n.lambda = 20 # Varying the penalty parameter to get different FDRs
    res$precisions.glasso = matrix(0,N,n.lambda)
    res$recalls.glasso = matrix(0,N,n.lambda)
    res$TPR.glasso = matrix(0,N,n.lambda)
    res$FPR.glasso = matrix(0,N,n.lambda)
    
    
    
    for(i in 1:N){
      y = mvtnorm::rmvnorm(n, rep(0, p), sigma.true)
      if(scale.data){
        y = scale(y)
      }
      res.glasso = huge(y, method = 'glasso', nlambda=n.lambda, verbose = F)
      res.glasso.omegas = res.glasso$icov # A list of precision matrices
      res$precisions.glasso[i,] = unlist(lapply(res.glasso.omegas, FUN = function(s) precision(abs(s)>1e-5, abs(omega.true)>1e-5)))
      res$recalls.glasso[i,] = unlist(lapply(res.glasso.omegas, FUN = function(s) recall(abs(s)>1e-5, abs(omega.true)>1e-5)))
      res$TPR.glasso[i,] = unlist(lapply(res.glasso.omegas, FUN = function(s) TPR(abs(s)>1e-5, abs(omega.true)>1e-5)))
      res$FPR.glasso[i,] = unlist(lapply(res.glasso.omegas, FUN = function(s) FPR(abs(s)>1e-5, abs(omega.true)>1e-5)))
    }
    
    
    # Average over all N simulations (i.e., find the average performance for each lambda/PPI)
    # Glasso
    res$mean.precisions.glasso = colMeans(res$precisions.glasso)
    res$mean.recalls.glasso = colMeans(res$recalls.glasso)
    res$mean.TPR.glasso = colMeans(res$TPR.glasso)
    res$mean.FPR.glasso = colMeans(res$FPR.glasso) 
  }
  
  
  if(include.ssl){
    n.ppi.thresh = 50 # Threshold the posterior inclusion probability to get different FDR (i.e. not thresholding matrix elements)
    #ppi.thresh = seq(0,1,n.ppi.thresh)
    res$precisions.ssl = matrix(0,N,n.ppi.thresh)
    res$recalls.ssl = matrix(0,N,n.ppi.thresh)
    res$TPR.ssl = matrix(0,N,n.ppi.thresh)
    res$FPR.ssl = matrix(0,N,n.ppi.thresh)
    
    pb <- progress_bar$new(format = "[:bar] :percent ETA: :eta", total = N*100)
    
    thresh_lo <- 0
    thresh_hi <- 1
    for(i in 1:N){
      pb$tick()
      adj_mat <- adj_gen(0.02,p)
      while(all(adj_mat ==0)){
        print('¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬')
        adj_mat <- adj_gen(0.02,p)
      }
      new.omega.true <- prec_from_adj(adj_mat)
      new.sigma.true <- solve(new.omega.true)
      #new.omega.true <- GM$Omega
      y = mvtnorm::rmvnorm(n, rep(0, p), new.sigma.true)
      if(scale.data){
        y = scale(y)
      }
      
      
      res.ssl <- GM(y, list_hyper, list_init) # Must specify the hyper params as well
      res.ssl.omega <- res.ssl$Omega
      res.ssl.m_delta <- res.ssl$m_delta
      
      thresh_lo <- min(res.ssl.m_delta[lower.tri(res.ssl.m_delta)])
      thresh_hi <- max(res.ssl.m_delta[lower.tri(res.ssl.m_delta)])
      
      
      ppi.thresh = seq(thresh_lo, thresh_hi, length.out = n.ppi.thresh)
      
      for(x in length(ppi.thresh)){
      
        res$precisions.ssl[i,x] <- precision( abs( res.ssl.m_delta ) < ppi.thresh[x], abs( new.omega.true ) > 0)
        res$recalls.ssl[i,x] <- recall( abs( res.ssl.m_delta ) < ppi.thresh[x], abs( new.omega.true ) >0)
        res$TPR.ssl[i,x] <- TPR( abs( res.ssl.m_delta ) < ppi.thresh[x], abs( new.omega.true ) >0)
        res$FPR.ssl[i,x] <- FPR( abs( res.ssl.m_delta ) < ppi.thresh[x], abs( new.omega.true )>0)
      }
    #browser()
      
    }
    
    
    
    res$mean.precisions.ssl = colMeans(res$precisions.ssl)
    res$mean.recalls.ssl = colMeans(res$recalls.ssl)
    res$mean.TPR.ssl = colMeans(res$TPR.ssl)
    res$mean.FPR.ssl = colMeans(res$FPR.ssl) 
  
  }
  
  bool_up <- upper.tri(res.ssl$m_delta)
  pred <- prediction(res.ssl$m_delta[bool_up], adj_mat[bool_up])
  #browser()
  #range(pred@cutoffs[[1]])
  #quantile(pred@cutoffs[[1]])
  perf <- performance(pred, measure = "prec", x.measure = "rec")
  plot(perf, main = '', xlim = c(0, 1), ylim = c(0, 1), lwd = 2, colorize = TRUE)
  perf <- performance(pred, measure = "tpr", x.measure = "fpr")
  plot(perf, main = '', xlim = c(0, 1), ylim = c(0, 1), lwd = 2, colorize = TRUE)
  abline(a= 0, b = 1, col = "gray80")
  table(adj_mat[bool_up],res.ssl$m_delta[bool_up] > 0.5)

  
  return(res)
}

plot_ROC = function(sim.obj, include.glasso=T, include.ssl=T, cutoff=NULL){
  include.both = include.glasso & include.ssl
  if(include.glasso){
    df.plot = data.frame(TPR=sim.obj$mean.TPR.glasso, FPR=sim.obj$mean.FPR.glasso)[-1,]
  }
  if(include.ssl){
    if(include.both){
      df.plot2 = data.frame(TPR=sim.obj$mean.TPR.ssl, FPR=sim.obj$mean.FPR.ssl)
      df.plot=rbind(df.plot, df.plot2)
      df.plot$method = c(rep('Glasso',length(sim.obj$mean.TPR.glasso)-1), rep('SSL',length(sim.obj$mean.TPR.ssl)))
    }
    else{
      df.plot = data.frame(TPR=sim.obj$mean.TPR.ssl, FPR=sim.obj$mean.FPR.ssl) 
    }
  }
  # If TPR = 1 is achieved, can extrapolate
  if(! include.both & max(df.plot$TPR==1, na.rm=T)){
    df.plot.extra = data.frame(TPR=1, FPR=1)
    df.plot = rbind(df.plot, df.plot.extra)
  }
  if(include.both){
    if(max(df.plot$TPR==1 & df.plot$method=='Glasso', na.rm=T)){
      df.plot.extra = data.frame(TPR=1, FPR=1, method='Glasso')
      df.plot = rbind(df.plot, df.plot.extra)
    }
    if(max(df.plot$TPR==1 & df.plot$method=='SSL', na.rm=T)){
      df.plot.extra = data.frame(TPR=1, FPR=1, method='SSL')
      df.plot = rbind(df.plot, df.plot.extra)
    }
  }
  if(is.null(cutoff)){
    cutoff=max(df.plot$FPR, na.rm=T)
  }
  if(include.both){
    ggplot(df.plot, aes(x=FPR, y=TPR, colour=method))+geom_point()+geom_line()+theme_bw()+ylim(0,1)+xlim(0,1)+geom_abline(slope=1, linetype='dashed',color='grey')
  }
  else{
    ggplot(df.plot, aes(x=FPR, y=TPR))+geom_point(colour='red')+geom_line(colour='red')+theme_bw()+ylim(0,1)+xlim(0,cutoff)+geom_abline(slope=1, linetype='dashed',color='grey')
  }
}

plot_PRC = function(sim.obj, include.glasso=T, include.ssl=T, cutoff=NULL){
  # Plot precision-recall curve
  include.both = include.glasso & include.ssl
  if(include.glasso){
    df.plot = data.frame(Precision=sim.obj$mean.precisions.glasso, Recall=sim.obj$mean.recalls.glasso)[-1,]
  }
  if(include.ssl){
    if(include.both){
      df.plot2 = data.frame(Precision=sim.obj$mean.precisions.ssl, Recall=sim.obj$mean.recalls.ssl)
      df.plot=rbind(df.plot, df.plot2)
      df.plot$method = c(rep('Glasso',length(sim.obj$mean.precisions.glasso)-1), rep('SSL',length(sim.obj$mean.precisions.ssl)))
    }
    else{
      df.plot = data.frame(Precision=sim.obj$mean.precisions.ssl, Recall=sim.obj$mean.recall.ssl)
    }
  }
  if(is.null(cutoff)){
    cutoff=max(df.plot$Recall, na.rm=T)
  }
  if(include.both){
    ggplot(df.plot, aes(x=Recall, y=Precision, colour=method))+geom_point()+geom_line()+theme_bw()+ylim(0,1)+xlim(0,cutoff)+geom_abline(slope=-1, intercept = 1,linetype='dashed',color='grey')
  }
  else{
    ggplot(df.plot, aes(x=Recall, y=Precision))+geom_point(colour='red')+geom_line(colour='red')+theme_bw()+ylim(0,1)+xlim(0,cutoff)+geom_abline(slope=-1,intercept = 1, linetype='dashed',color='grey')
  }
}

get_AUC = function(sim.obj, method='Glasso', cutoff=NULL){
  # Adding the point (0,0)
  if(method=='Glasso'){
    x = c(0,sim.obj$mean.FPR.glasso[-1])
    y = c(0,sim.obj$mean.TPR.glasso[-1])
  }
  else if(method=='SSL'){
    x = c(0,sim.obj$mean.FPR.ssl)
    y = c(0,sim.obj$mean.TPR.ssl)   
  }
  else if(method=='GM')
  # Extrapolate if TPR = 1 is achieved
  if(max(y)==1){
    x = c(x,1)
    y = c(y,1)
  }
  if(is.null(cutoff)){
    cutoff=max(x)
  }
  res= list(AUC=flux::auc(x, y), cutoff=cutoff, complete_AUC = 1-cutoff + flux::auc(x, y))
  return(res)
}

get_PRC_AUC = function(sim.obj, method='Glasso', cutoff=NULL){
  # Adding the point (0,0)
  if(method=='Glasso'){
    x = c(0,sim.obj$mean.recalls.glasso[-1])
    y = c(0,sim.obj$mean.precisions.glasso[-1])
  }
  else if(method=='SSL'){
    x = c(0,sim.obj$mean.recalls.ssl)
    y = c(0,sim.obj$mean.precisions.ssl)   
  }
  if(is.null(cutoff)){
    cutoff=max(x)
  }
  res= list(AUC=flux::auc(x, y), cutoff=cutoff)
  return(res)
}



