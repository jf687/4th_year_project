source('~/4th_year_project/Sims/ROC_calcs.R')

source('~/4th_year_project/GM/GM.R')
source('~/4th_year_project/GM/fun_GM.R')
source('~/4th_year_project/GM/fun_utils.R')



generate_GM_preds_labels <- function(list_hyper, list_init,n= 200, p = 100, scale.data = T){

  ggm.sf = huge::huge.generator(n=n, d=p,graph = 'scale-free', prob = 0.02) 

  predictions <- c()
  labels <- c()
    

  adj_mat <- ggm.sf$theta
  omega.true <- ggm.sf$omega
  sigma.true <- ggm.sf$sigma
  
  y = mvtnorm::rmvnorm(n, rep(0, p), sigma.true)
  if(scale.data){
    y = scale(y)
  }
  
  res.ssl <- GM(y, list_hyper, list_init)
  
  bool_up <- upper.tri(res.ssl$m_delta)
  
  labels <- c(labels, adj_mat[bool_up])
  predictions <- c(predictions, res.ssl$m_delta[bool_up])
    

  pred <- ROCR::prediction(predictions, labels)
  perf <- ROCR::performance(pred, measure = "tpr", x.measure = "fpr")
  
  auc_value <- ROCR::performance(pred, "auc")@y.values[[1]]
  prec_value <- ROCR::performance(pred, "prec")@y.values[[1]]
  #aic_value <- performance(pred, "aic")@y.values[[1]]
  
  
  return(list(predictions = predictions, labels = labels, pred = pred, perf = perf, auc_value = auc_value, prec_value = prec_value, v0 = list_hyper$v0, v1 = list_hyper$v1))#, aic_value = aic_value))
}


plot_ssl_ROC <- function(res.ssl){
  auc_value = res.ssl$auc_value
  perf = res.ssl$perf
  plot(perf, main = '', xlim = c(0, 1), ylim = c(0, 1), lwd = 2, colorize = TRUE)
  title(main = paste('ROC plot of SSL method. AUC = ',auc_value), sub = paste('v0 = ',list_hyper$v0, 'v1 =', list_hyper$v1))
  abline(a= 0, b = 1, col = "gray80")
  #table(predictions,labels > 0.5)
}


get_ave_AUC <- function(list_hyper, list_init, count = 50, n= 200, p = 100){
  pb <- progress_bar$new(format = "[:bar] :percent ETA: :eta", total = count)
  AUC_vals <- c()
  PREC_vals <- c()
  REC_vals <- c()
  for(p in 1:count){
    pb$tick()
    res.ssl <- generate_GM_preds_labels(list_hyper, list_init)
    AUC_vals <- c(AUC_vals, res.ssl$auc_value)
    PREC_vals <- c(PREC_vals, res.ssl$prec_value)
    REC_vals <- c(REC_vals, res.ssl$rec_value)
    
  }
  return(list(AUC = mean(AUC_vals), PREC = mean(PREC_vals), REC = mean(REC_vals) ))
}

v0_grid <- function(list_hyper, list_init, v0_list = list(0.05, 0.1, 0.5, 1, 5, 10, 50, 100), plot = F, save.data = F,count = 50, n = 200, p = 100){
  
  v0.list = list()
  auc.list = list()
  prec.list = list()
  rec.list = list()
  
  for(v0 in v0_list){
    list_hyper$v0 <- v0
    auc_value <- get_ave_AUC(list_hyper, list_init,count, n,p)$AUC
    prec_value <- get_ave_AUC(list_hyper, list_init,count, n,p)$PREC
    rec_value <- get_ave_AUC(list_hyper, list_init,count, n,p)$REC
    
    if(plot==T){
      res.ssl <- generate_GM_preds_labels(list_hyper, list_init)
      plot_ssl_ROC(res.ssl)
    }
    op <- c(op, list(v0,v1,auc_value))
    v0.list = c(v0.list, v0)
    auc.list = c(auc.list, auc_value)
    prec.list = c(prec.list, prec_value)
    rec.list = c(rec.list, rec_value)
  }
  list_hyper$v0 <- 0.5
  
  op <- list(n, p , v0.list, auc.list, prec.list, rec.list)
  if(save.data == T){
  save(op, file = 'SSL_run.rda')
  }
  return(op)
}

GM_ROC_sim <- function(n = 100, p = 50, plot = F){
  
  # create simulated data Y
  
  L = huge.generator(n=n, d = p, graph="hub")
  Y = mvtnorm::rmvnorm(n, rep(0, p), L$sigma)
  
  # estimate a precision matrix based on Y using mb, ct & glasso methods.
  
  out.mb = huge(L$data)
  out.ct = huge(L$data, method = "ct")
  out.glasso = huge(L$data, method = "glasso")
  
  #model selection using stars. For now we will use out.glasso
  out.select = huge.select(out.glasso, criterion = "stars", stars.thresh = 0.05,rep.num=10)
  plot(out.select)
  
  #set up grid of v0
  v0_v <- seq(1e-2, 1, length.out = 32)
  v1 <- 100
  
  #initialise list parameters
  list_hyper <- list(
    lambda = 2,
    v0_v = v0_v,
    v1 = v1,
    a_tau = 2,
    b_tau = 2,
    a_rho = 1,
    b_rho = p)
  
  
  list_init <-
    list(
      alpha_tau = 1,
      beta_tau = 1,
      alpha_rho = 1,
      beta_rho = p-1
    )
  
  # Perform inference using Graphical SSL model
  res1 <- navigm(Y, method = 'GM', 
                 version = 1,
                 list_hyper = list_hyper,
                 list_init = list_init, 
                 tol = 1e-3,
                 maxit = 1e5)
  #AIC_min <- min(res1$model_criterion$value)
  
  bool_up <- upper.tri(res1$estimates$m_delta)
  
  labels <- L$theta[bool_up]
  predictions <- res1$estimates$m_delta[bool_up]
  
  
  pred <- ROCR::prediction(predictions, labels)
  perf <- ROCR::performance(pred, measure = "tpr", x.measure = "fpr")
  
  auc_value <- ROCR::performance(pred, "auc")@y.values[[1]]
  prec_value <- ROCR::performance(pred, "prec")@y.values[[1]][613]
  browser()
  rec_vaue <- ROCR::performance(pred, "recall")@y.values[[1]][613]
  
  if(plot){
    plot(v0_v[2:16],res1$model_criterion$value[2:16], xlab = 'value of spike variance v0', ylab = 'AIC', type = 'l')
    plot(perf, main = '', xlim = c(0, 1), ylim = c(0, 1), lwd = 2, colorize = TRUE)
  }
  # an example of cutoffs 
  browser()
  cm <- confusion.matrix(labels[bool_up], predictions[bool_up] > 0.5)
  return(list(n=n,p=p,auc_value = auc_value, prec_value = prec_value, rec_value = rec_value, AIC = res1$model_criterion$value, 
              table = table, predictions = predictions, labels = labels))
  
  
}