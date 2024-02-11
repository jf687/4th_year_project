source('libraries.R')
source('GM.R')
source('fun_GM.R')
source('ROC_calcs.R')


generate_GM_preds_labels <- function(list_hyper, list_init,n= 200, p = 100, scale.data = T){

  predictions <- c()
  labels <- c()
    

  adj_mat <- adj_gen(0.02,p)
  while(FALSE){
    while(all(adj_mat ==0)){
      print('¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬')
      adj_mat <- adj_gen(0.02,p)
    }
  }
  new.omega.true <- prec_from_adj(adj_mat, p)
  new.sigma.true <- solve(new.omega.true)
  
  y = mvtnorm::rmvnorm(n, rep(0, p), new.sigma.true)
  if(scale.data){
    y = scale(y)
  }
  
  res.ssl <- GM(y, list_hyper, list_init)
  
  bool_up <- upper.tri(res.ssl$m_delta)
  
  labels <- c(labels, adj_mat[bool_up])
  predictions <- c(predictions, res.ssl$m_delta[bool_up])
    

  pred <- prediction(predictions, labels)
  perf <- performance(pred, measure = "tpr", x.measure = "fpr")
  
  auc_value <- performance(pred, "auc")@y.values[[1]]
  #aic_value <- performance(pred, "aic")@y.values[[1]]
  
  
  return(list(predictions = predictions, labels = labels, pred = pred, perf = perf, auc_value = auc_value, v0 = list_hyper$v0, v1 = list_hyper$v1))#, aic_value = aic_value))
}


plot_ssl_ROC <- function(res.ssl){
  auc_value = res.ssl$auc_value
  perf = res.ssl$perf
  plot(perf, main = '', xlim = c(0, 1), ylim = c(0, 1), lwd = 2, colorize = TRUE)
  title(main = paste('ROC plot of SSL method. AUC = ',auc_value), sub = paste('v0 = ',list_hyper$v0, 'v1 =', list_hyper$v1))
  abline(a= 0, b = 1, col = "gray80")
  #table(predictions,labels > 0.5)
}


get_ave_AUC <- function(list_hyper, list_init, count = 50){
  pb <- progress_bar$new(format = "[:bar] :percent ETA: :eta", total = count)
  AUC_vals <- c()
  for(p in 1:count){
    pb$tick()
    res.ssl <- generate_GM_preds_labels(list_hyper, list_init)
    AUC_vals <- c(AUC_vals, res.ssl$auc_value)
    
  }
  return(mean(AUC_vals))
}

v0_grid <- function(list_hyper, list_init, v0_list = list(0.05, 0.1, 0.5, 1, 5, 10, 50, 100), plot = F){
  op <- list()
  for(v0 in v0_list){
    list_hyper$v0 <- v0
    auc_value <- get_ave_AUC(list_hyper, list_init)
    cat(paste0('v0 = ',list_hyper$v0, '   AUC = ', auc_value,'\n'))
    if(plot==T){
      res.ssl <- generate_GM_preds_labels(list_hyper, list_init)
      plot_ssl_ROC(res.ssl)
    }
    op <- c(op, list(v0,v1,auc_value))
  }
  list_hyper$v0 <- 0.5
  
  return(op)
}