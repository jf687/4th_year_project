source('~/4th_year_project/Sims/ROC_calcs.R')

source('~/4th_year_project/GM/GM.R')
source('~/4th_year_project/GM/fun_GM.R')
source('~/4th_year_project/GM/fun_utils.R')



generate_GM_preds_labels <- function(list_hyper, list_init,n= 200, p = 100, scale.data = T){

  ggm.sf = huge::huge.generator(n=n, d=p,graph = 'scale-free', prob = 0.02)
    

  adj_mat <- ggm.sf$theta
  Omega.true <- ggm.sf$omega
  Sigma.true <- ggm.sf$sigma
  
  y = mvtnorm::rmvnorm(n, rep(0, p), Sigma.true)
  if(scale.data){
    y = scale(y)
  }
  
  res.ssl <- GM(y, list_hyper, list_init)
  
  bool_up <- upper.tri(res.ssl$m_delta)
  
  labels <- adj_mat[bool_up]
  predictions <- res.ssl$m_delta[bool_up]
    
  pred <- ROCR::prediction(predictions, labels)
  perf <- ROCR::performance(pred, measure = "tpr", x.measure = "fpr")
  auc_value <- ROCR::performance(pred, "auc")@y.values[[1]]
  thresh <- return_closest_threshold(pred@cutoffs[[1]])
  
  prec_value <- ROCR::performance(pred, "prec")@y.values[[1]][thresh]
  rec_value <- ROCR::performance(pred, "rec")@y.values[[1]][thresh]
  
  AIC <- aic_ssl(y, res.ssl$Omega, n,p)
  BIC <- bic_ssl(y, res.ssl$Omega, n, p)
  return(list(n = n, p = p , predictions = predictions, labels = labels, pred = pred, perf = perf, 
              auc_value = auc_value, prec_value = prec_value, rec_value = rec_value, 
                v0 = list_hyper$v0, v1 = list_hyper$v1, AIC = AIC, BIC = BIC, sparsity = ggm.sf$sparsity))
}

aic_ssl <- function(Y, omega.est, n,p){
  S = cov(Y)
  if (mean(dim(S) == dim(omega.est)) != 1) 
    stop("matrices must have the same dimension")
  if (det(omega.est) <= 0) 
    stop("precision matrix must be positive definite.")
  if (!isSymmetric(S)) 
    stop("sample covariance matrix must be symmetric")
  if (n <= 0) 
    stop("number of observations n must be positive")
  p <- nrow(omega.est)
  omega.est2 <- omega.est
  diag(omega.est2) <- rep(0, p)
  d <- sum(omega.est != 0)/2
  loglik <- n * log(det(omega.est))/2- n * sum(diag(S%*% omega.est))/2
  aic = -2 * loglik + 2 * d
  return(aic)
}

bic_ssl <- function(Y, omega.est, n,p){
  S = cov(Y)
  if (mean(dim(S) == dim(omega.est)) != 1) 
    stop("matrices must have the same dimension")
  if (det(omega.est) <= 0) 
    stop("precision matrix must be positive definite.")
  if (!isSymmetric(S)) 
    stop("sample covariance matrix must be symmetric")
  if (n <= 0) 
    stop("number of observations n must be positive")
  p <- nrow(omega.est)
  omega.est2 <- omega.est
  diag(omega.est2) <- rep(0, p)
  d <- sum(omega.est != 0)/2
  loglik <- n * log(det(omega.est))/2- n * sum(diag(S%*% omega.est))/2
  bic = -2 * loglik + d*log(n)
  return(bic)
}

plot_ssl_ROC <- function(res.ssl){
  auc_value = res.ssl$auc_value
  perf = res.ssl$perf
  plot(perf, main = '', xlim = c(0, 1), ylim = c(0, 1), lwd = 2, colorize = TRUE)
  title(main = paste('ROC plot of SSL method. AUC = ',auc_value), sub = paste('v0 = ',list_hyper$v0, 'v1 =', list_hyper$v1))
  abline(a= 0, b = 1, col = "gray80")
  #table(predictions,labels > 0.5)
}


get_ave_AUC <- function(list_hyper, list_init, count = 50, n_= 200, p_ = 100, v.length = 1){
  pb <- progress_bar$new(format = "[:bar] :percent ETA: :eta", total = count*v.length)
  AUC_vals <- c()
  PREC_vals <- c()
  REC_vals <- c()
  AIC_vals <- c()
  BIC_vals <- c()
  SPARS_vals <- c()
  
  for(p in 1:count){
    pb$tick()
    res.ssl <- generate_GM_preds_labels(list_hyper, list_init, n = n_ , p = p_)
    AUC_vals <- c(AUC_vals, res.ssl$auc_value)
    PREC_vals <- c(PREC_vals, res.ssl$prec_value)
    REC_vals <- c(REC_vals, res.ssl$rec_value)
    AIC_vals <- c(AIC_vals, res.ssl$AIC)
    BIC_vals <- c(BIC_vals, res.ssl$BIC)
    SPARS_vals <- c(SPARS_vals, res.ssl$sparsity)
    
  }
  return(list(AUC = mean(AUC_vals), PREC = mean(PREC_vals), REC = mean(REC_vals), AIC = mean(AIC_vals), BIC = mean(BIC_vals), sparsity = mean(SPARS_vals)))
}

v0_grid <- function(list_hyper, list_init, v0_list = list(0.05, 0.1, 0.5, 1, 5, 10, 50, 100), plot = F, save.data = F,count = 30, n = 200, p = 100){
  
  v.length = length(v0_list)
  v0.list = list()
  auc.list = list()
  prec.list = list()
  rec.list = list()
  aic.list = list()
  bic.list = list()
  sparsity.list = list()
  
  for(v0 in v0_list){
    list_hyper$v0 <- v0
    out <- get_ave_AUC(list_hyper, list_init,count, n,p,v.length = v.length)
    auc_value <- out$AUC
    prec_value <- out$PREC
    rec_value <- out$REC
    aic_value <- out$AIC
    bic_value <- out$BIC
    spars_value <- out$sparsity
    if(plot==T){
      res.ssl <- generate_GM_preds_labels(list_hyper, list_init)
      plot_ssl_ROC(res.ssl)
    }
    
    v0.list = c(v0.list, v0)
    auc.list = c(auc.list, auc_value)
    prec.list = c(prec.list, prec_value)
    rec.list = c(rec.list, rec_value)
    aic.list = c(aic.list, aic_value)
    bic.list = c(bic.list, bic_value)
    sparsity.list = c(sparsity.list, spars_value)
    
  }
  list_hyper$v0 <- 0.5
  
  op <- list(n = n, p = p , v0 =  v0.list, auc = auc.list, prec = prec.list, rec = rec.list, aic = aic.list, bic = bic.list, sparsity = sparsity.list)
  if(save.data){
  save(op, file = 'SSL_run.rda')
  }
  return(op)
}

GM_ROC_sim <- function(n = 400, p = 50, plot = F){
  
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
  v0_v <- seq(1e-6, 9e-2, length.out = 32)
  v1 <- 100
  
  #initialise list parameters
  list_hyper <- list(
    lambda = 2,
    v0_v = v0_v,
    v1 = v1,
    a_tau = 2,
    b_tau = 2,
    a_rho = 1,
    b_rho = p/3)
  
  
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
  
  ## THRESHOLD ON 0.5
  
  prec_value <- ROCR::performance(pred, "prec")@y.values[[1]][613]

  rec_value <- ROCR::performance(pred, "rec")@y.values[[1]][613]
  
  if(plot){
    plot(v0_v[2:32],res1$model_criterion$value[2:32], xlab = 'value of spike variance v0', ylab = 'AIC', type = 'l')
    plot(perf, main = '', xlim = c(0, 1), ylim = c(0, 1), lwd = 2, colorize = TRUE)
  }
  # an example of cutoffs 
  hist(res1$estimates$m_delta) 
  browser()
  cm <- confusion.matrix(labels[bool_up], predictions[bool_up] > 0.5)
  return(list(n=n,p=p,auc_value = auc_value, prec_value = prec_value, rec_value = rec_value, AIC = res1$model_criterion$value, 
              table = table, predictions = predictions, labels = labels))
  
  
}