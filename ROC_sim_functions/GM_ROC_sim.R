source('~/4th_year_project/aux_sim_funcs/ROC_calcs.R')
source('~/4th_year_project/aux_sim_funcs/select.v0.bic.R')

source('~/4th_year_project/aux_method_funcs/GM/GM.R')
source('~/4th_year_project/aux_method_funcs/GM/fun_GM.R')
source('~/4th_year_project/aux_method_funcs/GM/fun_utils.R')

if(F){
  ggm.sf = huge::huge.generator(n=n, d=p,graph = 'scale-free', prob = 0.02)
  y = mvtnorm::rmvnorm(n, rep(0, p), ggm.sf$sigma)
  y = scale(y)
}


GM_ROC_sim <- function(ggm.sf,y, list_hyper, list_init, v0.vals, t.vals){
  n <- dim(ggm.sf$data)[[1]]
  p <- dim(ggm.sf$data)[[2]]
  
  
  ### v0 and t selection here
  
  v0.selection <- select.v0.bic(y, list_hyper, list_init, v0.vals,t.vals,n,p)
  list_hyper$v0 <- v0.selection$v0.opt
  thresh <- v0.selection$t.opt
  

  adj_mat <- ggm.sf$theta
  Omega.true <- ggm.sf$omega
  Sigma.true <- ggm.sf$sigma

  
  res.ssl <- GM(y, list_hyper, list_init)
  
  bool_up <- upper.tri(res.ssl$m_delta)
  
  labels <- adj_mat[bool_up]
  predictions <- res.ssl$m_delta[bool_up]
    
  pred <- ROCR::prediction(predictions, labels)
  perf <- ROCR::performance(pred, measure = "tpr", x.measure = "fpr")
  auc_value <- ROCR::performance(pred, "auc")@y.values[[1]]

  prec_value <- precision( Omega.true != 0, res.ssl$m_delta > thresh)
  rec_value <- recall(Omega.true != 0, res.ssl$m_delta > thresh)
  
  AIC <- aic_ssl(y, res.ssl$Omega, n,p)
  BIC <- bic_ssl(y, res.ssl$Omega, n, p)
  
  sparsity <- sparsity(res.ssl$m_delta > thresh)
  
  return(list(n = n, p = p , predictions = predictions, labels = labels, pred = pred, perf = perf, 
              auc_value = auc_value, prec_value = prec_value, rec_value = rec_value, 
                v0 = list_hyper$v0, v1 = list_hyper$v1, AIC = AIC, BIC = BIC, sparsity = sparsity))
}

GM_average <- function(ggm.sf, N, list_hyper, list_init, v0.vals, t.vals){
  
  seed <- 0
  res = list()
  
  n <- dim(ggm.sf$data)[[1]]
  p <- dim(ggm.sf$data)[[2]]
  list_hyper$br <- p
  
  AUC.list = c()
  REC.list = c()
  PREC.list = c()
  SPARS.list = c()
  v0.list = c()
  
  for(i in 1:N){
    
    set.seed(seed)
    y = mvtnorm::rmvnorm(n, rep(0, p), ggm.sf$sigma)
    y = scale(y)
    
    op <- GM_ROC_sim(ggm.sf, y,  list_hyper, list_init, v0.vals, t.vals )
    
    AUC.list = c(AUC.list, op$auc_value)
    REC.list = c(REC.list, op$rec_value)
    PREC.list = c(PREC.list, op$prec_value)
    SPARS.list = c(SPARS.list, op$sparsity)
    v0.list = c(v0.list, op$v0)
    
    seed <= seed + 1
    
  }

  res$AUC = list(mean = mean(AUC.list), SE = sd(AUC.list)/sqrt(N))
  res$REC = list(mean = mean(REC.list), SE = sd(REC.list)/sqrt(N))
  res$PREC = list(mean = mean(PREC.list), SE = sd(PREC.list)/sqrt(N))
  res$SPARS = list(mean = mean(SPARS.list), SE = sd(SPARS.list)/sqrt(N))
  res$AUC = list(mean = mean(AUC.list), SE = sd(AUC.list)/sqrt(N))
  res$v0 = list(mean = mean(v0.list), SE = sd(v0.list)/sqrt(N))
  
  res$n = n
  res$p = p
  res$v1 = list_hyper$v1
  
  return(res)
  
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



