source('libraries.R')

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
  rec_value <- ROCR::performance(pred, "recall")@y.values[[1]][613]
  
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