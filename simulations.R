library(huge)
library(MASS)
library(matrixcalc)
library(pROC)
library(PRROC)
library(ROCR)

set.seed(123)

adj_gen <- function(prob,n,p){
  adjacency_matrix <- matrix(rbinom(n*p, size = 1, prob = prob), nrow = n, ncol = p)
  return(adjacency_matrix)
}

prec_from_adj <- function(A){
  E <- A
  nb_edges <- sum(E == 1)
  vec_magnitude <- c(-1,1)
  
  E[A == 1] <- runif(nb_edges, min = vec_magnitude[1], max = vec_magnitude[2])
  
  E_bar <- (E + t(E)) / 2
  
  msign <- matrix(1, nrow = nrow(E), ncol = ncol(E))
  msign[upper.tri(msign)] <- sample(c(-1,1), size = sum(upper.tri(msign)),  prob = c(0.5, 0.5), replace = TRUE)
  msign[lower.tri(msign)] <- t(msign)[lower.tri(msign)]
  E_bar <- E_bar * msign
  
  # minimum eigenvalue
  min_eigen <- min(eigen(E_bar, only.values = TRUE)$values)
  if (min_eigen < 0) {
    Omega <- E_bar + (0.1 - min_eigen) * diag(p)
  } else{
    Omega <- E_bar + 0.1 * diag(p)
  }
  
  return(Omega)
  
}

# function that generates input data from a given covariance matrix
data_gen <- function(n,p, sigma){
  Y <- matrix(0,n,p)
  for (i in 1:n){
    y <- mvrnorm(p,mu = numeric(p),Sigma = sigma )
    Y[i,] <- y
  }
  return(Y)
}


#############################################
confusion.matrix = function (g, g.hat) {
  # A matrix containing the number of true positives (TP), false positives (FP), true negatives (TN) and false negatives (FN).
  if (mean(dim(g[, ]) == dim(g.hat[, ])) != 1) 
    stop("matrices must have the same dimension")
  if (mean((g[, ] + 0) %in% c(0, 1)) != 1 | mean((g.hat[, ] + 0) %in% c(0, 1)) != 1) 
    stop("g and g.hat must be adjacency matrices with elements in {0,1}")
  p = nrow(g[, ])
  g = g[, ]
  g.hat = g.hat[, ]
  diag(g) = rep(0, p)
  diag(g.hat[, ]) = rep(0, p)
  tp = sum(g.hat[, ] == 1 & g[, ] == 1)/10
  fp = sum(g.hat[, ] == 1 & g[, ] == 0)/10
  tn = (sum(g.hat[, ] == 0 & g[, ] == 0) - p)/10
  fn = sum(g.hat[, ] == 0 & g[, ] == 1)/10
  return(matrix(10 * c(tp, fp, fn, tn), nrow = 2, byrow = T)/2)
}

precision = function (g, g.hat) {
  # The proportion of predicted edges that are true
  conf.mat = confusion.matrix(g, g.hat)
  if (conf.mat[1, 1] == 0 & conf.mat[2, 1] == 0) {
    return(1)
  }
  else if (conf.mat[1, 1] == 0 & conf.mat[1, 2] == 0) {
    return(1)
  }
  else {
    return(conf.mat[1, 1]/(conf.mat[1, 1] + conf.mat[1, 2]))
  }
}

recall = function (g, g.hat) {
  # The proportion of true edges that were identified by the estimated graph
  conf.mat = confusion.matrix(g, g.hat)
  if (conf.mat[1, 1] == 0 & conf.mat[2, 1] == 0) {
    return(1)
  }
  else {
    return(conf.mat[1, 1]/(conf.mat[1, 1] + conf.mat[2, 1]))
  }
}

##########################################

#function to calculate the sparsity of a matrix (note that when huge generator is used, use $sparsity)
sparsity <- function(omega, strict = TRUE, threshold = 1e-15){
  # Calculate the total number of elements in the matrix
  total_elements <- length(omega)
  
  # Calculate the number of zero elements in the matrix, or the number below a given threshold
  if (strict){
    zero_elements <- sum(omega == 0)
  } else {
    zero_elements <- sum(omega < threshold)
  }
  
  # Calculate the sparsity as the proportion of zero elements
  sparsity <- 1 - zero_elements / total_elements
  
  return(sparsity)
}

#function to return a treated GGM example from HUGE
GGM_gen <- function(n,p){
  ggm.sf = huge::huge.generator(n=n, d=p,graph = 'scale-free') 
  omega.true = ggm.sf$omega # The true precision matrix
  #ggm.sf$sparsity # True sparsity
  omega.true[which(abs(omega.true)<1e-15)]=0 # IMPORTANT: ensure precision matrix has entries equal to exactly zero (huge gives some very small nonzero elements for some reason)
  return(list(ggm.sf = ggm.sf, omega.true = omega.true, sigma.true = ggm.sf$sigma, A = ggm.sf$theta))
}

# function that generates the inverse precision matrix based on the GM priors
GM_gen <-function(n,p,list_hyper, list_init, thresh = 0.02){
  omega <- matrix(0,nrow = p, ncol = p)
  
  lambda <- list_hyper$lambda
  
  a <- list_hyper$a
  b <- list_hyper$b
  
  ar <- list_hyper$ar
  br <-list_hyper$br
  
  v0 <- list_hyper$v0
  v1 <- list_hyper$v1
  
  tau <- rgamma(1,a,b)
  rho <- thresh
  
  var1 <- v1^2/tau
  var0 <- v0^2/tau
  
  for (i in 1:p){
    omega[i,i] <- rexp(1,lambda/2)
    for (j in 1:i-1){
      if (rbinom(1, size =1, prob =  rho) == 1){
        omega[i,j] <- rnorm(1,0,var1)
        omega[j,i] <- omega[i,j]
      } else {
        omega[i,j] <- rnorm(1,0,var0)
        omega[j,i] <- omega[i,j]
      }
    }
  }
  
  min_eigen <- min(eigen(omega, only.values = TRUE)$values)
  if (min_eigen < 0) {
    omega <- omega + (0.1 - min_eigen) * diag(p)
  } else{
    omega <- omega + 0.1 * diag(p)
  }
  
  
  sigma <- solve(omega)
  
  return(list(Sigma = sigma, Omega = omega))
}

# function that generates the inverse precision matrix based on the EXTENDED GM priors


# example format for list_hyper_x:
#' list_hyper_x <- list(lambda = 2,
#'                    v0 = 100,
#'                    v1 = 0.5,
#'                    a = 2,
#'                    b = 2, 
#'                    ar = 1,
#'                    br = p, 
#'                    asig = 1,
#'                    bsig = 1,
#'                    n0 = -3,
#'                    t0 = 1)
#'            

GMx_gen <- function(n,p,list_hyper_x, X = list(1,2,3,4)){
  
  if (!is.numeric(n) || !is.numeric(p)) {
    stop("Input arguments must be numeric.")
  }
  
  a <- list_hyper_x$a
  b <- list_hyper_x$br
  
  ar <- list_hyper_x$ar
  br <-list_hyper_x$br
  
  v0 <- list_hyper_x$v0
  v1 <- list_hyper_x$v1
  
  tau <- rgamma(1,a,b)
  
  var1 <- v1^2/tau
  var0 <- v0^2/tau
  
  asig <- list_hyper_x$asig
  bsig <- list_hyper_x$bsig
  
  n0 <- list_hyper_x$n0
  t0 <- list_hyper_x$t0
  
  zeta <- rnorm(1,mean = n0, sd = t0)
  cat('$$$',zeta,'$$$')
  
  sig_neg2 <- rgamma(1, asig, bsig)
  sig <- sig_neg2^-0.5
  
  BETA <- matrix(0, p, p)
  for (i in 1:p){
    for (j in 1:p){
      BETA[i,j] <- rnorm(1,0,sig)
    }
  }
  
  
  RHO <- array(0 , dim = c(p,p,length(X)))
  for (k in 1:length(X)){
    a.ord <- X[[k]]
    for (i in 1:p){
      for (j in 1:i){
        RHO[i,j,k] <- pnorm(zeta + a.ord*BETA[i,j])
        if (i!=j){
          RHO[j,i,k] <- RHO[i,j,k]  
        }
      }
    }
  }
  
  DELTA <- array(0 , dim = c(p,p,length(X)))
  OMEGA <- array(0 , dim = c(p,p,length(X)))
  
  for (k in 1:length(X)){
    for (i in 1:p){
      for (j in 1:i){
        DELTA[i,j,k] <- rbinom(1, size =1, prob = RHO[i,j,k]) 
      }
    }
    for (i in 1:p){
      if (i ==j){
        OMEGA[i,i,k] <- rexp(1,lambda/2)
      }
      else{
        for (j in 1:i){
          if (DELTA[i,j,k] == 1){
            OMEGA[i,j,k] <- rnorm(1,0,var1)
            OMEGA[j,i,k] <- OMEGA[i,j,k]
          } else {
            OMEGA[i,j,k] <- rnorm(1,0,var0)
            OMEGA[j,i,k] <- OMEGA[i,j,k]
          }
        }
      }
    }
    
  }
  
  for (k in 1:length(X)){
    
    
    omega_bar <- (OMEGA[,,k] + t(OMEGA[,,k])) / 2
    
    min_eigen <- min(eigen(omega_bar, only.values = TRUE)$values)
    if (min_eigen < 0) {
      OMEGA[,,k] <- omega_bar + (0.1 - min_eigen) * diag(p)
    } else{
      OMEGA[,,k] <- omega_bar + 0.1 * diag(p)
    }
    
  }
  
  SIGMA = array(0 , dim = c(p,p,length(X)))
  for (k in length(X)){
    SIGMA[,,k] <- solve(OMEGA[,,k])
  }
  
  return(list(Sigma = SIGMA, Omega = OMEGA))
}

# function that creates an example covariate vector
#create_X <- function(a_max,n) {
#  random_X <- sample(a_max, size= n, replace = TRUE)
#  return(random_X)
#}


glasso_sim <- function(n = 100,p=200, omega.true){
  
  # Input Validation
  if (!is.numeric(n) || !is.numeric(p)) {
    stop("Input arguments must be numeric.")
  }
  
  # Simulating input data
  Y <- mvtnorm::rmvnorm(n, rep(0, p), solve(omega.true))
  
  #running graphical Lasso
  result <- huge(Y, method = 'glasso')
  result= huge.select(result,criterion = 'stars',stars.thresh = 0.05)
  omega.est.glasso <- result$opt.icov
  omega.est.glasso[which(abs(omega.est.glasso)<1e-15)]=0
  omega.true <- matrix(omega.true, p , p)
  #plotting results  
  plot(result, layout = 'spring')
  #browser()
  #returning the precision and recall values
  prec <- precision(omega.true!=0, omega.est.glasso!=0)
  rec <- recall(omega.true!=0, omega.est.glasso!=0)
  
  return(result) # = result, precision = prec , recall = rec))
}

GM_sim <- function(n = 100,p=200, omega.true){
  
  # Input Validation
  if (!is.numeric(n) || !is.numeric(p)) {
    stop("Input arguments must be numeric.")
  }
  Y <- mvtnorm::rmvnorm(n, rep(0, p), solve(omega.true))
  omega.est.gm <- GM(Y, list_hyper, list_init)
  
  plot(result, layout = 'spring')
  
  # returning the precision and recall values
  # prec <- precision(omega.true!=0, omega.est.glasso!=0)
  # rec <- recall(omega.true!=0, omega.est.glasso!=0)
  
  
}



file_example <- function(n,p,X,list_hyper,list_hyper_x){
  if (FALSE){
    print('=== Example 1: GGM generation using HUGE ===')
    ggm.sf <- GGM_gen(n,p)
    omega.true <- ggm.sf$omega.true
    
    glasso_sim <- glasso_sim(n,p,omega.true)
    result <- glasso_sim$result
    print('=== Precision ===')
    print(glasso_sim$precision)
    print('=== Recall ===')
    print(glasso_sim$recall)
    cat('Sparsity = ',sparsity(omega.true),'\n')
    
    print('=== Example 2: GGM generation using GM model ===')
    gm <- GM_gen(n,p,list_hyper, list_init)
    omega.true <- gm$Omega
    
    glasso_sim <- glasso_sim(n,p,omega.true)
    result <- glasso_sim$result
    print('=== Precision ===')
    print(glasso_sim$precision)
    print('=== Recall ===')
    print(glasso_sim$recall)
    
    print('=== Example 3: GGM generation using GMx model ===')
    for (k in 1:length(X)){
      cat("***Covariate = ",k,"***")
      gm <- GMx_gen(n,p,list_hyper_x, list_init)
      omega.true <- gm$Omega[,,k]
      
      
      glasso_sim <- glasso_sim(n,p,omega.true)
      result <- glasso_sim$result
      
      print('=== Precision ===')
      print(glasso_sim$precision)
      print('=== Recall ===')
      print(glasso_sim$recall)
    }
  }
  print("### Example 4: GGM generation using 5% sparsity ###")
  A <- adj_gen(0.05,n,p)
  omega.true <- prec_from_adj(A)
  glasso_sim <- glasso_sim(n,p,omega.true)
  result <- glasso_sim$result
  print('=== Precision ===')
  print(glasso_sim$precision)
  print('=== Recall ===')
  print(glasso_sim$recall)
  
  
}

n <- 50
p <- 100
X <- list(1,2,3,4)

list_hyper_x <- list(lambda = 2,
                     v0 = 100,
                     v1 = 0.5,
                     a = 2,
                     b = 2, 
                     ar = 1,
                     br = p, 
                     asig = 1,
                     bsig = 1,
                     n0 = -2,
                     t0 = 1)

list_hyper <- list(lambda = 2,
                   v0 = 100,
                   v1 = 0.5,
                   a = 2,
                   b = 2, 
                   ar = 1,
                   br = p)

binarise_mat <- function(mat, t){
  for (index in 1:length(mat)){
    if (abs(mat[index]) > t){
      mat[index] <- 1
    }else{
      mat[index] <- 0
    }
  }
  return(mat)
}






file_run_1 <- function(n,p,X,list_hyper,N, seed=123){
  set.seed(seed)
  preds <- list()
  labels <- list()
  
  for (iteration in 1:N){
    
    ggm.sf <- GGM_gen(n,p)
    omega.true <- ggm.sf$omega.true
    result <- glasso_sim(n,p,omega.true)
    preds[[iteration]] <- abs(result$opt.icov[upper.tri(result$opt.icov)])
    labels[[iteration]] <- omega.true[upper.tri(omega.true)] != 0
    
  }
  
  # sparsity check
  sum(preds[[iteration]] > 1e-5) / length(preds[[iteration]])
  
  # plot_pr(labels = labels, predictions = preds)
  plot_pr(labels = labels[[iteration]], predictions = preds[[iteration]])
  
  # check plot_pr
  pred <- ROCR::prediction(preds, labels)
  perf <- performance(pred, "tpr", "fpr")
  plot(perf)
}

file_run_2 <- function(n,p,X,list_hyper,N, seed=123){
  set.seed(seed)
  preds <- list()
  labels <- list()
  
  for (iteration in 1:N){
    
    ggm.sf <- GM_gen(n, p, list_hyper, list_init)
    omega.true <- ggm.sf$Omega
    result <- glasso_sim(n,p,omega.true)
    preds[[iteration]] <- abs(result$opt.icov[upper.tri(result$opt.icov)])
    labels[[iteration]] <- omega.true[upper.tri(omega.true)] != 0
    
  }
  
  # sparsity check
  sum(preds[[iteration]] > 1e-5) / length(preds[[iteration]])
  
  # plot_pr(labels = labels, predictions = preds)
  plot_pr(labels = labels[[iteration]], predictions = preds[[iteration]])
  
  # check plot_pr
  pred <- ROCR::prediction(preds, labels)
  perf <- performance(pred, "tpr", "fpr")
  plot(perf)
}

file_run_3 <- function(){
  ROC_data <- form_ROC(nvec = 1, p,n, list_hyper, list_init, plot = FALSE, t = 0.5, x_lim = 0.02, y_lim = 0, title = 'hello', main_ = 'hello')
  pred <- ROCR::prediction(ROC_data$predictions, ROC_data$labels)
  perf <- performance(pred, "tpr", "fpr")
  plot(perf)
}

#file_run_1(n,p,X,list_hyper,N, seed=123)
