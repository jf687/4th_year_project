GGM_gen <- function(n,p){
  ggm.sf = huge::huge.generator(n=n, d=p,graph = 'scale-free') 
  omega.true = ggm.sf$omega # The true precision matrix
   # True sparsity
  omega.true[which(abs(omega.true)<1e-15)]=0 # IMPORTANT: ensure precision matrix has entries equal to exactly zero (huge gives some very small nonzero elements for some reason)
  return(list(ggm.sf = ggm.sf, omega.true = omega.true, sigma.true = ggm.sf$sigma, A = ggm.sf$theta, sparsity = ggm.sf$sparsity))
}


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

TPR = function(g, g.hat){
  p <- nrow(g[, ])
  g <- g[, ]!=0
  g.hat <- g.hat[, ]!=0
  diag(g) <- rep(0, p)
  diag(g.hat[, ]) <- rep(0, p)
  tp = sum(g.hat[, ] == 1 & g[, ] == 1)/2
  pos = sum(g[, ] == 1)/2
  return(tp/pos)
}

FPR = function(g, g.hat){
  p <- nrow(g[, ])
  g <- g[, ]!=0
  g.hat <- g.hat[, ]!=0
  diag(g) <- rep(0, p)
  diag(g.hat[, ]) <- rep(0, p)
  fp = sum(g.hat[, ] == 1 & g[, ] == 0)/2
  neg = (sum(g[, ] == 0)-p)/2
  return(fp/neg)
}

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


get_AUC = function(sim.obj, method='Glasso', cutoff=NULL){
  # Adding the point (0,0)
  if(method=='Glasso'){
    x = c(0,sim.obj$mean.FPR.glasso[-1])
    y = c(0,sim.obj$mean.TPR.glasso[-1])
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

adj_gen <- function(prob,p){
  adjacency_matrix <- matrix(0, nrow = p, ncol = p)
  for(i in 1:p){
    for(j in i:p-1){
      adjacency_matrix[i,j+1] <- rbinom(1,1,prob)
      adjacency_matrix[j+1,i] <- adjacency_matrix[i,j+1]
    }
  }
  return(adjacency_matrix)
}

prec_from_adj <- function(A,p){
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

# Function to calculate AIC. N 
AIC_GSS <- function(estimates, N){
  Omega <- estimates$Omega
  
  Omega[estimates$m_delta <= 0.5] <- 0 
  d <- sum(estimates$m_delta > 0.5)
  
  diag(Omega) <- diag(estimates$Omega) 
  
  det.Omega <- determinant(Omega, logarithm = TRUE)$modulus[1] 
  
  return( sum(diag(estimates$Sigma %*% Omega)) - N * det.Omega + 2 * d)
  }

# Function to calculate BIC
BIC_GSS <- function(estimates, N){
  Omega <- estimates$Omega
  
  Omega[estimates$m_delta <= 0.5] <- 0 
  d <- sum(estimates$m_delta > 0.5)
  
  diag(Omega) <- diag(estimates$Omega) 
  
  det.Omega <- determinant(Omega, logarithm = TRUE)$modulus[1] 
  return( sum(diag(estimates$Sigma %*% Omega)) - N * det.Omega + log(N)* d)
}

return_closest_threshold <- function(thresholds, target = 0.5) {
  op <- -1
  min.diff <- 1
  for(i in 1:length(thresholds)){
    if(abs(thresholds[[i]]-target) < min.diff){
      min.diff <- abs(thresholds[[i]]-target) 
      op <- i
    }
  }
  return(op)
}