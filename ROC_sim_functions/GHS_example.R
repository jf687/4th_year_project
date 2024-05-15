
source('~/4th_year_project/aux_method_funcs/GHS/GHS.R')
source('~/4th_year_project/aux_sim_funcs/ROC_calcs.R')
library(huge)
library(progress)

test_ghs = function(){
  set.seed(999)
  sf = huge::huge.generator(n=100, d=50,graph = 'scale-free', prob = 0.02)
  y = mvtnorm::rmvnorm(100, rep(0, 50), sf$sigma)
  y = scale(y)
  opGHS <- GHS_ROC_sim(sf,y)
  plot.GHS(opGHS)
}

# Example of how to use GHS (graphical horseshoe)
GHS_ROC_sim <- function(data.sf, y.sf, x = 0.25){
  if(x >0.5){
    x <- 0.25
  }
  theta.true = data.sf$omega # The precision matrix
  theta.true[which(theta.true<10e-5,arr.ind=T)]=0
  
  n <- dim(data.sf$data)[[1]]
  p <- dim(data.sf$data)[[2]]
  
  #data.sf$sparsity # True sparsity: 0.02
  start.time <- Sys.time()
  ghs.res = GHS(t(y.sf)%*%y.sf,n,burnin=100,nmc=1000)
  end.time <- Sys.time()
  theta.est.ghs = apply(ghs.res$thetas.sampled, c(1,2), mean) # Posterior mean of MCMC samples
  
  lambda = mean(ghs.res$lambdas.sampled)

  # To threshold, Li et al. (2017) use the symmetric central 50% posterior credible intervals for variable selection. 
  # That is, if the 50% posterior credible interval of an off-diagonal element of does not contain zero, that element is considered a discovery, and vice versa. 
  theta.ghs.credible.lower = apply(ghs.res$thetas.sampled, c(1,2), FUN = function(s) quantile(s, 0.25))
  theta.ghs.credible.upper = apply(ghs.res$thetas.sampled, c(1,2), FUN = function(s) quantile(s, 0.75))
  theta.ghs.selected = theta.est.ghs
  theta.ghs.selected[(theta.ghs.credible.lower < 0 & theta.ghs.credible.upper > 0)] = 0
  
  
  diag(theta.ghs.selected) <- 0
  diag(theta.true) <- 0

  
  sparsity = sparsity(theta.ghs.selected)
  precision = precision(theta.true!= 0, theta.ghs.selected!= 0)
  recall = recall(theta.true!= 0, theta.ghs.selected!= 0)
  
  # To vary the threshold for ROC curves etc, you can increase/decrease the width of the credible intervals used to determine nonzero effects. 
  # I.e., replace 0.25 and 0.75 by x and 1.0-x for any 0<x<0.5.
  quant.thresh = seq(0,0.48,length.out=50)
  theta.path.ghs = list()
  TPR.ghs = rep(0,length(quant.thresh))
  FPR.ghs = rep(0,length(quant.thresh))
  
  pb <- progress_bar$new(format = "[:bar] :percent ETA: :eta", total = length(quant.thresh))
  
  for(l in 1:length(quant.thresh)){
    
    pb$tick()
    theta.ghs.credible.lower = apply(ghs.res$thetas.sampled, c(1,2), FUN = function(s) quantile(s, quant.thresh[l]))
    theta.ghs.credible.upper = apply(ghs.res$thetas.sampled, c(1,2), FUN = function(s) quantile(s, 1-quant.thresh[l]))
    theta.ghs.selected = theta.est.ghs
    theta.ghs.selected[(theta.ghs.credible.lower < 0 & theta.ghs.credible.upper > 0)] = 0
    diag(theta.ghs.selected) = 1
    theta.path.ghs[[l]] = theta.ghs.selected
    TPR.ghs[l] = TPR(theta.true!=0, theta.ghs.selected!=0)
    FPR.ghs[l] = FPR(theta.true!=0, theta.ghs.selected!=0)
  }
  
  x = c(0,FPR.ghs)
  y = c(0,TPR.ghs)
  AUC=flux::auc(c(0, FPR.ghs), c(0,TPR.ghs))

  return(list(n = n, p = p, AUC = AUC, precision = precision, lambda = lambda,
              recall = recall, sparsity = sparsity, FPR.ghs = FPR.ghs, TPR.ghs = TPR.ghs,
              time = end.time - start.time))

}

GHS_average <- function(ghs.sf, N){
  
  res = list()
  seed <- 0
  
  n <- dim(ghs.sf$data)[[1]]
  p <- dim(ghs.sf$data)[[2]]
  
  AUC.list = c()
  REC.list = c()
  PREC.list = c()
  SPARS.list = c()
  lambda.list = c()
  time.list = c()
  
  for(i in 1:N){
    print(i)
     
    set.seed(seed)
    y = mvtnorm::rmvnorm(n, rep(0, p), ghs.sf$sigma)
    y = scale(y)
    
    
    op <- GHS_ROC_sim(ghs.sf, y)
    
    AUC.list = c(AUC.list, op$AUC)
    REC.list = c(REC.list, op$recall)
    PREC.list = c(PREC.list, op$precision)
    SPARS.list = c(SPARS.list, op$sparsity)
    lambda.list = c(lambda.list, op$lambda)
    time.list = c(time.list, op$time)
    seed <- seed + 1
  }
  
  res$AUC = list(mean = mean(AUC.list), SE = sd(AUC.list)/sqrt(N))
  res$REC = list(mean = mean(REC.list), SE = sd(REC.list)/sqrt(N))
  res$PREC = list(mean = mean(PREC.list), SE = sd(PREC.list)/sqrt(N))
  res$SPARS = list(mean = mean(SPARS.list), SE = sd(SPARS.list)/sqrt(N))
  res$lambda = list(mean = mean(lambda.list), SE = sd(lambda.list)/sqrt(N))
  res$time = list(mean = mean(time.list), SE = sd(time.list)/sqrt(N))
  
  res$n = n
  res$p = p
  
  return(res)
  
}

plot.GHS <- function(op){
  plot(c(op$FPR.ghs,1),c(op$TPR.ghs,1),type='l', col='red', 
       xlab = "False Positive Rate", ylab = "True Positive Rate", xlim = c(0,1) , 
       main = "ROC plot using GHS method with data size n = 200, p = 100.")
  abline(a= 0, b = 1, col = "gray80")
  i <- gmean_maximisaton(op$FPR.ghs, op$TPR.ghs)$i
  points(op$FPR.ghs[[i]], y = op$TPR.ghs[[i]])
}

GHS =  function(S,n,burnin=500,nmc=3000){
  # GHS MCMC sampler using data-augmented
  # block (column-wise) Gibbs sampler
  # Input:
  #     S = Y'*Y : sample covariance matrix * n
  #     n: sample size
  #     burnin, nmc : number of MCMC burnins and saved samples
  
  # Output:
  #     thetas.sampled: p by p by nmc arrays of saved posterior samples of
  #     precision matrix
  #     lambdas.sampled: p*(p-1)/2 by nmc vector of saved samples of lambda
  #     squared (local tuning parameter)
  #     taus.sampled: 1 by nmc vector of saved samples of tau squared (global
  #     tuning parameter)
  p = nrow(S)
  omega_save = array(0,c(p,p,nmc))
  lambda_sq_save = matrix(0,p*(p-1)/2,nmc)
  tau = rep(0,nmc)
  tau_sq_save = rep(0,nmc)
  
  ind_all = matrix(0,p-1,p)
  for (i in 1:p){
    if (i==1) ind = matrix(2:p,ncol=1)
    else if (i==p) ind = matrix(1:(p-1),ncol=1)
    else ind = matrix(c(1:(i-1),(i+1):p),ncol=1)
    
    ind_all[,i] = ind
  }
  
  # set initial values
  Omega = diag(1,p) 
  Sigma = diag(1,p)
  Nu = matrix(1,p,p)
  Lambda_sq = matrix(1,p,p)
  Nu[1:p,1:p] = 1
  tau_sq = 1
  xi = 1
  
  for (iter in 1:(burnin+nmc)){  
    ### sample Sigma and Omega=inv(Sigma)
    for (i in 1:p){
      
      ind = ind_all[,i]     
      Sigma_11 = Sigma[ind,ind]
      sigma_12 = Sigma[ind,i]
      sigma_22 = Sigma[i,i]
      s_21 = S[ind,i]
      s_22 = S[i,i]
      lambda_sq_12 = Lambda_sq[ind,i]
      nu_12 = Nu[ind,i]
      ## sample gamma and beta
      gamma = rgamma(1,shape=(n/2+1),scale=2/s_22)    # random gamma with shape=n/2+1, rate=s_22/2. Sample 1. 
      #gamma= 0.06757892
      inv_Omega_11 = Sigma_11 - sigma_12%*%t(sigma_12)/sigma_22
      inv_C = s_22*inv_Omega_11+diag(1/(lambda_sq_12*tau_sq))
      #if(i!=1) return(inv_C)
      inv_C_chol = chol(inv_C) 
      #cat(iter,' ', i, '\n')
      mu_i= tryCatch({
        -solve(inv_C,s_21)
      }, error = function(e) {
        return(list(inv_C,s_21))
      })
      #mu_i = -solve(inv_C,s_21)
      beta = mu_i+ solve(inv_C_chol,rnorm(p-1)) # This is where it fails for too high dim data
      #kkk=c(0.1156077, 1.7215140, 0.10147592)
      #beta = mu_i+ solve(inv_C_chol,kkk)
      omega_12 = beta
      omega_22 = gamma + t(beta)%*%inv_Omega_11%*%beta;
      ## sample lambda_sq and nu
      rate = omega_12^2/(2*tau_sq)+1/nu_12
      lambda_sq_12 = 1/rgamma(length(omega_12), shape=1, scale = 1/rate)   # random inv gamma with shape=1, rate=rate. 
      nu_12 = 1/rgamma(length(lambda_sq_12), shape=1, scale = 1/(1+1/lambda_sq_12))   # random inv gamma with shape=1, rate=1+1/lambda_sq_12
      #lambda_sq_12 = 1/c(4.059683, 1.923486, 1.999748)
      #nu_12 = 1/c(1.5817724, 0.4459362, 0.9780721)
      ## update Omega, Sigma, Lambda_sq, Nu
      Omega[i,ind] = omega_12
      Omega[ind,i] = omega_12
      Omega[i,i] = omega_22
      temp = inv_Omega_11%*%beta
      Sigma_11 = inv_Omega_11 + temp%*%t(temp)/gamma;
      sigma_12 = -temp/gamma
      sigma_22 = 1/gamma
      Sigma[ind,ind] = Sigma_11
      Sigma[i,i] = sigma_22
      Sigma[i,ind] = sigma_12
      Sigma[ind,i] = sigma_12
      Lambda_sq[i,ind] = lambda_sq_12
      Lambda_sq[ind,i] = lambda_sq_12
      Nu[i,ind] = nu_12
      Nu[ind,i] = nu_12
    }
    # sample tau_sq and xi
    omega_vector =  Omega[lower.tri(Omega)] # Only elements in Omega below the diagonal
    lambda_sq_vector = Lambda_sq[lower.tri(Lambda_sq)]
    rate = 1/xi + sum(omega_vector^2/(2*lambda_sq_vector))
    tau_sq = 1/rgamma(1, shape=(p*(p-1)/2+1)/2, scale = 1/rate) # inv gamma w/ shape=(p*(p-1)/2+1)/2, rate=rate. Only sample one. 
    xi = 1/rgamma(1, shape=1, scale = 1/(1+1/tau_sq))   # inv gamma w/ shape=1, rate=1+1/tau_sq
    
    # save Omega, lambda_sq, tau_sq
    if (iter > burnin){         
      omega_save[,,(iter-burnin)] = Omega
      lambda_sq_save[,(iter-burnin)] = lambda_sq_vector
      tau_sq_save[(iter-burnin)] = tau_sq
    }
  }
  return(list(thetas.sampled = omega_save,lambdas.sampled = lambda_sq_save,
              taus.samples = tau_sq_save))
}


