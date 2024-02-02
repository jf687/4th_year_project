library(navigm) # https://github.com/XiaoyueXI/navigm (our package which implements GM using VBECM)
set.seed(1234)

# see attached
load("data_seed1.rda")

#
lambda <- 2; v0 <- 0.1; v1 <- 100; P <- ncol(Y)
# not grid search for spike variance here; but in practice we should do so
list_hyper <- list(
  lambda = lambda,
  v0 = v0,
  v1 = v1,
  a_tau = 2,
  b_tau = 2,
  a_rho = 1,
  b_rho = P)
list_init <-
  list(
    alpha_tau = 1,
    beta_tau = 1,
    alpha_rho = 1,
    beta_rho = P-1
  )

#
out <-
  navigm_core(Y = net$Y,
              V = V,
              method = "GM",
              list_hyper = list_hyper,
              list_init = list_init, 
              tol = 1e-3,
              maxit = 1e5,
              debug = T,
              version = 1)

# an example of cutoffs 
library(ROCR)
hist(out$estimates$m_delta) 
bool_up <- upper.tri(net$A)
pred <- prediction(out$estimates$m_delta[bool_up], net$A[bool_up])
range(pred@cutoffs[[1]])
quantile(pred@cutoffs[[1]])

# 
table(net$A[bool_up],out$estimates$m_delta[bool_up] > 0.5)

