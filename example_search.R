if(!require(remotes)) install.packages("remotes")
remotes::install_github("XiaoyueXI/navigm")
library(navigm)

seed <- 123; set.seed(seed)
N <- 100; P <- 50;
Y <- matrix(rnorm(N*P), nrow = N, ncol = P)

# estimate a precision matrix based on Y 
lambda <- 2
v0_v <- s0_v <- seq(1e-2, 1, length.out = 16)
v1 <- s1 <- 100
P <- ncol(V)

#
list_hyper <- list(
  lambda = lambda,
  v0_v = v0_v,
  v1 = v1,
  s0_v = s0_v, 
  s1 = s1,
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

res1 <- navigm(Y, method = 'GM', 
                     version = 1,
                     list_hyper = list_hyper,
                     list_init = list_init, 
                     tol = 1e-3,
                     maxit = 1e5)

res1$estimates$m_delta