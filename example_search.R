
source('libraries.R')
seed <- 123; set.seed(seed)
N <- 100; P <- 50;
Y <- matrix(rnorm(N*P), nrow = N, ncol = P)

# estimate a precision matrix based on Y 
lambda <- 2
#<<<<<<< HEAD
v0_v <- s0_v <- seq(1e-2, 1, length.out = 16)
v1 <- s1 <- 100
P <- ncol(Y)# NOT RUN {
#generate data
L = huge.generator(d = 20, graph="hub")
out.mb = huge(L$data)
out.ct = huge(L$data, method = "ct")
out.glasso = huge(L$data, method = "glasso")

#model selection using ric
out.select = huge.select(out.mb)
plot(out.select)

#model selection using stars
out.select = huge.select(out.ct, criterion = "stars", stars.thresh = 0.05,rep.num=10)
plot(out.select)

#model selection using ebic
out.select = huge.select(out.glasso,criterion = "ebic")
plot(out.select)
#=======
v0_v <- seq(1e-2, 1, length.out = 16)
v1 <- 1000
P <- ncol(Y)
#>>>>>>> ab795f0b846821b23a836809263fc49ad232f80b
0
#
list_hyper <- list(
  lambda = lambda,
  v0_v = v0_v,
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

res1 <- navigm(Y, method = 'GM', 
                     version = 1,
                     list_hyper = list_hyper,
                     list_init = list_init, 
                     tol = 1e-3,
                     maxit = 1e5)

#res1$estimates$m_delta
plot(v0_v[2:16],res1$model_criterion$value[2:16], xlab = 'value of spike variance v0', ylab = 'AIC', type = 'l')
