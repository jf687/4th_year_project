<<<<<<< HEAD

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
=======
if(!require(remotes)) install.packages("remotes")
remotes::install_github("XiaoyueXI/navigm")
library(navigm)
library(huge)

>>>>>>> f713fc7ea2885f62707bf63bb1e62ec597a092c8
#generate data

<<<<<<< HEAD
#model selection using stars
out.select = huge.select(out.ct, criterion = "stars", stars.thresh = 0.05,rep.num=10)
plot(out.select)
=======
set.seed(123)
L = huge.generator(d = 50, graph="hub")
# L = huge.generator(d = 50, graph="scale-free")
par(mfrow = c(1,2))
image(L$omega)
>>>>>>> f713fc7ea2885f62707bf63bb1e62ec597a092c8

out.glasso = huge(L$data, method = "glasso")
out.select = huge.select(out.glasso,criterion = "ebic")
<<<<<<< HEAD
plot(out.select)
#=======
v0_v <- seq(1e-2, 1, length.out = 16)
v1 <- 1000
P <- ncol(Y)
#>>>>>>> ab795f0b846821b23a836809263fc49ad232f80b
0
=======

# 
v0_v <- seq(1e-2, 1, length.out = 16)
v1 <- 100
P <- ncol(L$data)

>>>>>>> f713fc7ea2885f62707bf63bb1e62ec597a092c8
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
    beta_rho = P
  )

res1 <- navigm(L$data, 
               method = 'GM', 
               version = 1,
               list_hyper = list_hyper,
               list_init = list_init, 
               tol = 1e-3,
               maxit = 1e5,
               full_output = T,
               debug = T)

res1$model_criterion

pred <- prediction(res1$estimates$m_delta[bool_up], as.matrix(L$theta)[bool_up])
perf <- performance(pred,"tpr","fpr")
plot(perf)
pred <- prediction(out.select$opt.icov[bool_up], as.matrix(L$theta)[bool_up])
perf <- performance(pred,"tpr","fpr")
plot(perf, add = T ,col = "red")

# image(res1$estimates$m_delta > 0.5)
# image(out.select$opt.icov)

<<<<<<< HEAD
#res1$estimates$m_delta
plot(v0_v[2:16],res1$model_criterion$value[2:16], xlab = 'value of spike variance v0', ylab = 'AIC', type = 'l')
=======
>>>>>>> f713fc7ea2885f62707bf63bb1e62ec597a092c8
