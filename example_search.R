if(!require(remotes)) install.packages("remotes")
remotes::install_github("XiaoyueXI/navigm")
library(navigm)
library(huge)

#generate data

set.seed(123)
L = huge.generator(d = 50, graph="hub")
# L = huge.generator(d = 50, graph="scale-free")
par(mfrow = c(1,2))
image(L$omega)

out.glasso = huge(L$data, method = "glasso")
out.select = huge.select(out.glasso,criterion = "ebic")

# 
v0_v <- seq(1e-2, 1, length.out = 16)
v1 <- 100
P <- ncol(L$data)

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

