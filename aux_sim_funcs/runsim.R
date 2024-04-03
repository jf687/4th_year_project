
source('~/4th_year_project/Sims/ROC_calcs.R')
source('~/4th_year_project/Sims/Glasso_ROC_sim.R')
source('~/4th_year_project/GM/GM.R')
source('~/4th_year_project/GM/fun_GM.R')
source('~/4th_year_project/GM/fun_utils.R')
source('~/4th_year_project/GHS/GHS_example.R')
source('~/4th_year_project/BGlasso/BGlasso_example.R')



seed <- 123; 
set.seed(seed)
n <- 100
p <- 50
N <- 1
Y <- matrix(rnorm(n*p), nrow = n, ncol = p)
list_hyper <- list(lambda = 2,                   
                   v0 = 0.5,                  
                   v1 = 100,                   
                   a = 2,                   
                   b = 2,                   
                   ar = 1,                   
                   br = p)
list_init <- list(a_rho = 1,                  
                  b_rho = 1,                  
                  a_tau = 1,                  
                  b_tau = 1)


sf = huge::huge.generator(n=n, d=p,graph = 'scale-free', prob = 0.02)

adj_mat <- sf$theta
omega.true <- sf$omega
sigma.true <- sf$sigma

y = mvtnorm::rmvnorm(n, rep(0, p), sigma.true)
y = scale(y)



ns <- seq(from = 50, to = 1000, length.out = 5)
ps <- seq(from = 40, to = 500, length.out = 5)

BGlasso.op <- list()
GHS.op <- list()
GM.op <- list()
t.vals = seq(0,1,length.out = 100)
v0.vals = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8)
res <- select.v0.bic(y, list_hyper, list_init, v0.vals = seq(0,1,length.out = 100), t.vals = 
                       seq(0,1,length.out = 100),n,p)

for(n in ns){
  for(p in ps){
    BGlasso <- multiple_BGlasso(20,n,p) #BGlasso_ROC_sim(n = n, p = p)
    BGlasso.op <- c(BGlasso.op,BGlasso)
  }
}

save(BGlasso.op, file='BGlasso_op.rda')
#save(GHS.op, file='GHS_op.rda')
#save(GM.op, file='GM_op.rda')


n <- 100
p <- 50

v0s <- list(0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4)

v0.op.1 <- v0_grid(list_hyper, list_init, v0_list = v0s, save.data = F, n=200, p = 100)
save(v0.op.1, file='SSL_1.rda')
v0.op.2 <- v0_grid(list_hyper, list_init, v0_list = v0s, save.data = F, n=400, p = 200)
save(v0.op.2, file='SSL_2.rda')
v0.op.3 <- v0_grid(list_hyper, list_init, v0_list = v0s, save.data = F, n=750, p = 300)
save(v0.op.3, file='SSL_3.rda')


v0opt <- select.v0.bic(y, list_hyper, list_init,v0.vals,t.vals,n,p)
v0opt <- v0opt$v0.opt

list_hyper$v0 <- v0opt

res.ssl <- generate_GM_preds_labels(list_hyper, list_init)






