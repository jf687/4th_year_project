source('~/4th_year_project/Sims/ROC_calcs.R')

source('~/4th_year_project/GM/GM.R')
source('~/4th_year_project/GM/fun_GM.R')
source('~/4th_year_project/GM/fun_utils.R')

# Test the code, for only the graphical lasso ------------------------------------

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

ns <- seq(from = 50, to = 1000, by = 25)
ps <- seq(from = 40, to = 500, by = 20)

freq_run(ns = ns, ps = ps, save.data = T)

n <- 100
p <- 50

v0_grid(list_hyper, list_init, v0_list = v0s, save.data = T, n=200, p = 100)
v0_grid(list_hyper, list_init, v0_list = v0s, save.data = T, n=400, p = 200)
v0_grid(list_hyper, list_init, v0_list = v0s, save.data = T, n=750, p = 300)







