source('cutoff_curves.R')
source('GM.R')
source('fun_GM.R')
source('fun_utils.R')

# Test the code, for only the graphical lasso ------------------------------------

seed <- 123; 
set.seed(seed)
n <- 100
p <- 50
N <- 50
Y <- matrix(rnorm(n*p), nrow = n, ncol = p)
list_hyper <- list(lambda = 2,                   
                   v0 = 100,                  
                   v1 = 0.5,                   
                   a = 2,                   
                   b = 2,                   
                   ar = 1,                   
                   br = p)
list_init <- list(a_rho = 1,                  
                  b_rho = 1,                  
                  a_tau = 1,                  
                  b_tau = 1)



ggm.sf = GGM_gen(n,p)
omega.true = ggm.sf$omega.true
res.sim <- perform_ROC_simulation(omega.true, n, list_hyper, list_init, include.ssl=T,N)
# Plot ROC curve
cutoff <- get_AUC(res.sim, method = 'Glasso')
cutoff <- cutoff$cutoff
#plot_ROC(res.sim, include.ssl = T)
#plot_PRC(res.sim, include.ssl = T)

# Get cut-off AUC
#get_AUC(res.sim, method='Glasso')
#get_PRC_AUC(res.sim, method='Glasso')

#get_AUC(res.sim, method='SSL')
#get_PRC_AUC(res.sim, method='SSL')