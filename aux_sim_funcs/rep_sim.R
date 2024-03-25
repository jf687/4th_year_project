source('~/4th_year_project/Sims/ROC_calcs.R')
source('~/4th_year_project/Sims/Glasso_ROC_sim.R')
source('~/4th_year_project/GM/GM.R')
source('~/4th_year_project/GM/fun_GM.R')
source('~/4th_year_project/GM/fun_utils.R')
source('~/4th_year_project/GHS/GHS_example.R')
source('~/4th_year_project/BGlasso/BGlasso_example.R')

##TEST MARK##
list_hyper <- list(lambda = 2,                   
                   v0 = 0.5,                  
                   v1 = 100,                   
                   a = 2,                   
                   b = 2,                   
                   ar = 1,                   
                   br = NULL)
list_init <- list(a_rho = 1,                  
                  b_rho = 1,                  
                  a_tau = 1,                  
                  b_tau = 1)

ns <- seq(from = 100, to = 400, length.out = 2)
ps <- seq(from = 50, to = 100, length.out = 2)


v0.vals <- seq(from = 0.01, to = 1.3, length.out = 20)
t.vals <- seq(from = 0.01, to = 1, length.out = 100)

N <- 50

SSL.op <- list()
SSL.runtimes <- list()
GLasso.op <- list()
GLasso.runtimes <- list()
BGLasso.op <- list()
BGLasso.runtimes <- list()

GHS.op <- list()
GHS.runtimes <- list()

SSL.bool <- T
GLasso.bool <- T
BGLasso.bool <- T
GHS.bool <- T

for(n in list(100)){
  for(p in list(50)){
    runtime <- list()
    if(n>p){
      
      cat('=========n:',n,'p:',p,'==========')
       
      ###### data generation ######
      
      sf = huge::huge.generator(n=n, d=p,graph = 'scale-free', prob = 0.02)
      
      adj_mat <- sf$theta
      omega.true <- sf$omega
      sigma.true <- sf$sigma
      
      runtime$n <- n
      runtime$p <- p
      
      
      
      ###### SSL ######
      if(SSL.bool){
        
        start.time <- Sys.time()
  
        SSL <- GM_average(sf, N, list_hyper, list_init, v0.vals, t.vals)
        SSL.op <- c(SSL.op, SSL)
        
        end.time <- Sys.time()
        
        runtime$time <- end.time - start.time
        
        SSL.runtimes <- c(SSL.runtimes, runtime)
        
        print('#### SSL simulation complete ####')
      }
      ##### GLasso #####
      
      if(GLasso.bool){
        
        start.time <- Sys.time()
        
        GLasso <- GLasso_average(sf, N)
        GLasso.op <- c(GLasso.op, GLasso)
        
        end.time <- Sys.time()
        
        runtime$time <- end.time - start.time
        
        GLasso.runtimes <- c(GLasso.runtimes, runtime)
        
        print('#### GLasso simulation complete ####')
      }
      ##### BGLasso #####
      
      if(BGLasso.bool){
        start.time <- Sys.time()
        
        BGLasso <- BGLasso_average(sf, N)
        BGLasso.op <- c(BGLasso.op, BGLasso)
        
        end.time <- Sys.time()
        
        runtime$time <- end.time - start.time
        
        BGLasso.runtimes <- c(BGLasso.runtimes, runtime)
        
        print('#### BGLasso simulation complete ####')
      }
      
      if(GHS.bool){
        start.time <- Sys.time()
        
        GHS <- GHS_average(sf, N)
        GHS.op <- c(GHS.op, GHS)
        
        end.time <- Sys.time()
        
        runtime$time <- end.time - start.time
        
        GHS.runtimes <- c(GHS.runtimes, runtime)
        
        print('#### GHS simulation complete ####')
      }
 
    }
    
  }
}

save(SSL.op, file = 'SSL.outputs')
save(SSL.runtimes, file = 'SSL.runtimes')
save(GLasso.op, file = 'GLasso.outputs')
save(GLasso.runtimes, file = 'GLasso.runtimes')
save(BGLasso.op, file = 'BGLasso.outputs')
save(BGLasso.runtimes, file = 'BGLasso.runtimes')
save(GHS.op, file = 'GHS.outputs')
save(GHS.runtimes, file = 'GHS.runtimes')
