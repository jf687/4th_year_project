source('~/4th_year_project/aux_sim_funcs/ROC_calcs.R')
source('~/4th_year_project/aux_sim_funcs/libraries.R')
source('~/4th_year_project/aux_sim_funcs/select.v0.bic.R')

source('~/4th_year_project/ROC_sim_functions/Glasso_ROC_sim.R')
source('~/4th_year_project/ROC_sim_functions/BGlasso_example.R')
source('~/4th_year_project/ROC_sim_functions/GHS_example.R')
source('~/4th_year_project/ROC_sim_functions/GM_ROC_sim.R')

source('~/4th_year_project/aux_method_funcs/GM/GM.R')
source('~/4th_year_project/aux_method_funcs/GM/fun_GM.R')
source('~/4th_year_project/aux_method_funcs/GM/fun_utils.R')
source('~/4th_year_project/aux_method_funcs/GHS/GHS.R')


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


v0.vals <- seq(from = 0.9, to = 1.9, length.out = 20)
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

SSL.bool <- F
GLasso.bool <- F
BGLasso.bool <- T
GHS.bool <- F

load('~/4th_year_project/HUGE_data/one.RData')
load('~/4th_year_project/HUGE_data/two.RData')
load('~/4th_year_project/HUGE_data/three.RData')

for(n in list(200)){
  for(p in list(100)){
    if(n>p){
      
      cat('=========n:',n,'p:',p,'==========')
       
      ###### data generation ######
      
      sf <-  sf200_100
      
      adj_mat <- sf$theta
      omega.true <- sf$omega
      sigma.true <- sf$sigma
      
      

      ###### SSL ######
      if(SSL.bool){
  
        SSL <- GM_average(sf, N, list_hyper, list_init, v0.vals, t.vals, plot = TRUE)
        SSL.op <- c(SSL.op, SSL)
            
        print('#### SSL simulation complete ####')
      }
      ##### GLasso #####
      
      if(GLasso.bool){

        GLasso <- GLasso_average(sf, N, plot = TRUE)
        GLasso.op <- c(GLasso.op, GLasso)

        
        print('#### GLasso simulation complete ####')
      }
      ##### BGLasso #####
      
      if(BGLasso.bool){

        
        BGLasso <- BGLasso_average(sf, N)
        BGLasso.op <- c(BGLasso.op, BGLasso)
        
        
        print('#### BGLasso simulation complete ####')
      }
      
      if(GHS.bool){

        GHS <- GHS_average(sf, N)
        GHS.op <- c(GHS.op, GHS)

        
        print('#### GHS simulation complete ####')
      }
 
    }
    
  }
}

save(SSL.op, file = 'results/SSL200_100_outputs.RData')

#save(GLasso.op, file = 'results/GLasso_outputs.RData')

save(BGLasso.op, file = 'results/BGLasso200_100_outputs.RData')

save(GHS.op, file = 'results/GHS200_100_outputs.RData')

for(n in list(200)){
  for(p in list(150)){
    if(n>p){
      
      cat('=========n:',n,'p:',p,'==========')
      
      ###### data generation ######
      
      sf <-  sf200_150
      
      adj_mat <- sf$theta
      omega.true <- sf$omega
      sigma.true <- sf$sigma
      
      
      
      ###### SSL ######
      if(SSL.bool){
        
        SSL <- GM_average(sf, N, list_hyper, list_init, v0.vals, t.vals, plot = TRUE)
        SSL.op <- c(SSL.op, SSL)
        
        print('#### SSL simulation complete ####')
      }
      ##### GLasso #####
      
      if(GLasso.bool){
        
        GLasso <- GLasso_average(sf, N, plot = TRUE)
        GLasso.op <- c(GLasso.op, GLasso)
        
        
        print('#### GLasso simulation complete ####')
      }
      ##### BGLasso #####
      
      if(BGLasso.bool){
        
        
        BGLasso <- BGLasso_average(sf, N)
        BGLasso.op <- c(BGLasso.op, BGLasso)
        
        
        print('#### BGLasso simulation complete ####')
      }
      
      if(GHS.bool){
        
        GHS <- GHS_average(sf, N)
        GHS.op <- c(GHS.op, GHS)
        
        
        print('#### GHS simulation complete ####')
      }
      
    }
    
  }
}

save(SSL.op, file = 'results/SSL200_150_outputs.RData')

#save(GLasso.op, file = 'results/GLasso_outputs.RData')

save(BGLasso.op, file = 'results/BGLasso200_150_outputs.RData')

save(GHS.op, file = 'results/GHS200_150_outputs.RData')

for(n in list(200)){
  for(p in list(200)){
    if(n>p){
      
      cat('=========n:',n,'p:',p,'==========')
      
      ###### data generation ######
      
      sf <-  sf200_200
      
      adj_mat <- sf$theta
      omega.true <- sf$omega
      sigma.true <- sf$sigma
      
      
      
      ###### SSL ######
      if(SSL.bool){
        
        SSL <- GM_average(sf, N, list_hyper, list_init, v0.vals, t.vals, plot = TRUE)
        SSL.op <- c(SSL.op, SSL)
        
        print('#### SSL simulation complete ####')
      }
      ##### GLasso #####
      
      if(GLasso.bool){
        
        GLasso <- GLasso_average(sf, N, plot = TRUE)
        GLasso.op <- c(GLasso.op, GLasso)
        
        
        print('#### GLasso simulation complete ####')
      }
      ##### BGLasso #####
      
      if(BGLasso.bool){
        
        
        BGLasso <- BGLasso_average(sf, N)
        BGLasso.op <- c(BGLasso.op, BGLasso)
        
        
        print('#### BGLasso simulation complete ####')
      }
      
      if(GHS.bool){
        
        GHS <- GHS_average(sf, N)
        GHS.op <- c(GHS.op, GHS)
        
        
        print('#### GHS simulation complete ####')
      }
      
    }
    
  }
}

save(SSL.op, file = 'results/SSL200_200_outputs.RData')

#save(GLasso.op, file = 'results/GLasso_outputs.RData')

save(BGLasso.op, file = 'results/BGLasso200_200_outputs.RData')

save(GHS.op, file = 'results/GHS200_200_outputs.RData')



