

##load in methods
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

source('~/4th_year_project/extension/gmt.R')
source('~/4th_year_project/extension/functions.R')
##load in data
source('~/4th_year_project/Genomic_application/format_data.R')


## install packages for plotting networks
if(F){
  install.packages("network")
  install.packages("sna")
}
library(network)
library(sna)
## create function to plot graphs with the PPI matrix as input and threshold of 0.5

plot_graph <- function(stage_data, prec_matr, t = 0.5, title = 'Proteomic graphical network for \n Stage I breast cancer patients'){
  bin_matrix <- abs(prec_matr) > t
  net <- network(bin_matrix, directed=FALSE)
  network.vertex.names(net) <- colnames(stage_data)
  plot(net, displaylabels=TRUE, main= title, label.cex = 0.55)
  
}

## Begin by running each stage seperately using the SSL model and plotting the resulting graph

stage_example_SSL <- function(stage_data = brca_dat_stageIII){
  p <- dim(stage_data)[[2]]
  n <- dim(stage_data)[[2]]
  list_hyper <- list(lambda = 2,                   
                     v0 = 0.2,                  
                     v1 = 100,                   
                     a = 2,                   
                     b = 2,                   
                     ar = 1,                   
                     br = p)
  list_init <- list(a_rho = 1,                  
                    b_rho = 1,                  
                    a_tau = 1,                  
                    b_tau = 1)
    
  res.ssl <- GM(stage_data, list_hyper, list_init)
  browser()
  print(sparsity(res.ssl$m_delta, strict = F, 0.5))
  plot_graph(stage_data, res.ssl$m_delta, 0.5)
}

stage_example_GLasso <- function(stage_data = brca_dat_stageI, t = 0){
  n.lambda = 15
  res.glasso = huge(stage_data, method = 'glasso', nlambda=n.lambda, verbose = F)
  
  sim.obj = huge.select(res.glasso, criterion = 'stars')
  res.glasso.omega.opt = sim.obj$opt.icov
  browser()
  print(sparsity(res.glasso.omega.opt))
  plot_graph(stage_data, res.glasso.omega.opt, t)
  #plot_graph(stage_data, res.glasso.omega.opt, t+0.02)
  #plot_graph(stage_data, res.glasso.omega.opt, t+0.04)
}

TCGA_GMT <- function(){
  
  Ys <- list(brca_dat_stageI, brca_dat_stageII, brca_dat_stageIII)
  ts <- 1:length(Ys)
  out <- gmt(Ys,
             ts,
             debug = T)
  browser()
  
  
  plot_graph(brca_dat_stageI, out$estimates$m_deltas[[1]] , 0.5, title = 'Proteomic network of stage I breast carcinoma patients')
  plot_graph(brca_dat_stageII, out$estimates$m_deltas[[2]],  0.5, title = 'Proteomic network of stage II breast carcinoma patients')
  plot_graph(brca_dat_stageIII, out$estimates$m_deltas[[3]], 0.5, title = 'Proteomic network of stage III breast carcinoma patients')
  
  
}

#Ys <- list(brca_dat_stageI, brca_dat_stageII, brca_dat_stageIII)
#ts <- 1:length(Ys)
#out <- gmt(Ys,
 #          ts,
  #         debug = T)

