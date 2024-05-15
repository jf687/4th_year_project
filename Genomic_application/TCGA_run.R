

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

##load in data
source('~/4th_year_project/Genomic_application/format_data.R')


## install packages for plotting networks
if(F){
  install.packages("network")
  install.packages("sna")
}
## create function to plot graphs with the PPI matrix as input and threshold of 0.5

plot_graph <- function(stage_data, prec_matr, t , title = 'Proteomic graphical network for \n Stage I breast cancer patients'){
  bin_matrix <- abs(prec_matr) > t
  net <- network(bin_matrix, directed=FALSE)
  network.vertex.names(net) <- colnames(stage_data)
  plot(net, displaylabels=TRUE, main= title, label.cex = 0.55)
  
}

## Begin by running each stage seperately using the SSL model and plotting the resulting graph

stage_example_SSL <- function(stage_data = brca_dat_stageIII){
  p <- dim(stage)[[1]]
  #n <- dim(stage)[[2]]
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
  
  stage_data <- as.matrix(stage_data)
    
  res.ssl <- GM(stage_data, list_hyper, list_init)
  plot_graph(res.ssl$m_delta)
}

stage_example_GLasso <- function(stage_data = brca_dat_stageI, t = 0){
  n.lambda = 15
  res.glasso = huge(stage_data, method = 'glasso', nlambda=n.lambda, verbose = F)
  
  sim.obj = huge.select(res.glasso, criterion = 'stars')
  res.glasso.omega.opt = sim.obj$opt.icov
  print(sparsity(res.glasso.omega.opt))
  plot_graph(stage_data, res.glasso.omega.opt, t)
  plot_graph(stage_data, res.glasso.omega.opt, t+0.02)
  plot_graph(stage_data, res.glasso.omega.opt, t+0.04)
}

TCGA_GMT <- function(){
  
  Ys <- list(brca_dat_stageI, brca_dat_stageII, brca_dat_stageIII)
  out <- gmt(Ys,
             ts,
             debug = T)

  
  plot_graph(out$estimates$m_deltas[[1]], brca_dat_stageI, 0.5)
  plot_graph(out$estimates$m_deltas[[2]], brca_dat_stageI, 0.5)
  plot_graph(out$estimates$m_deltas[[3]], brca_dat_stageI, 0.5)
  
  
}
