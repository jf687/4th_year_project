

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

stage_example_SSL <- function(stage_data = brca_dat_stageIII, filename = 'SSL_stageIII.RData'){
  p <- dim(stage_data)[[2]]
  n <- dim(stage_data)[[1]]
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
  v0.vals <- seq(from = 0.5, to = 2.5, length.out = 18)
  t.vals <- seq(from = 0.01, to = 1, length.out = 100)
  
  v0.selection <- select.v0.bic(stage_data, list_hyper, list_init, v0.vals,t.vals,n=140,p=131)
  list_hyper$v0 <- v0.selection$v0.opt
  thresh <- v0.selection$t.opt
  

  res.ssl <- GM(stage_data, list_hyper, list_init)
  #save(res.ssl, file = filename)
  plot_graph(stage_data, res.ssl$m_delta, thresh)
  browser()
  nedges <- which(res.ssl$m_delta<thresh)
  res.ssl$Omega[nedges]<- 0
  return(list(PPI = res.ssl$m_delta, omega.opt = res.ssl$Omega, thresh =  thresh))
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
  
  
  
  v0.selection <- ave_v0()
  thresh <- v0.selection$t.opt
  v0.opt <- v0.selection$v0.opt
  
  ts <- 1:length(Ys)
  
  out <- gmt(Ys,
             ts,
             debug = T,
             set_v0 = v0.opt)
  
  nedges <- lapply(out$estimates$m_deltas, function(x) which(x<thresh))
  browser()
  for(net in ts){
    out$estimates$Omegas[nedges[[net]]] <- 0
  }
  browser()
  
  omega_opts <- list(out$estimates$Omegas[[1]], out$estimates$Omegas[[2]], out$estimates$Omegas[[3]])
  browser()
  
  plot_graph(brca_dat_stageI, out$estimates$m_deltas[[1]] , thresh, title = 'Proteomic network of stage I breast carcinoma patients')
  plot_graph(brca_dat_stageII, out$estimates$m_deltas[[2]],  thresh, title = 'Proteomic network of stage II breast carcinoma patients')
  plot_graph(brca_dat_stageIII, out$estimates$m_deltas[[3]], thresh, title = 'Proteomic network of stage III breast carcinoma patients')
  
  save(out, file = 'SSLx_TCGA.RData')
  
  return(list(PPIs = out$estimates$m_deltas, omega_opts = omega_opts, thresh = thresh ))
}

all_stages_SSL <- function(){
  SSL_I <- stage_example_SSL(brca_dat_stageI, 'SSL_stageI.RData')
  SSL_II <- stage_example_SSL(brca_dat_stageII, 'SSL_stageII.RData')
  SSL_III <- stage_example_SSL(brca_dat_stageIII, 'SSL_stageIII.RData')
  SSL_x <- TCGA_GMT()
  save(SSL_I, file = "Genomic_application/data/SSL_stageI.RData")
  save(SSL_II, file = "Genomic_application/data/SSL_stageII.RData")
  save(SSL_III, file = "Genomic_application/data/SSL_stageIII.RData")
  save(SSL_x, file = "Genomic_application/data/SSLx_TCGA.RData")
}

ave_v0 <- function(){
  list_hyper <- list(lambda = 2,                   
                     v0 = 0.2,                  
                     v1 = 100,                   
                     a = 2,                   
                     b = 2,                   
                     ar = 1,                   
                     br = NULL)
  list_init <- list(a_rho = 1,                  
                    b_rho = 1,                  
                    a_tau = 1,                  
                    b_tau = 1)
  v0.vals <- seq(from = 0.5, to = 2.5, length.out = 25)
  t.vals <- seq(from = 0.01, to = 1, length.out = 100)
  
  list_hyper$br = dim(brca_dat_stageI)[[2]]
  v0I <- select.v0.bic(brca_dat_stageI, list_hyper, list_init, v0.vals,t.vals,n=140,p=dim(brca_dat_stageI)[[2]])
  
  list_hyper$br = dim(brca_dat_stageII)[[2]]
  v0II <- select.v0.bic(brca_dat_stageII, list_hyper, list_init, v0.vals,t.vals,n=140,p=dim(brca_dat_stageII)[[2]])
  
  list_hyper$br = dim(brca_dat_stageIII)[[2]]
  v0III <- select.v0.bic(brca_dat_stageIII, list_hyper, list_init, v0.vals,t.vals,n=140,p=dim(brca_dat_stageIII)[[2]])
  browser()
  return(list(v0.opt = mean(v0I$v0.opt, v0II$v0.opt, v0III$v0.opt), t.opt = mean(v0I$t.opt, v0II$t.opt, v0III$t.opt)))
  }
