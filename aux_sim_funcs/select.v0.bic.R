source('GM/GM.R')
source('GM/fun_utils.R')
source('GM/fun_GM.R')
select.v0.bic <- function(Y, list_hyper, list_init, v0.vals,t.vals,n,p){
  Omega.t.list <- list()
  Omega.opt.list <- list()
  
  for(j in 1:length(v0.vals)){
    cat('v0 = ',v0.vals[[j]])
    list_hyper$v0 <- v0.vals[j]
    out <- GM(Y, list_hyper, list_init)
    
    Omega <- out$Omega
    for(i in 1:length(t.vals)){
      Omega.t <- Omega
      Omega.t[out$m_delta < t.vals[i]] <- 0
      Omega.t.list[[i]] <-  Omega.t
    }
  
    
    BICs = unlist(lapply(Omega.t.list, FUN = function(s) bic_ssl(Y, s, n,p) ))
    t.opt = t.vals[which.min(BICs)]
    Omega.opt = Omega.t.list[[which.min(BICs)]]
    Omega.opt.list[[j]] <- Omega.opt
  }
  
  v0.BICs = unlist(lapply(Omega.opt.list, FUN = function(s) bic_ssl(Y, s, n,p) ))
  v0.opt = v0.vals[which.min(v0.BICs)]
  Omega.opt.opt = Omega.opt.list[[which.min(v0.BICs)]]
  
  return(list(v0.opt = v0.opt, Omega.opt = Omega.opt.opt, t.opt = t.opt))
}