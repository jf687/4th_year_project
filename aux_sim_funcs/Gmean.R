gmean_maximisaton = function(FPRs, TPRs){
  stopifnot(length(FPRs) == length(TPRs))
  max_g.mean = 0
  max_i = 0
  for(i in 1:length(FPRs)){
    g.mean <- sqrt(TPRs[[i]] * (1-FPRs[[i]]))
    if( g.mean >= max_g.mean){
      max_g.mean <- g.mean
      max_i <- i
    }
  }
  precision <- TPRs[[i]]/(TPRs[[i]] + FPRs[[i]])
  
  FNR <- 1 - TPRs[[i]]
  TNR <- 1 - FPRs[[i]]
  
  return(list(g.mean = g.mean, i = i, precision = precision, 
              FNR = FNR, TNR = TNR))
  
}