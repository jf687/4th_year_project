

path = dirname(rstudioapi::getSourceEditorContext()$path)
source(paste0(path, "/functions.R"))
source(paste0(path, "/gmt.R"))

#
Tn <- 5
ts <- 1:Tn
P <- 20
reps <- 1

# at the first time point
#
aucs <- c()
precs <-  c()
recs <-  c()

for(c in 1:reps){
  set.seed(c)
  A0 <- matrix(rbinom(P*P, 1, 0.03), nrow = P)
  A0[lower.tri(A0)] = t(A0)[lower.tri(A0)]
  diag(A0) <- 0
  print(sum(A0[upper.tri(A0)])/sum(upper.tri(A0)))
  
  #
  # 3 edges become present over time; 1 becomes absent over time
  #
  edge_id <- which(A0 == T & upper.tri(A0), arr.ind = T)
  no_edge_id <- which(A0 == F & upper.tri(A0), arr.ind = T)
  npos <- 3; nneg <- 1
  pos_id <- matrix(no_edge_id[sample(1:nrow(no_edge_id),npos),], ncol=2)
  neg_id <- matrix(edge_id[sample(1:nrow(edge_id),nneg),], ncol=2)
  pos_id_sign_change <- (ceiling(Tn/2)-1):(ceiling(Tn/2)+1)
  neg_id_sign_change <- ceiling(Tn/2)
  As <- replicate(Tn, A0, simplify=FALSE)
  
  for (i in 1:npos) {
    for (t in (pos_id_sign_change[i] + 1):Tn) {
      As[[t]][pos_id[i,1],pos_id[i,2]] <-  As[[t]][pos_id[i,2],pos_id[i,1]] <- T
    }
  }
  
  for (i in 1:nneg) {
    for (t in (neg_id_sign_change[i] + 1):Tn) {
      As[[t]][neg_id[i,1],neg_id[i,2]] <- As[[t]][neg_id[i,2],neg_id[i,1]] <- F
    }
  }
  sparsity <- sapply(As, function(x)sum(x[upper.tri(x)])/sum(upper.tri(x)))
  print(sparsity)
  
  
  # simulate data
  #
  Ns <- sample(c(100,150,200,250), Tn, replace = T)
  ans <- lapply(1:Tn, function(t)
    generate_network_tan(n = Ns[t],
                         p = P,
                         A = As[[t]], 
                         empirical = T))
  Ys <- lapply(ans, function(x)x$Y)
  
  # inference
  #
  out <- gmt(Ys,
             ts,
             debug = T)
  
  FP <- 0
  TP <- 0
  tot <- 0
  

  for(l in list(1,3,4,5)){
    r <- out$estimates$m_deltas[[l]]
    
    bool_up <- upper.tri(r)
    
    labels <- As[[l]][bool_up]
    predictions <- r[bool_up]
    
    pred <- ROCR::prediction(predictions, labels)
    perf <- ROCR::performance(pred, measure = "tpr", x.measure = "fpr")
    auc_value <- ROCR::performance(pred, "auc")@y.values[[1]]
    
    prec_value <- precision( As[[l]], r > 0.5)
    rec_value <- recall(As[[l]], r > 0.5)
    
    aucs[[l]] <- c(aucs[[l]], auc_value)
    precs[[l]] <- c(precs[[l]], prec_value)
    recs[[l]] <- c(recs[[l]], rec_value)
    
    plot(perf, main = '', xlim = c(0, 1), ylim = c(0, 1), lwd = 2, colorize = TRUE)

  }
}

par(mfrow = c(1,2))
plot(out$debugs$vec_ELBO_CM, xlab = "iterations", ylab = "ELBOs at M-step")
plot(unlist(out$debugs$list_ELBO), xlab = "iterations", ylab = "ELBOs")

hist(out$estimates$mu_zeta[upper.tri(out$estimates$mu_zeta)], labels = T, xlab = "zeta", main = "")
hist(out$estimates$mu_beta[upper.tri(out$estimates$mu_beta)], labels = T, xlab = "beta", main = "")

out$estimates$mu_beta[pos_id]
out$estimates$mu_beta[neg_id]

sapply(ans, function(x)x$A[pos_id])
sapply(out$estimates$m_deltas, function(x)x[pos_id])

sapply(ans, function(x)x$A[neg_id])
sapply(out$estimates$m_deltas, function(x)x[neg_id])
