form_ROC_GGM <- function(nvec, p,n, plot = TRUE, t = 1e-5, x_lim = 0.02, y_lim = 0, title = 'title'){
  predictions <- c()
  labels <- c()
  for (i in 1:nvec){
    cat('=============',i,'=============\n')
    ggm.sf <- GGM_gen(n,p)
    omega.true <- ggm.sf$omega.true
    result <- glasso_sim(n,p,omega.true)
    #preds[[iteration]] <- abs(result$opt.icov[upper.tri(result$opt.icov)])
    #labels[[iteration]] <- omega.true[upper.tri(omega.true)] != 0
    predictions <- c(predictions, abs(result$opt.icov[upper.tri(result$opt.icov)]))
    labels <- c(labels, abs(omega.true[upper.tri(omega.true)]))
  }
  labels = abs(labels) > 1e-5
  #labels <- binarise_mat(labels, t)
  print(title)
  pred <- prediction(predictions, labels)
  auc_value <- performance(pred, "auc")@y.values[[1]]
  if (plot == TRUE){
    
    #Cairo(file=title,
    #      type="png",
    #      units="px", 
    #      width=500, 
    #      height=500, 
    #      pointsize=12, 
    #      dpi="auto")
    
    perf <- performance(pred, measure = "tpr", x.measure = "fpr")
    #x_lim <- perf@fp[-1]/ (perf@fp[-1] + perf@tn[-1])
    plot(perf, main = 'TPR vs FPR', xlim = c(0, x_lim), ylim = c(y_lim, 1), lwd = 2, colorize = TRUE)
    abline(0,1,lty = 'dashed')
    #labs(subtitle = paste("Area Under Curve (AUC):",auc_value)) #,'.  n = ',n, ', v0 = ',list_hyper$v0,', v1 = ',list_hyper$v1 ), side = 3, line = 0.5, cex = 0.8)
    
    
    perf.2 <- ROCR::performance(pred, measure = "prec", x.measure = "rec")
    plot(perf.2, main = 'Precision-Recall', lwd = 2, colorize = TRUE)
    #mtext(paste("Area Under Curve (AUC):",auc_value)) #,'.  n = ',n, ', v0 = ',list_hyper$v0,', v1 = ',list_hyper$v1 ), side = 3, line = 0.5, cex = 0.8)
    
    plot_pr(labels, predictions)
  }
  
  return(list(predictions, labels, auc_value, perf, perf.2))
}



plot_pr <- function(labels, predictions){
  
  # Assuming 'predictions' are the predicted probabilities and 'labels' are the true labels (0 or 1)
  roc_curve <- PRROC::roc.curve(labels, predictions,curve=TRUE)
  pr_curve <- PRROC::pr.curve(labels, predictions,curve=TRUE)
  # Plot Precision-Recall curve
  plot(roc_curve)
  plot(pr_curve, col = "blue", main = "Precision-Recall Curve", lwd = 2)
  
}

dat <- form_ROC_GGM(nvec = 1, p = 100,n = 200, plot = TRUE, t = 1e-5, x_lim = 0.05, y_lim = 0, title = 'title')
