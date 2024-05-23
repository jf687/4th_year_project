library(reshape2)
library(patchwork)
library(igraph)
library(ggplot2)
library(network)
library(sna)
library(GGally)

get_and_print_edges <- function(a.mat, col.names, theta.mat=NULL) {
  # Function for printing all the edges in a graph with a layout that can be inserted into a latex table.
  # Also returns a data frame containing the edges
  # a.mat:          the adjacency matrix
  # col.names:      the names of the nodes in the graph
  # theta.mat:      the precision matrix. If included, the size the partial correlations are included in the table as well
  a.mat[which(diag(rep(1, ncol(a.mat))) == 1, arr.ind = T)] <- 0 # Make diagonal zero
  pairs <- which(a.mat[, ] == 1, arr.ind = T)
  df <- data.frame(t(apply(pairs, 1, sort))) # Sort so that the node in the pair whose name is first in the alphabet is first.
  df <- unique(df)
  names <- cbind(col.names[df[, 1]], col.names[df[, 2]])
  if (!is.null(theta.mat)){
    effect = round(-cov2cor(theta.mat)[cbind(df[,1], df[,2])], 5)
    names = cbind(names, effect)
    return(names)
  }
  for (i in 1:nrow(names)) {
    cat(names[i, 1], " & ", names[i, 2], " \\\\ \n")
  }
  return(names)
}