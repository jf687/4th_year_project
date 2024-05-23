rm(list=ls())
source('aux_sim_funcs/ROC_calcs.R')
source('Genomic_application/help_funcs.R')
load('Genomic_application/data/BRCA_dat_formatted.RData')


# Investigate TCGA proteomic results

# Load and format output of SSL and SSLx -----------------------------
load('Genomic_application/data/SSL_stageI.RData')
omega.ssl.i = -cov2cor(res.sslI$omega.opt) # Convert to partial correlations
diag(omega.ssl.i) = 1
load('Genomic_application/data/SSL_stageII.RData')
omega.ssl.ii = -cov2cor(res.sslII$omega.opt) # Convert to partial correlations
diag(omega.ssl.ii) = 1
load('Genomic_application/data/SSL_stageIII.RData')
omega.ssl.iii = -cov2cor(res.sslIII$omega.opt) # Convert to partial correlations
diag(omega.ssl.iii) = 1
load('Genomic_application/data/SSLx_TCGA.RData')
omega.ssl.x.i = -cov2cor(omega_opts[[1]])
diag(omega.ssl.x.i) = 1
omega.ssl.x.ii = -cov2cor(omega_opts[[2]])
diag(omega.ssl.x.ii) = 1
omega.ssl.x.iii = -cov2cor(omega_opts[[3]])
diag(omega.ssl.x.iii) = 1

# Check properties of networks ---------------------------------------

# Sparsity
# SSL estimates
sparsity(omega.ssl.i) # 0.06564885
sparsity(omega.ssl.ii) # 0.09782736
sparsity(omega.ssl.iii) # 0.06353494
# SSLx estimates
sparsity(omega.ssl.x.i) # 0.03957722
sparsity(omega.ssl.x.ii) # 0.1000587
sparsity(omega.ssl.x.iii) # 0.04110393

# Edge agreement 
# SSL estimates
confusion.matrix(omega.ssl.i!=0, omega.ssl.ii!=0)
#      [,1] [,2]
#[1,]  222  611
#[2,]  337 7345
confusion.matrix(omega.ssl.i!=0, omega.ssl.iii!=0)
#[,1] [,2]
#[1,]  170  371
#[2,]  389 7585
confusion.matrix(omega.ssl.ii!=0, omega.ssl.iii!=0)
#[,1] [,2]
#[1,]  213  328
#[2,]  620 7354

# SSLx estimates 
confusion.matrix(omega.ssl.x.i!=0, omega.ssl.x.ii!=0)
#[,1] [,2]
#[1,]  315  537
#[2,]   22 7641
confusion.matrix(omega.ssl.x.i!=0, omega.ssl.x.iii!=0)
#[,1] [,2]
#[1,]  315   35
#[2,]   22 8143
confusion.matrix(omega.ssl.x.ii!=0, omega.ssl.x.iii!=0)
#[,1] [,2]
#[1,]  315   35
#[2,]  537 7628

# Is it the same edges they all agree on?
(sum(omega.ssl.x.i!=0 & omega.ssl.x.ii!=0 & omega.ssl.x.iii!=0)-nrow(omega.ssl.i))/2
# 315
#Yes

# More similarity for SSLx

# Plot Matthews correlation coefficient to visualise similarity ---------------------

res.list = list(omega.ssl.i, omega.ssl.ii, omega.ssl.iii, omega.ssl.x.i, omega.ssl.x.ii, omega.ssl.x.iii)
MCC.mat = sapply(res.list, function(x) sapply(res.list, function(y) MCC(x!=0,y!=0))) # A matrix of pairwise MCC scores
MCC.mat[!lower.tri(MCC.mat)]=NA
rownames(MCC.mat)=colnames(MCC.mat) = c('Stage I (SSL)', 'Stage II (SSL)', 'Stage III (SSL)', 'Stage I (SSLx)', 'Stage II (SSLx)', 'Stage III (SSLx)')

data_melt <- melt(MCC.mat,na.rm=T)  
ggp <- ggplot(data_melt, aes(Var1, Var2)) + geom_tile(aes(fill = value))+
  scale_fill_gradient2(low = "white", high = "red", limit = c(0,1), space = "Lab",name="MCC") +
  theme_minimal()+ #theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))+
  coord_fixed()+xlab('')+ylab('')+theme(text = element_text(size = 20),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                                        plot.title=ggplot2::element_text(size=25,hjust=0.5)) #+ggtitle('stabJGL')
ggp 

pdf("Genomic_application/plots/MCC_plot.pdf", 8,5)
ggp 
dev.off()

# Look at top hub proteins of SSLx estimates ------------------------------------------

# Make igraph objects for easier computations
g.stage.i = igraph::graph.adjacency(omega.ssl.x.i!=0, mode='undirected', diag=F)
g.stage.ii = igraph::graph.adjacency(omega.ssl.x.ii!=0, mode='undirected', diag=F)
g.stage.iii = igraph::graph.adjacency(omega.ssl.x.iii!=0, mode='undirected', diag=F)
p=ncol(omega.ssl.i)

# Make data frames with the node degree of each protein, for each stage
df.degree = data.frame(degree=c(igraph::degree(g.stage.i), igraph::degree(g.stage.ii),igraph::degree(g.stage.iii)), 
                       Group=factor(c(rep('Stage I', p), rep('Stage II', p),rep('Stage III', p))))

df.degree.stage.i = data.frame(protein=colnames(brca_dat[[1]]),gene=mapping.frame$gene[match(colnames(brca_dat[[1]]),mapping.frame$protein)], 
                            degree=df.degree[which(df.degree$Group=='Stage I'),1])
df.degree.stage.ii = data.frame(protein=colnames(brca_dat[[2]]),
                            gene=mapping.frame$gene[match(colnames(brca_dat[[2]]),mapping.frame$protein)], 
                            degree=df.degree[which(df.degree$Group=='Stage II'),1])
df.degree.stage.iii = data.frame(protein=colnames(brca_dat[[3]]),gene=mapping.frame$gene[match(colnames(brca_dat[[3]]),mapping.frame$protein)], 
                                 degree=df.degree[which(df.degree$Group=='Stage III'),1])

# Order according to node degree
df.degree.stage.i.ordered = df.degree.stage.i[rev(order(df.degree.stage.i$degree)),]
df.degree.stage.iii.ordered = df.degree.stage.iii[rev(order(df.degree.stage.iii$degree)),]
df.degree.stage.ii.ordered = df.degree.stage.ii[rev(order(df.degree.stage.ii$degree)),]

# Make table of top 10% hubs
top.frac=0.9 
df.degree.stage.i.ordered.top = df.degree.stage.i.ordered[which(df.degree.stage.i.ordered$degree>quantile(df.degree.stage.i.ordered$degree,top.frac)),]
df.degree.stage.iii.ordered.top = df.degree.stage.iii.ordered[which(df.degree.stage.iii.ordered$degree>quantile(df.degree.stage.iii.ordered$degree,top.frac)),]
df.degree.stage.ii.ordered.top = df.degree.stage.ii.ordered[which(df.degree.stage.ii.ordered$degree>quantile(df.degree.stage.ii.ordered$degree,top.frac)),]

# Identify which hubs that are unique to that stage, to mark them in table.
which.unique.stage.i = which(! df.degree.stage.i.ordered.top$protein %in% c(df.degree.stage.iii.ordered.top$protein,df.degree.stage.ii.ordered.top$protein))
which.unique.stage.ii = which(! df.degree.stage.ii.ordered.top$protein %in% c(df.degree.stage.iii.ordered.top$protein,df.degree.stage.i.ordered.top$protein))
which.unique.stage.iii = which(! df.degree.stage.iii.ordered.top$protein %in% c(df.degree.stage.ii.ordered.top$protein,df.degree.stage.i.ordered.top$protein))
which.common.stage.i = which(df.degree.stage.i.ordered.top$protein %in% df.degree.stage.iii.ordered.top$protein & df.degree.stage.i.ordered.top$protein %in% df.degree.stage.ii.ordered.top$protein)
which.common.stage.ii = which(df.degree.stage.ii.ordered.top$protein %in% df.degree.stage.i.ordered.top$protein & df.degree.stage.ii.ordered.top$protein %in% df.degree.stage.iii.ordered.top$protein)
which.common.stage.iii= which(df.degree.stage.iii.ordered.top$protein %in% df.degree.stage.i.ordered.top$protein & df.degree.stage.iii.ordered.top$protein %in% df.degree.stage.ii.ordered.top$protein)

# Maximum number of top nodes across stages
dim.neighs.top = max(c(length(df.degree.stage.i.ordered.top$protein),length(df.degree.stage.iii.ordered.top$protein),length(df.degree.stage.ii.ordered.top$protein)))

# Print the top hub proteins, with the corresponding gene name, along with the node degree. Print for latex table formatting
for(i in 1:dim.neighs.top){
  if(is.na(df.degree.stage.i.ordered.top[i,1])) { cat(' & & && ') }
  else {
    if(i %in% which.common.stage.i){
      cat(paste0('\\textbf{',df.degree.stage.i.ordered.top$protein[i], '} & \\emph{ ',df.degree.stage.i.ordered.top$gene[i], '} & ', df.degree.stage.i.ordered.top$degree[i], ' && '))
    }
    else if(i %in% which.unique.stage.i){
      cat(paste0('\\color{red}{',df.degree.stage.i.ordered.top$protein[i], '} & \\emph{ ',df.degree.stage.i.ordered.top$gene[i], '} & ', df.degree.stage.i.ordered.top$degree[i], ' && '))
    }
    else { cat(paste0(df.degree.stage.i.ordered.top$protein[i], ' & \\emph{',df.degree.stage.i.ordered.top$gene[i], '} & ', df.degree.stage.i.ordered.top$degree[i], ' && ')) }
  }
  if(is.na(df.degree.stage.ii.ordered.top[i,1])) { cat(' & & && ')  }
  else {
    if(i %in% which.common.stage.ii){
      cat(paste0('\\textbf{',df.degree.stage.ii.ordered.top$protein[i], '} & \\emph{ ',df.degree.stage.ii.ordered.top$gene[i], '} & ', df.degree.stage.ii.ordered.top$degree[i], ' && '))
    }
    else if(i %in% which.unique.stage.ii){
      cat(paste0('\\color{red}{',df.degree.stage.ii.ordered.top$protein[i], '} & \\emph{ ',df.degree.stage.ii.ordered.top$gene[i], '} & ', df.degree.stage.ii.ordered.top$degree[i], ' && '))
    }
    else { cat(paste0(df.degree.stage.ii.ordered.top$protein[i], ' & \\emph{',df.degree.stage.ii.ordered.top$gene[i], '} & ', df.degree.stage.ii.ordered.top$degree[i], ' && ')) }
  }
  if(is.na(df.degree.stage.iii.ordered.top[i,1])) { cat(' & & \\\\ \n') }
  else {
    if(i %in% which.common.stage.iii){
      cat(paste0('\\textbf{',df.degree.stage.iii.ordered.top$protein[i], '} & \\emph{ ',df.degree.stage.iii.ordered.top$gene[i], '} & ', df.degree.stage.iii.ordered.top$degree[i],' \\\\ \n  '))
    }
    else if(i %in% which.unique.stage.iii){
      cat(paste0('\\color{red}{',df.degree.stage.iii.ordered.top$protein[i], '} & \\emph{ ',df.degree.stage.iii.ordered.top$gene[i], '} & ', df.degree.stage.iii.ordered.top$degree[i], ' \\\\ \n  '))
    }
    else { cat(paste0(df.degree.stage.iii.ordered.top$protein[i], ' & \\emph{',df.degree.stage.iii.ordered.top$gene[i], '} & ', df.degree.stage.iii.ordered.top$degree[i], ' \\\\ \n  ')) }
  }
}


# Write edge lists to files --------------------------------

edges.stage.i = get_and_print_edges(omega.ssl.x.i!=0,colnames(brca_dat[[1]]),omega.ssl.x.i)
colnames(edges.stage.i) =  c("Protein1", "Protein2", "PartialCor")
edges.stage.iii = get_and_print_edges(omega.ssl.x.iii!=0,colnames(brca_dat[[2]]),omega.ssl.x.iii)
colnames(edges.stage.iii) =  c("Protein1", "Protein2", "PartialCor")
edges.stage.ii = get_and_print_edges(omega.ssl.x.ii!=0,colnames(brca_dat[[3]]),omega.ssl.x.ii)
colnames(edges.stage.ii) =  c("Protein1", "Protein2", "PartialCor")
write.csv(edges.stage.i, file = "Genomic_application/edge_lists/edges_stageI.csv", row.names = F,quote=F)
write.csv(edges.stage.iii, file = "Genomic_application/edge_lists/edge_stageIII.csv", row.names = F,quote=F)
write.csv(edges.stage.ii, file = "Genomic_application/edge_lists/edges_stageII.csv", row.names = F,quote=F)


# Plot resulting networks ------------------------------

# Plot with common edges marked
unique.list = list()
unique.list[[1]] = (omega.ssl.x.i!=0) + 0  # If edge in stage.i
unique.list[[1]][which(omega.ssl.x.i!=0 & omega.ssl.x.ii!=0 & omega.ssl.x.iii!=0)] = 3 # Present in all
unique.list[[2]] = (omega.ssl.x.ii!=0) + 0  # If edge in stage.ii
unique.list[[2]][which(omega.ssl.x.i!=0 & omega.ssl.x.ii!=0 & omega.ssl.x.iii!=0)] = 3 # Present in all
unique.list[[3]] = (omega.ssl.x.iii!=0) + 0  # If edge in stage.iii
unique.list[[3]][which(omega.ssl.x.i!=0 & omega.ssl.x.ii!=0 & omega.ssl.x.iii!=0)] = 3 # Present in all

# Get layout for Stage I
set.seed(123)
net.layout = network::network(unique.list[[1]],directed=F, ignore.eval=F,names.eval='weights')
x.layout = sna::gplot.layout.fruchtermanreingold(net.layout, NULL)

nets.layout=list()
net.layout.stage.i = network::network(unique.list[[1]],directed=F, ignore.eval=F,names.eval='weights')
network::set.edge.attribute(net.layout.stage.i, "color", c("black", "grey75","red", "aquamarine4")[(net.layout.stage.i %e% "weights")+1])
net.layout.stage.i %v% "x" = x.layout[, 1]
net.layout.stage.i %v% "y" = x.layout[, 2]
nets.layout[[1]] = GGally::ggnet2(net.layout.stage.i,node.size = 1.5, edge.size = 0.25,alpha=0.9,mode = c('x','y'),color = 'royalblue1', edge.color = 'color')+
  ggplot2::ggtitle('Stage I') + ggplot2::theme(plot.title=ggplot2::element_text(size=15,hjust=0.5))
net.layout.stage.ii = network::network(unique.list[[2]],directed=F, ignore.eval=F,names.eval='weights')
network::set.edge.attribute(net.layout.stage.ii, "color", c("black", "grey75","red", "aquamarine4")[(net.layout.stage.ii %e% "weights")+1])
net.layout.stage.ii %v% "x" = x.layout[, 1]
net.layout.stage.ii %v% "y" = x.layout[, 2]
nets.layout[[2]] = GGally::ggnet2(net.layout.stage.ii,node.size = 1.5, edge.size = 0.25,alpha=0.9,mode = c('x','y'),color = 'royalblue1', edge.color = 'color')+
  ggplot2::ggtitle('Stage II') + ggplot2::theme(plot.title=ggplot2::element_text(size=15,hjust=0.5))
net.layout.stage.iii = network::network(unique.list[[3]],directed=F, ignore.eval=F,names.eval='weights')
network::set.edge.attribute(net.layout.stage.iii, "color", c("black", "grey75","red", "aquamarine4")[(net.layout.stage.iii %e% "weights")+1])
net.layout.stage.iii %v% "x" = x.layout[, 1]
net.layout.stage.iii %v% "y" = x.layout[, 2]
nets.layout[[3]] = GGally::ggnet2(net.layout.stage.iii,node.size = 1.5, edge.size = 0.25,alpha=0.9,mode = c('x','y'),color = 'royalblue1', edge.color = 'color')+
  ggplot2::ggtitle('Stage III') + ggplot2::theme(plot.title=ggplot2::element_text(size=15,hjust=0.5))

pdf("Genomic_application/plots/TCGA_BRCA_networks.pdf", 10,4)
print(gridExtra::grid.arrange(grobs=nets.layout,ncol=3))
dev.off()



# Check evidence in STRING  -----------------------------

# Get edges with gene names
edges.stage.i.gene.all = get_and_print_edges(omega.ssl.x.i!=0,mapping.frame$gene[match(colnames(brca_dat[[1]]), mapping.frame$protein)])
edges.stage.ii.gene.all = get_and_print_edges(omega.ssl.x.ii!=0,mapping.frame$gene[match(colnames(brca_dat[[2]]), mapping.frame$protein)])
edges.stage.iii.gene.all = get_and_print_edges(omega.ssl.x.iii!=0,mapping.frame$gene[match(colnames(brca_dat[[3]]), mapping.frame$protein)])
edges.stage.i.gene = edges.stage.i.gene.all[-which(duplicated(edges.stage.i.gene.all)),]
edges.stage.ii.gene = edges.stage.ii.gene.all[-which(duplicated(edges.stage.ii.gene.all)),]
edges.stage.iii.gene = edges.stage.iii.gene.all[-which(duplicated(edges.stage.iii.gene.all)),]
colnames(edges.stage.i.gene) =  c("Gene1", "Gene2")
colnames(edges.stage.iii.gene) =  c("Gene1", "Gene2")
colnames(edges.stage.ii.gene) =  c("Gene1", "Gene2")


# Print unique names for STRING
cat(unique(mapping.frame$gene), sep='\n')
# We get the network that you can view in STRING here: https://string-db.org/cgi/network?taskId=bcS0en9EYBdt&sessionId=bm3Yk1Wzg2Mb
string.output = 'https://string-db.org/cgi/generatetaskspecificdownloadfile?taskId=bcS0en9EYBdt&downloadDataFormat=tsv_short&cpnonce=bzVALeqOthCt&downloadFileName=string_interactions_short.tsv'
string.rppa = utils::read.csv(string.output, sep = "\t")
# Sort so that gene pairs are given in alphabetical order.
string.rppa[, 1:2] = t(apply(string.rppa[, 1:2], 1, sort)) # First gene in alphabet is always in col 1
colnames(string.rppa)[1:2] =  c("Gene1", "Gene2")

# Combine edges from STRING and SSLx into one table, so that we can check how many edges occur twice in the table 
# Stage I
df.all.edges.stage.i =  rbind(edges.stage.i.gene, string.rppa[, 1:2])
ids.in.string = which(duplicated(df.all.edges.stage.i)) - nrow(edges.stage.i.gene) # True the second time an edge occurs (which is the duplicate). We subtract nrow to get the index of the edge in STRING df.
# Count all pairs with these gene pairs XX HERE
edges.stage.i.gene.evidence = string.rppa[ids.in.string,1:2]
n.extra.pairs.evidence.stage.i = 0
for(i in 1:nrow(edges.stage.i.gene.evidence)){
  n.rep = 0
  for(j in 1:nrow(edges.stage.i.gene.all)){
    if(all(edges.stage.i.gene.all[j,] == edges.stage.i.gene.evidence[i,])){
      n.rep = n.rep + 1
    }
  }
  n.extra.pairs.evidence.stage.i = n.extra.pairs.evidence.stage.i + (n.rep-1)
}
# Number of edges with evidence in STRING database
length(ids.in.string)+n.extra.pairs.evidence.stage.i # 13
(length(ids.in.string)+n.extra.pairs.evidence.stage.i )/nrow(edges.stage.i.gene.all) # 0.03857567
# Also the amount of unique edges
length(ids.in.string)/nrow(edges.stage.i.gene) # 0.04

# Stage II
df.all.edges.stage.ii =  rbind(edges.stage.ii.gene, string.rppa[, 1:2])
ids.in.string = which(duplicated(df.all.edges.stage.ii)) - nrow(edges.stage.ii.gene) # True the second time an edge occurs (which is the duplicate). We subtract nrow to get the index of the edge in STRING df.
# Count all pairs with these gene pairs
edges.stage.ii.gene.evidence = string.rppa[ids.in.string,1:2]
n.extra.pairs.evidence.stage.ii = 0
for(i in 1:nrow(edges.stage.ii.gene.evidence)){
  n.rep = 0
  for(j in 1:nrow(edges.stage.ii.gene.all)){
    if(all(edges.stage.ii.gene.all[j,] == edges.stage.ii.gene.evidence[i,])){
      n.rep = n.rep + 1
    }
  }
  n.extra.pairs.evidence.stage.ii = n.extra.pairs.evidence.stage.ii + (n.rep-1)
}
# Number of edges with evidence in STRING database
length(ids.in.string) + n.extra.pairs.evidence.stage.ii # 23
(length(ids.in.string)+n.extra.pairs.evidence.stage.ii)/nrow(edges.stage.ii.gene.all) # 0.02699531
# Also the amount of unique edges
length(ids.in.string)/nrow(edges.stage.ii.gene) # 0.02417303

# Stage III
df.all.edges.stage.iii =  rbind(edges.stage.iii.gene, string.rppa[, 1:2])
ids.in.string = which(duplicated(df.all.edges.stage.iii)) - nrow(edges.stage.iii.gene) # True the second time an edge occurs (which is the duplicate). We subtract nrow to get the index of the edge in STRING df.
# Count all pairs with these gene pairs 
edges.stage.iii.gene.evidence = string.rppa[ids.in.string,1:2]
n.extra.pairs.evidence.stage.iii = 0
for(i in 1:nrow(edges.stage.iii.gene.evidence)){
  n.rep = 0
  for(j in 1:nrow(edges.stage.iii.gene.all)){
    if(all(edges.stage.iii.gene.all[j,] == edges.stage.iii.gene.evidence[i,])){
      n.rep = n.rep + 1
    }
  }
  n.extra.pairs.evidence.stage.iii = n.extra.pairs.evidence.stage.iii + (n.rep-1)
}
# Number of edges with evidence in STRING database
length(ids.in.string)+n.extra.pairs.evidence.stage.iii # 13
(length(ids.in.string)+n.extra.pairs.evidence.stage.iii)/nrow(edges.stage.iii.gene.all) # 0.03714286
# Also the amount of unique edges
length(ids.in.string)/nrow(edges.stage.iii.gene) # 0.03846154


# Find the amount of evidence for naive SSL ----------------------------------------

# Find graphs
g.stage.i.naive= igraph::graph.adjacency(omega.ssl.i!=0, mode='undirected', diag=F)
g.stage.ii.naive= igraph::graph.adjacency(omega.ssl.ii!=0, mode='undirected', diag=F)
g.stage.iii.naive= igraph::graph.adjacency(omega.ssl.iii!=0, mode='undirected', diag=F)
df.degree.naive = data.frame(degree=c(igraph::degree(g.stage.i.naive), igraph::degree(g.stage.ii.naive),igraph::degree(g.stage.iii.naive)), 
                             Group=factor(c(rep('Stage I', p), rep('Stage II', p),rep('Stage III', p))))


# Printing top list
top.frac=0.9
df.degree.stage.i.naive = data.frame(protein=colnames(brca_dat[[1]]),gene=mapping.frame$gene[match(colnames(brca_dat[[1]]),mapping.frame$protein)], 
                                     degree=df.degree.naive[which(df.degree.naive$Group=='Stage I'),1])
df.degree.stage.ii.naive = data.frame(protein=colnames(brca_dat[[2]]),
                                      gene=mapping.frame$gene[match(colnames(brca_dat[[2]]),mapping.frame$protein)], 
                                      degree=df.degree.naive[which(df.degree.naive$Group=='Stage II'),1])
df.degree.stage.iii.naive = data.frame(protein=colnames(brca_dat[[3]]),gene=mapping.frame$gene[match(colnames(brca_dat[[3]]),mapping.frame$protein)], 
                                       degree=df.degree.naive[which(df.degree.naive$Group=='Stage III'),1])

df.degree.stage.i.ordered.naive = df.degree.stage.i.naive[rev(order(df.degree.stage.i.naive$degree)),]
df.degree.stage.iii.ordered.naive = df.degree.stage.iii.naive[rev(order(df.degree.stage.iii.naive$degree)),]
df.degree.stage.ii.ordered.naive = df.degree.stage.ii.naive[rev(order(df.degree.stage.ii.naive$degree)),]
df.degree.stage.i.ordered.top.naive = df.degree.stage.i.ordered.naive[which(df.degree.stage.i.ordered.naive$degree>quantile(df.degree.stage.i.ordered.naive$degree,top.frac)),]
df.degree.stage.iii.ordered.top.naive = df.degree.stage.iii.ordered.naive[which(df.degree.stage.iii.ordered.naive$degree>quantile(df.degree.stage.iii.ordered.naive$degree,top.frac)),]
df.degree.stage.ii.ordered.top.naive = df.degree.stage.ii.ordered.naive[which(df.degree.stage.ii.ordered.naive$degree>quantile(df.degree.stage.ii.ordered.naive$degree,top.frac)),]

# Identify which hubs that are unique to that stage, to mark them in table.
which.unique.stage.i.naive = which(! df.degree.stage.i.ordered.top.naive$protein %in% c(df.degree.stage.iii.ordered.top.naive$protein,df.degree.stage.ii.ordered.top.naive$protein))
which.unique.stage.ii.naive = which(! df.degree.stage.ii.ordered.top.naive$protein %in% c(df.degree.stage.iii.ordered.top.naive$protein,df.degree.stage.i.ordered.top.naive$protein))
which.unique.stage.iii.naive = which(! df.degree.stage.iii.ordered.top.naive$protein %in% c(df.degree.stage.ii.ordered.top.naive$protein,df.degree.stage.i.ordered.top.naive$protein))
which.common.stage.i.naive = which(df.degree.stage.i.ordered.top.naive$protein %in% df.degree.stage.iii.ordered.top.naive$protein & df.degree.stage.i.ordered.top.naive$protein %in% df.degree.stage.ii.ordered.top.naive$protein)
which.common.stage.ii.naive = which(df.degree.stage.ii.ordered.top.naive$protein %in% df.degree.stage.i.ordered.top.naive$protein & df.degree.stage.ii.ordered.top.naive$protein %in% df.degree.stage.iii.ordered.top.naive$protein)
which.common.stage.iii.naive = which(df.degree.stage.iii.ordered.top.naive$protein %in% df.degree.stage.i.ordered.top.naive$protein & df.degree.stage.iii.ordered.top.naive$protein %in% df.degree.stage.ii.ordered.top.naive$protein)

# Maximum table length
dim.neighs.top.naive = max(c(length(df.degree.stage.i.ordered.top.naive$protein),length(df.degree.stage.ii.ordered.top.naive$protein),length(df.degree.stage.iii.ordered.top.naive$protein)))

# Also print gene name
for(i in 1:dim.neighs.top.naive){
  if(is.na(df.degree.stage.i.ordered.top.naive[i,1])) { cat(' & & && ') }
  else {
    if(i %in% which.common.stage.i.naive){
      cat(paste0('\\textbf{',df.degree.stage.i.ordered.top.naive$protein[i], '} & \\emph{',df.degree.stage.i.ordered.top.naive$gene[i], '} & ', df.degree.stage.i.ordered.top.naive$degree[i], ' && '))
    }
    else if(i %in% which.unique.stage.i.naive){
      cat(paste0('\\color{red}{',df.degree.stage.i.ordered.top.naive$protein[i], '} & \\emph{',df.degree.stage.i.ordered.top.naive$gene[i], '} & ', df.degree.stage.i.ordered.top.naive$degree[i], ' && '))
    }
    else { cat(paste0(df.degree.stage.i.ordered.top.naive$protein[i], ' & \\emph{',df.degree.stage.i.ordered.top.naive$gene[i], '} & ', df.degree.stage.i.ordered.top.naive$degree[i], ' && ')) }
  }
  if(is.na(df.degree.stage.ii.ordered.top.naive[i,1])) { cat(' & & && ') }
  else {
    if(i %in% which.common.stage.ii.naive){
      cat(paste0('\\textbf{',df.degree.stage.ii.ordered.top.naive$protein[i], '} & \\emph{',df.degree.stage.ii.ordered.top.naive$gene[i], '} & ', df.degree.stage.ii.ordered.top.naive$degree[i], ' && '))
    }
    else if(i %in% which.unique.stage.ii.naive){
      cat(paste0('\\color{red}{',df.degree.stage.ii.ordered.top.naive$protein[i], '} & \\emph{',df.degree.stage.ii.ordered.top.naive$gene[i], '} & ', df.degree.stage.ii.ordered.top.naive$degree[i], ' && '))
    }
    else { cat(paste0(df.degree.stage.ii.ordered.top.naive$protein[i], ' & \\emph{',df.degree.stage.ii.ordered.top.naive$gene[i], '} & ', df.degree.stage.ii.ordered.top.naive$degree[i], ' && ')) }
  }
  if(is.na(df.degree.stage.iii.ordered.top.naive[i,1])) { cat(' & & \\\\ \n') }
  else {
    if(i %in% which.common.stage.iii.naive){
      cat(paste0('\\textbf{',df.degree.stage.iii.ordered.top.naive$protein[i], '} & \\emph{',df.degree.stage.iii.ordered.top.naive$gene[i], '} & ', df.degree.stage.iii.ordered.top.naive$degree[i], ' \\\\ \n  '))
    }
    else if(i %in% which.unique.stage.iii.naive){
      cat(paste0('\\color{red}{',df.degree.stage.iii.ordered.top.naive$protein[i], '} & \\emph{',df.degree.stage.iii.ordered.top.naive$gene[i], '} & ', df.degree.stage.iii.ordered.top.naive$degree[i], ' \\\\ \n  '))
    }
    else { cat(paste0(df.degree.stage.iii.ordered.top.naive$protein[i], ' & \\emph{',df.degree.stage.iii.ordered.top.naive$gene[i], '} & ', df.degree.stage.iii.ordered.top.naive$degree[i], ' \\\\ \n ')) }
  }
}


# Get edges with gene names
edges.stage.i.gene.naive.all = get_and_print_edges(omega.ssl.i!=0,mapping.frame$gene[match(colnames(brca_dat[[1]]), mapping.frame$protein)])
edges.stage.iii.gene.naive.all = get_and_print_edges(omega.ssl.iii!=0,mapping.frame$gene[match(colnames(brca_dat[[3]]), mapping.frame$protein)])
edges.stage.ii.gene.naive.all = get_and_print_edges(omega.ssl.ii!=0,mapping.frame$gene[match(colnames(brca_dat[[2]]), mapping.frame$protein)])
edges.stage.i.gene.naive = edges.stage.i.gene.naive.all[-which(duplicated(edges.stage.i.gene.naive.all)),]
edges.stage.iii.gene.naive = edges.stage.iii.gene.naive.all[-which(duplicated(edges.stage.iii.gene.naive.all)),]
edges.stage.ii.gene.naive = edges.stage.ii.gene.naive.all[-which(duplicated(edges.stage.ii.gene.naive.all)),]
colnames(edges.stage.i.gene.naive) =  c("Gene1", "Gene2")
colnames(edges.stage.ii.gene.naive) =  c("Gene1", "Gene2")
colnames(edges.stage.iii.gene.naive) =  c("Gene1", "Gene2")

# Check with STRING 
# Stage I
df.all.edges.stage.i.naive =  rbind(edges.stage.i.gene.naive, string.rppa[, 1:2])
ids.in.string.naive = which(duplicated(df.all.edges.stage.i.naive)) - nrow(edges.stage.i.gene.naive) # True the second time an edge occurs (which is the duplicate). We subtract nrow to get the index of the edge in STRING df.
# Count all pairs with evidence
edges.stage.i.gene.evidence.naive = string.rppa[ids.in.string.naive,1:2]
n.extra.pairs.evidence.stage.i.naive = 0
for(i in 1:nrow(edges.stage.i.gene.evidence.naive)){
  n.rep = 0
  for(j in 1:nrow(edges.stage.i.gene.naive.all)){
    if(all(edges.stage.i.gene.naive.all[j,] == edges.stage.i.gene.evidence.naive[i,])){
      n.rep = n.rep + 1
    }
  }
  n.extra.pairs.evidence.stage.i.naive = n.extra.pairs.evidence.stage.i.naive + (n.rep-1)
}
# Number of edges with evidence in STRING database
length(ids.in.string.naive) + n.extra.pairs.evidence.stage.i.naive # 7
(length(ids.in.string.naive)+ n.extra.pairs.evidence.stage.i.naive)/nrow(edges.stage.i.gene.naive.all) # 0.01252236
# Also the amount of unique edges
length(ids.in.string.naive)/nrow(edges.stage.i.gene.naive) # 0.01323251

# Stage II
df.all.edges.stage.ii.naive =  rbind(edges.stage.ii.gene.naive, string.rppa[, 1:2])
ids.in.string.naive = which(duplicated(df.all.edges.stage.ii.naive)) - nrow(edges.stage.ii.gene.naive) # True the second time an edge occurs (which is the duplicate). We subtract nrow to get the index of the edge in STRING df.
# Count all pairs with evidence
edges.stage.ii.gene.evidence.naive = string.rppa[ids.in.string.naive,1:2]
n.extra.pairs.evidence.stage.ii.naive = 0
for(i in 1:nrow(edges.stage.ii.gene.evidence.naive)){
  n.rep = 0
  for(j in 1:nrow(edges.stage.ii.gene.naive.all)){
    if(all(edges.stage.ii.gene.naive.all[j,] == edges.stage.ii.gene.evidence.naive[i,])){
      n.rep = n.rep + 1
    }
  }
  n.extra.pairs.evidence.stage.ii.naive = n.extra.pairs.evidence.stage.ii.naive + (n.rep-1)
}
# Number of edges with evidence in STRING database
length(ids.in.string.naive)+n.extra.pairs.evidence.stage.ii.naive # 23
(length(ids.in.string.naive)+n.extra.pairs.evidence.stage.ii.naive)/nrow(edges.stage.ii.gene.naive.all) # 0.02761104
# Also the amount of unique edges
length(ids.in.string.naive)/nrow(edges.stage.ii.gene.naive) # 0.02393617

# Stage III
df.all.edges.stage.iii.naive =  rbind(edges.stage.iii.gene.naive, string.rppa[, 1:2])
ids.in.string.naive = which(duplicated(df.all.edges.stage.iii.naive)) - nrow(edges.stage.iii.gene.naive) # True the second time an edge occurs (which is the duplicate). We subtract nrow to get the index of the edge in STRING df.
# Count all pairs with evidence
edges.stage.iii.gene.evidence.naive = string.rppa[ids.in.string.naive,1:2]
n.extra.pairs.evidence.stage.iii.naive = 0
for(i in 1:nrow(edges.stage.iii.gene.evidence.naive)){
  n.rep = 0
  for(j in 1:nrow(edges.stage.iii.gene.naive.all)){
    if(all(edges.stage.iii.gene.naive.all[j,] == edges.stage.iii.gene.evidence.naive[i,])){
      n.rep = n.rep + 1
    }
  }
  n.extra.pairs.evidence.stage.iii.naive = n.extra.pairs.evidence.stage.iii.naive + (n.rep-1)
}
# Number of edges with evidence in STRING database
length(ids.in.string.naive)+n.extra.pairs.evidence.stage.iii.naive # 15
(length(ids.in.string.naive)+n.extra.pairs.evidence.stage.iii.naive)/nrow(edges.stage.iii.gene.naive.all) # 0.02772643
# Also the amount of unique edges
length(ids.in.string.naive)/nrow(edges.stage.iii.gene.naive) # 0.02554028


# Have a look at stage-specific edges and their partial correlation ---------------------------------------------------------------------

omega.unique.stage.i = (omega.ssl.x.i!=0 & omega.ssl.x.ii==0 & omega.ssl.x.iii==0)
get_and_print_edges(omega.unique.stage.i!=0,colnames(brca_dat[[1]]), theta.mat = omega.ssl.x.i)

omega.unique.stage.ii = (omega.ssl.x.ii!=0 & omega.ssl.x.i==0 & omega.ssl.x.iii==0)
get_and_print_edges(omega.unique.stage.ii!=0,colnames(brca_dat[[1]]), theta.mat = omega.ssl.x.ii)

omega.unique.stage.iii = (omega.ssl.x.iii!=0 & omega.ssl.x.ii==0 & omega.ssl.x.i==0)
get_and_print_edges(omega.unique.stage.iii!=0,colnames(brca_dat[[1]]), theta.mat = omega.ssl.x.iii)

# Remember these are protein names and not gene names

# Discuss some of these. For example strong negative partial cor (-0.28) between "CAVEOLIN1" and "MEK1" in stage III only: interesting?

# For example edge between "CYCLINE1" and "MAPKPT202Y204 proteins with strong partial correlation (0.25), 
# not much on this edge in literature and no evidence in STRING - a potential new finding? Highlight relevance of both proteins/genes in BRCA and discuss implications

# Print unique Stage III subnetwork -------------------------------------------------------------------

subnet.unique.stage.ii= (omega.ssl.x.iii!=0 & omega.ssl.x.ii== 0 & omega.ssl.x.i==0) + 0
sparsity(subnet.unique.stage.ii) # 0.004110393

set.seed(111)
net.common = network::network(subnet.unique.stage.ii ,directed=F)
net.plot = GGally::ggnet2(net.common, edge.size = 0.3,alpha=0.8,mode = "fruchtermanreingold",color = 'royalblue1',label=T,
                          size='degree',size.min = 1)+
  ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5, size=0.6))


pdf("Genomic_application/plots/net_stageIII_unique.pdf", 12, 12)
net.plot
dev.off()

# Can discuss this a bit - something that can be linked to high disease severity? Relates to the edge mentioned above

# Plot some effects that vary with stage ---------------------------------------------------------

subnet.common= (omega.ssl.x.iii!=0 & omega.ssl.x.ii!= 0 & omega.ssl.x.i!=0) + 0
edges.common.all = get_and_print_edges(subnet.common!=0, mapping.frame$gene[match(colnames(brca_dat[[2]]), mapping.frame$protein)])
edges.common = edges.common.all[-which(duplicated(edges.common.all)),]
colnames(edges.common) =  c("Gene1", "Gene2")

# Which edges are present in all three conditions and have evidence in STRING
df.comb = rbind(edges.common ,string.rppa[, 1:2])
df.comb.merged = apply(df.comb,1,FUN = function(s) paste0(sort(s)[1],sort(s)[2]))
df.evidence.common = df.comb[which(duplicated(df.comb.merged)),]
df.evidence.common
nrow(df.evidence.common) # 17


# Look at edges between proteins with in string
protein1 = mapping.frame$protein[match(df.evidence.common[,1 ],mapping.frame$gene)]
protein2 =mapping.frame$protein[match(df.evidence.common[,2 ],mapping.frame$gene)]
protein1.ind = match(protein1,colnames(brca_dat[[1]]))
protein2.ind = match(protein2,colnames(brca_dat[[1]]))
edges.listed = apply(cbind(protein1,protein2), 1, FUN = function(s) paste(sort(s)[1],'-', sort(s)[2]))
partialcors.i = omega.ssl.x.i[cbind(protein1.ind,protein2.ind)]
partialcors.ii = omega.ssl.x.ii[cbind(protein1.ind,protein2.ind)]
partialcors.iii = omega.ssl.x.iii[cbind(protein1.ind,protein2.ind)]

df.plot = data.frame(partialcor = c(partialcors.i,partialcors.ii,partialcors.iii), 
                     stage = c(rep('I', nrow(df.evidence.common)), rep('II', nrow(df.evidence.common)),rep('III', nrow(df.evidence.common))), 
                     Edge = rep(edges.listed,3))

# only include edges that have a notable change in partial correlation
ind.include = abs(partialcors.i-partialcors.iii)>0.05 
df.plot = df.plot[c(ind.include,ind.include,ind.include),]

gg=ggplot2::ggplot(df.plot,aes(x=stage,y=partialcor, group=Edge)) + geom_line(aes(color=Edge))+geom_point(aes(color=Edge))+theme_bw()+ylab('Partial correlation')+xlab('Tumour stage')+
  geom_hline(yintercept=0, linetype='dashed', color='grey')

pdf("Genomic_application/plots/TCGA_partialcor_change.pdf", 6, 4)
gg
dev.off()


