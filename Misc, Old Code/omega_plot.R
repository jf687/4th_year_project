library(ggplot2)


#x_values <- seq(-3, 3, length.out = 100)

# Calculate the corresponding probability density function (pdf) values for the normal distribution
#y_values <- dnorm(x_values, 0, 0.5)

# Plot the normal distribution
#plot(x_values, y_values, type = "l", lty = 1, col = "blue", xlab = "X", ylab = "Density",
     #main = "Normal Distribution")
plot_omega <- function(){

  sparsity <- 0.35
  M <- 20000
  spikes <- M*(1-sparsity)
  slabs <- M*sparsity
  values <- c()
  dist <- c(rep('sp', spikes), rep('sl',slabs))
  
  Omega <- c()
  v1 <- 1
  v0 <- 10
  for (n in 1:spikes){
    omega <- rnorm(1,0,v1)
    Omega <- c(Omega, omega)
  }
  for (n in 1:slabs){
    omega <- rnorm(1,0,v0)
    Omega <- c(Omega, omega)
  }
  
  df <- data.frame(Omega, dist)
  ggplot(df, aes(x = Omega)) + geom_histogram(aes(y=..density..), bins = 60, color = 'blue',fill = "blue", color = "black", alpha = 0.7) + labs(title = "Spike-and-Slab Prior Distribution",
         x = "Precision Matrix Entry Values",
         y = "Density")+
    geom_density(alpha=0.1, color = 'black', fill="blue", adjust = 3, lty = 7) 
}

