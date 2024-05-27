generate_huge_data <- function(){
  n <- 200
  p <- 100
  sf200_100 = huge::huge.generator(n=, d=p,graph = 'scale-free', prob = 0.02)
  p<- 150
  sf200_150 = huge::huge.generator(n=, d=p,graph = 'scale-free', prob = 0.02)
  p <- 200
  sf200_200 = huge::huge.generator(n=, d=p,graph = 'scale-free', prob = 0.02)
  save(sf200_100, file = 'HUGE_data/one.RData')
  save(sf200_150, file = 'HUGE_data/two.RData')
  save(sf200_200, file = 'HUGE_data/three.RData')
}