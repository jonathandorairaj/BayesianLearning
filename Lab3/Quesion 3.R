library(ggplot2)
rm(list=ls())
mu <- 13
sigma_sqr <- 3
t <- 300
phi <- seq(from = -1 ,to = 1, by =0.25 )

ar_func <- function(phi,mu,sigma_sqr,t){
  counter <- 0
  result_vector <- rep(0,t)
  result_vector[1] <- mu 
  for(i in 2:t){
    epislon <- rnorm(1,0,sqrt(sigma_sqr))
    x_i <- mu+phi*(result_vector[i-1]-mu)+ epislon
    result_vector[i] <- x_i 
  }
  
  return(result_vector)
}


test_phi_func <- function(phi,mu,sigma_sqr,t){
  
  phi_test_df <- data.frame(matrix(0, nrow = t, ncol = length(phi)))
  colnames(phi_test_df) <- phi
  
  for (j in 1:length(phi)) {
    phi_test <- ar_func(phi[j], mu, sigma_sqr, t)
    phi_test_df[, j] <- phi_test
  }
  return(phi_test_df)
}
phi_df <- test_phi_func(phi,mu,sigma_sqr,t)


for(k in 1:length(phi)){
  plot_data <- phi_df[,k]
  plot(x=1:300, plot_data, type = "l", main = paste(" phi = ", phi[k]))
  Sys.sleep(1)
}

