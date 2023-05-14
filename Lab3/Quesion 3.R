rm(list=ls())
library(ggplot2)
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


#####2nd exercise#####
library(rstan)
x <- ar_func(phi=0.2, mu, sigma_sqr, t)
y <- ar_func(phi=0.95, mu, sigma_sqr, t)

stan_code <- "
data {
  int<lower=0> T;         // Number of time points
  vector[T] x;            
  vector[T] y;            
}
parameters {
  real mu_x;
  real mu_y;
  real phi_x;
  real phi_y;
  real sigma_x;
  real sigma_y;
}
model {
  // After some research, it is common to use a flat prior or a vague prior as 
  // non-informative prior
  // Therefore, we pick normal distribution with higher variance as prior.
  mu_x ~ normal(0, 50);
  mu_y ~ normal(0, 50);
  phi_x ~ normal(0, 10);
  phi_y ~ normal(0, 10);
  
  sigma_x  ~ normal(0, 50);
  sigma_y  ~ normal(0, 50);
  

  x[2:T] ~ normal(mu_x + phi_x * (x[1:(T - 1)] - mu_x), sigma_x);
  y[2:T] ~ normal(mu_y + phi_y * (y[1:(T - 1)] - mu_y), sigma_y);
}
"



data_list <- list(
  T = t,
  x = x,
  y = y
)

# Set the MCMC settings
niter <- 2500
warmup <- 500

# Compile the Stan model
model <- stan_model(model_code = stan_code)

# Fit the Stan model to the data
# Set the conture to avoid too many divergent after warmup
control <- list(adapt_delta = 0.90, stepsize = 0.0001)
fit<- sampling(model, data = data_list, warmup = warmup, iter = niter, chains = 6,control=control)
summary(fit)$summary


posterior_x <- extract(fit, pars = c("mu_x", "phi_x"))
posterior_df_x <- data.frame(posterior_x)
ggplot(data = posterior_df_x, aes(x = mu_x, y = phi_x)) +
  stat_density_2d() +
  xlab("mu_x") +
  ylab("phi_x") +
  ggtitle("Joint distribution of mu_x and phi_x")


posterior_y <- extract(fit, pars = c("mu_y", "phi_y"))
posterior_df_y <- data.frame(posterior_y)
ggplot(data = posterior_df_y, aes(x = mu_y, y = phi_y)) +
  stat_density_2d()+
  xlab("mu_y") +
  ylab("phi_y") +
  ggtitle("Joint distribution of mu_y and phi_y")
