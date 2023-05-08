rm(list=ls())
library(ggplot2)
library(mvtnorm)
#library(extraDistr)

precipitation <- readRDS("Precipitation.rds")
log_prec <- log(precipitation)

mu_post <- function(mu_0,tau_0_sqr,sigma_sqr_current,y,n){
  
  tau_n_sqr = 1 / (1/tau_0_sqr + n/sigma_sqr_current)
  mu_n <- tau_n_sqr * (mu_0/tau_0_sqr + sum(y)/sigma_sqr_current)
  mupost <- rnorm(1, mu_n, sqrt(tau_n_sqr))
  return(mupost)
  
}

sigma_sqr_post <- function(mu_current, v_0, sigma_sqr_0, y,n, use_myinvchi=TRUE) {
  
  v_n <- v_0 + n
  elem1 <- v_0*sigma_sqr_0
  elem2 <- sum((y - mu_current)^2)
  elem3 <- n+v_0
  elem_comb <- (elem1+elem2)/elem3
  if (use_myinvchi){
    sigma_sqr_post <- my_inv_chi(v_n,elem_comb)
  }
  else {
    sigma_sqr_post <- rinvchisq(1, v_n, elem_comb) #Requires package extraDistr
  }
  
  

  return(sigma_sqr_post)
}

my_inv_chi<- function(df,tau_sqr) {
  X <- rchisq(1,df)
  inv_chi <- (df*tau_sqr)/X
  return(inv_chi)
}

#init
mu_0 <- 0
tau_0_sqr <- 1
sigma_sqr_0 <- 1 # These are also not the same
v_0 <- 1 #degree of freedom for chi square

gibbs_sampler <- function(nstep, data, mu_0, tau_0_sqr, v_0, sigma_sqr_0) {
  # Initialize the parameters
  mu_current <- 0
  sigma_sqr_current <- 1

  mu_samples <- rep(0,nstep)
  sigma_sqr_samples <- rep(0,nstep)
  
  for (i in 1:nstep) {
   
    mu_current <- mu_post(mu_0, tau_0_sqr, sigma_sqr_current, y=data, length(data))
    #print(mu_current)
    sigma_sqr_current <- sigma_sqr_post(mu_current, v_0, sigma_sqr_0, y=data, length(data),use_myinvchi=TRUE)
    #print(sigma_sqr_current)

    mu_samples[i] <- mu_current
    sigma_sqr_samples[i] <- sigma_sqr_current
  }
  
  output_df <- data.frame(mu_sample = mu_samples, sigma_sample = sigma_sqr_samples)
  
  # Return the posterior samples
  return(output_df)
}
sample_gibbs <- gibbs_sampler(nstep=10000, data=log_prec, mu_0, tau_0_sqr, v_0, sigma_sqr_0)
ggplot(data=sample_gibbs, aes(x = 1:length(mu_sample), y = mu_sample)) +
  geom_line() 

ggplot(data=sample_gibbs, aes(x = 1:length(sigma_sample), y = sigma_sample)) +
  geom_line() 


#in lecture slide IF=1 + 2 ∑^∞_{k=1} ρ_k where  ρ_k is autocorrelation at lag k
my_acf <- acf(sample_gibbs)
if_mu <- 1 + 2*sum(my_acf)
if_sigma <- 1 + 2*sum(acf[,2])

