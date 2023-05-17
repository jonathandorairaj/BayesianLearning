
######Q1######


library(ggplot2)
library(mvtnorm)
library(gridExtra)
#library(extraDistr)
precipitation <- readRDS("Precipitation.rds")
log_prec <- log(precipitation)

mu_post <- function(mu_0,tau_0_sqr,sigma_square_current,y,n){
  
  tau_n_sqr = 1 / (1/tau_0_sqr + n/sigma_square_current)
  mu_n <- tau_n_sqr * (mu_0/tau_0_sqr + sum(y)/sigma_square_current)
  mupost <- rnorm(1, mu_n, sqrt(tau_n_sqr))
  return(mupost)
  
}

sigma_square_post <- function(mu_current, v_0, sigma_square_0, y,n, use_myinvchi=TRUE) {
  
  v_n <- v_0 + n
  elem1 <- v_0*sigma_square_0
  elem2 <- sum((y - mu_current)^2)
  elem3 <- n+v_0
  elem_comb <- (elem1+elem2)/elem3
  if (use_myinvchi){
    sigma_square_post <- my_inv_chi(v_n,elem_comb)
  }
  else {
    sigma_square_post <- rinvchisq(1, v_n, elem_comb) #Requires package extraDistr
  }
  
  
  
  return(sigma_square_post)
}

my_inv_chi<- function(df,tau_sqr) {
  X <- rchisq(1,df)
  inv_chi <- (df*tau_sqr)/X
  return(inv_chi)
}

#init
mu_0 <- 0
tau_0_sqr <- 1
sigma_square_0 <- 1 
v_0 <- 1 #degree of freedom for chi square

sigma_square_0 <- 1
v_0 <- 1 #degree of freedom for chi square
gibbs_sampler <- function(nstep, data, mu_0, tau_0_sqr, v_0, sigma_square_0) {
  # Init parameters
  mu_current <- 0
  sigma_square_current <- 1
  mu_samples <- rep(0,nstep)
  sigma_square_samples <- rep(0,nstep)
  for (i in 1:nstep) {
    mu_current <- mu_post(mu_0, tau_0_sqr, sigma_square_current, y=data, length(data))
    #print(mu_current)
    sigma_square_current <- sigma_square_post(mu_current, v_0, sigma_square_0, y=data,
                                              length(data),use_myinvchi=TRUE)
    #print(sigma_square_current)
    mu_samples[i] <- mu_current
    sigma_square_samples[i] <- sigma_square_current
  }
  output_df <- data.frame(mu_sample = mu_samples, sigma_sample = sigma_square_samples)
  return(output_df)
}
sample_gibbs <- gibbs_sampler(nstep=10000, data=log_prec, mu_0, tau_0_sqr, v_0, sigma_square_0)
my_acf <- acf(sample_gibbs$mu_sample,plot = F) # getting the autocorrelation
if_mu <- 1 + 2 * sum(my_acf$acf[-1])
my_acf <- acf(sample_gibbs$sigma_sample,plot = F)
if_sigma <- 1 + 2 * sum(my_acf$acf[-1])
ggplot(data=sample_gibbs, aes(x = 1:length(mu_sample), y = mu_sample)) +
  geom_line()+labs(title = "mu sample")
ggplot(data=sample_gibbs, aes(x = 1:length(sigma_sample), y = sigma_sample)) +
  geom_line() +labs(title = "sigma sample")
# Using Gibbs sample's mu and sigma to sample posterior prediction
post_pred_samples <- rnorm(n = 10000, mean = sample_gibbs$mu_sample,
                           sd = sqrt(sample_gibbs$sigma_sample))
ggplot(data = data.frame(y = log_prec), aes(x = y)) +
  geom_histogram(aes(y =after_stat(density)), color = "black", binwidth = 0.2) +
  geom_density(data = data.frame(y = post_pred_samples), aes(x = y, y = after_stat(density)),
               color = "red",linewidth = 2) +
  labs(x = "Log- Daily Precipitation ", y = "Density")

######Q2######

ebay_data <- read.table("eBayNumberOfBidderData.dat", header = TRUE)
model1 <- glm(formula = nBids ~ .,data = ebay_data[,-2],family = 'poisson')
summary(model1)
# Initialize values
n_cols <- ncol(ebay_data[,-1])
#remove 1st column since that is target variable and convert to matrix
# matrix of features
covariates <- as.matrix(ebay_data[,-1])
labels <- as.matrix(ebay_data[,1])
mu <- rep(0, n_cols)
initVal <- matrix(0, n_cols, 1)
Sigma <- as.matrix(100 * solve(t(covariates)%*%covariates))
LogPosteriorFunc <- function(betas, X, y, mu, Sigma){
  log_prior <- dmvnorm(betas, mu, Sigma, log=TRUE)
  log_likelihood <- sum(X%*%betas * y - exp(X%*%betas) -log(factorial(y)))
  res <- log_prior + log_likelihood
  return(res)
}
# Optimizer
OptimRes <- optim(initVal, LogPosteriorFunc, gr = NULL, y = labels, X = covariates,
                  mu = mu, Sigma = Sigma, method=c("BFGS"),
                  control=list(fnscale=-1), hessian=TRUE)
beta_mode <- OptimRes$par
jacobian <- OptimRes$hessian
inv_jacobian <- -solve(jacobian)
beta_draws <- as.matrix(rmvnorm(10000,mean = beta_mode,sigma = inv_jacobian))
beta_estimate <- colMeans(beta_draws)
hist(beta_draws,breaks = 50,main = 'Histogram of Posterior Draws',xlab = 'Betas')
MetHas_RandomWalk <- function(nDraws,fun,mu,Sigma,c){
  #initialize matrix
  draw_matrix <- matrix(0,nrow = nDraws,ncol = n_cols)
  #initialize first row to mu
  draw_matrix[1,] <- mu
  for(i in 2:nDraws){
    # sample from multivariate normal distribution
    proposed_sample <- as.vector(rmvnorm(n = 1,mean = draw_matrix[i-1,],
                                         sigma = c*as.matrix(Sigma)))
    #print(proposed_sample)
    # IMPORTANT : the log is inside the posterior function
    log_acceptance_prob <- exp(fun(proposed_sample)- fun(draw_matrix[i-1,]))
    #random sample
    u <- runif(1)
    # calculate acceptance probability
    a <- min(1,log_acceptance_prob)
    if(u <= a){
      #accept sample
      draw_matrix[i,] <- proposed_sample
    }
    else{
      # stay at same values from previous draw
      draw_matrix[i,] <- draw_matrix[i-1,]
    }
  }
  return(draw_matrix)
}
# this function we pass to MetHas Algorithm, can be changed to another posterior density
logPostFunc <- function(theta){
  res <- dmvnorm(theta,mean = beta_estimate,sigma = inv_jacobian,log = TRUE)
  if(is.na(res)){
    print(theta)
  }
  return(res)
}
df <- MetHas_RandomWalk(nDraws = 10000,fun = logPostFunc,mu = rep(0,n_cols),
                        Sigma = inv_jacobian,c = 1)
# assign colnames
colnames(df) <- colnames(ebay_data)[2:10]
## plotting
plot_list <- list()
for (col in colnames(df)) {
  # Plot iterations vs every column
  p <- ggplot(data = as.data.frame(df), aes_string(x = 1:nrow(df), y = col)) +
    geom_line(col = 'blue') +
    labs(x = "Iterations", y = col)
  #show plot
  #print(p)
  plot_list[[col]] <- p
}
# Arranging in 1 fig
grid.arrange(grobs = plot_list, ncol = 2)
# extra 1 at the start for the intercept - Const
new_data <- c(1,1,0,1,0,1,0,1.2,0.8)
# lambda = eˆBeta*x
# discarding first 1500 samples as burn-in period
lambda <- exp(df[-c(1:1500),] %*% new_data)
samples <- c()
for (i in 1:nrow(df[-c(1:1500),])) {
  #sample from each row of df to get the predictive distribution based on the
  #posterior betas
  samples[i] <- rpois(1,lambda = lambda[i])
}
hist(samples)
res <- length(samples[samples == 0])/length(samples)


######Q3######

mu <- 13
sigma_square <- 3
t <- 300
phi <- seq(from = -1 ,to = 1, by =0.25 )

ar_func <- function(phi,mu,sigma_square,t){
  counter <- 0
  result_vector <- rep(0,t)
  result_vector[1] <- mu
  for(i in 2:t){
    epislon <- rnorm(1,0,sqrt(sigma_square))
    x_i <- mu+phi*(result_vector[i-1]-mu)+ epislon
    result_vector[i] <- x_i
  }
  return(result_vector)
}
test_phi_func <- function(phi,mu,sigma_square,t){
  phi_test_df <- data.frame(matrix(0, nrow = t, ncol = length(phi)))
  colnames(phi_test_df) <- phi
  for (j in 1:length(phi)) {
    phi_test <- ar_func(phi[j], mu, sigma_square, t)
    phi_test_df[, j] <- phi_test
  }
  return(phi_test_df)
}
phi_df <- test_phi_func(phi,mu,sigma_square,t)
for(k in 1:length(phi)){
  plot_data <- phi_df[,k]
  plot(x=1:300, plot_data, type = "l", main = paste(" phi = ", phi[k]))
  #Sys.sleep(1)
}
library(rstan)
set.seed(12345)
x <- ar_func(phi=0.2, mu, sigma_square, t)
y <- ar_func(phi=0.95, mu, sigma_square, t)
stan_code <- "
data {
int<lower=0> T; // Number of time points
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
sigma_x ~ normal(0, 50);
sigma_y ~ normal(0, 50);
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
niter <- 5000
warmup <- 500
# Compile the Stan model
model <- stan_model(model_code = stan_code)
# Fit the Stan model to the data
# Set the control to avoid too many divergent after warmup
control <- list(adapt_delta = 0.90, stepsize = 0.0001)
fit<- sampling(model, data = data_list, warmup = warmup, iter = niter, chains = 5,control=control,refresh = 0)
#refresh = 0 to mute the printout
             
summary(fit)$summary
post_mean <- summary(fit)$summary[, "mean"]
interval_025 <- summary(fit)$summary[, "25%"]
interval_975 <- summary(fit)$summary[, "97.5%"]
n_eff <- summary(fit)$summary[, "n_eff"]
Rhat <- summary(fit)$summary[, "Rhat"]
post_mean
interval_025
interval_975
n_eff
Rhat

posterior_x <- extract(fit, pars = c("mu_x", "phi_x"))
posterior_df_x <- data.frame(posterior_x)
ggplot(data = posterior_df_x, aes(x = mu_x, y = phi_x)) +
  stat_density_2d() +
  xlab("mu_x") +
  ylab("phi_x") +
  ggtitle("Join tposterior of mu_x and phi_x")

posterior_y <- extract(fit, pars = c("mu_y", "phi_y"))
posterior_df_y <- data.frame(posterior_y)
ggplot(data = posterior_df_y, aes(x = mu_y, y = phi_y)) +
  stat_density_2d()+
  xlab("mu_y") +
  ylab("phi_y") +
  ggtitle("Joint posterior of mu_y and phi_y")