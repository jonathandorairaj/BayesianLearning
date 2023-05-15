
# 2a) 
library(mvtnorm)
library(ggplot2)

ebay_data <- read.table("eBayNumberOfBidderData.dat", header = TRUE)

model <- glm(formula = nBids ~ .,data = ebay_data,family = 'poisson')

summary(model)

#2b)

LogPosteriorFunc <- function(beta, X, y, mu, Sigma){
  log_prior <- dmvnorm(beta, mu, Sigma, log=TRUE)
  log_likelihood <- sum(X%*%beta * y - exp(X%*%beta) -log(factorial(y)))
  return(log_likelihood + log_prior)
}

# Initialize values
n_cols <- ncol(ebay_data[,-1])
covariates <- as.matrix(ebay_data[,-1])
labels <- as.matrix(ebay_data[,1])
mu <- rep(0, n_cols)
initVal <- matrix(0, n_cols, 1)
Sigma <- as.matrix(100 * solve(t(covariates)%*%covariates))

# Optimizer
OptimRes <- optim(initVal, LogPosteriorFunc, gr = NULL, y = labels, X = covariates, mu = mu, Sigma = Sigma, method=c("BFGS"), control=list(fnscale=-1), hessian=TRUE)

beta_mode <- OptimRes$par
jacobian <- OptimRes$hessian
inv_jacobian <- -solve(jacobian)

beta_draws <- as.matrix(rmvnorm(10000,mean = beta_mode,sigma = inv_jacobian))
beta_estimates <- colMeans(beta_draws)


### 

MetHas_RandomWalk <- function(nDraws,fun,mu,Sigma,c){

  draw_matrix <- matrix(0,nrow = nDraws,ncol = n_cols)
  #initialize first row to mu 
  draw_matrix[1,] <- mu
  
  for(i in 2:nDraws){
    #
    proposed_sample <- as.vector(rmvnorm(n = 1,mean = draw_matrix[i-1,],sigma = c*as.matrix(Sigma)))
    print(proposed_sample)
    # the log is inside the posterior function
    log_acceptance_prob <- exp(fun(proposed_sample)- fun(draw_matrix[i-1,]))
    
    u <- runif(1)
    a <- min(1,log_acceptance_prob)
    if(u <= a){
      draw_matrix[i,] <- proposed_sample
    }
    else{
      draw_matrix[i,] <- draw_matrix[i-1,]
    }
    
  }
  return(draw_matrix)
}
  
logPostFunc <- function(theta){
  res <- dmvnorm(theta,mean = beta_estimates,sigma = inv_jacobian,log = TRUE)
  if(is.na(res)){
    print(theta)
  }
  return(res)
}

df <- MetHas_RandomWalk(nDraws = 10000,fun = logPostFunc,mu = rep(0,n_cols),Sigma = inv_jacobian,c = 1)

colnames(df) <- colnames(ebay_data)[2:10]
## plotting 
plot_list <- list()
for (col in colnames(df)) {
  # Create a new plot for each column
  p <- ggplot(data = as.data.frame(df), aes_string(x = 1:nrow(df), y = col)) +
    geom_line(col = 'blue') +
    labs(x = "Iterations", y = col)  # Add axis labels
  
  # Display the plot
  print(p)
  plot_list[[col]] <- p
}

p


library(gridExtra)

# Arrange the four plots in a 3x3 grid
grid.arrange(grobs = plot_list, ncol = 3)


new_data <- c(1,1,0,1,0,1,0,1.2,0.8)
lambda <- exp(df[-c(1:1000),] %*% new_data)

samples <- c()
for (i in 1:nrow(df[-c(1:1000),])) {
  samples[i] <- rpois(1,lambda = lambda)
  
}

hist(samples)

head(samples)

length(samples[samples == 0])/length(samples)
