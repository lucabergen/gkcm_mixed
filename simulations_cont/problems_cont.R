
dgp_null1 <- function(n, error, ...){
  
  eps_X <- rnorm(n, 0, 1)
  eps_Y <- rnorm(n, 0, 1)
  Z <- matrix(rnorm(n*7,0,1),ncol = 7)
  X <- tanh(0.5*Z[,1] - 0.9*Z[,2] + Z[,3] + eps_X)
  Y <- exp(-0.8*Z[,4]*Z[,5] + 0.6*Z[,6]*Z[,7] + eps_Y)
  
  list(X = X, Y = Y, Z = Z)
}

dgp_null2 <- function(n, error, ...){
  
  eps_X <- rnorm(n, 0, 1)
  eps_Y <- rnorm(n, 0, 1)
  Z <- matrix(rnorm(n*7,0,1),ncol = 7)
  X <- sin(2*pi*Z[,1]) + 0.1*eps_X
  Y <- sin(2*pi*Z[,1]) + eps_Y
  
  list(X = X, Y = Y, Z = Z)
}

dgp_alt1 <- function(n, error, ...){
  
  eps_X <- rnorm(n, 0, 1)
  eps_Y <- rnorm(n, 0, 1)
  Z <- matrix(rnorm(n*7,0,1),ncol = 7)
  Y <- exp(-0.8*Z[,4]*Z[,5] + 0.6*Z[,6]*Z[,7] + eps_Y)
  X <- tanh(0.5*Z[,1] - 0.9*Z[,2] + 0.3*Z[,3]*Y + eps_X)
  
  list(X = X, Y = Y, Z = Z)
}

dgp_alt2 <- function(n, error, ...){
  
  eps_X <- rnorm(n, 0, 1)
  eps_Y <- rnorm(n, 0, 1)
  Z <- matrix(rnorm(n*7,0,1),ncol = 7)
  Y <- 0.2*Z[,2]^3 + tanh(Z[,4]) + eps_Y
  X <- sin(Z[,1]) - 0.4*Z[,2]^2 + cos(0.2*pi*Y) + eps_X

  list(X = X, Y = Y, Z = Z)
}


