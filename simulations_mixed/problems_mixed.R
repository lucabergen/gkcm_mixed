
dgp_null1 <- function(n, ...){
  
  Z <- matrix(rnorm(n*7,0,1),ncol = 7)
  eps <- rnorm(n,0,1)
  
  g1 <- function(Z){0.6*sin(Z[,1]*Z[,2])}
  g2 <- function(Z){rep(0, nrow(Z))}
  g3 <- function(Z){- 0.7*Z[,3]^2}
  
  Y <- 0.8*sin(Z[,1]*Z[,2]) - 0.5*Z[,3]^2 + 0.5*eps
  
  ex <- exp(cbind(g1(Z), g2(Z), g3(Z)))
  prob <- 0.1 + 0.7 * sweep(ex, 1, rowSums(ex), "/")   
  X <- apply(prob, 1, function(p) sample(1:3, size = 1, prob = p))
  X <- factor(X, levels = 1:3)

  list(X = X, Y = Y, Z = Z)
  
}

dgp_alt1 <- function(n, error, ...){
  
    Z <- matrix(rnorm(n*7,0,1),ncol = 7)
    eps <- rnorm(n,0,1)

    g1 <- function(Z){0.5*Z[,2] - 0.4*Z[,4]^2}
    g2 <- function(Z){0.6*exp(Z[,1]) + 0.3*Z[,2]*Z[,3]}
    g3 <- function(Z,Y){0.7*Z[,5] + 0.8*(sin(Y + Z[,6]))}

    Y <- 0.6*Z[,1] - 0.3*Z[,4]*Z[,5] + 0.2*Z[,6]^2 + eps

    ex <- exp(cbind(g1(Z), g2(Z), g3(Z,Y)))
    prob <- 0.1 + 0.7 * sweep(ex, 1, rowSums(ex), "/")
    X <- apply(prob, 1, function(p) sample(1:3, size = 1, prob = p))
    X <- factor(X, levels = 1:3)

    list(X = X, Y = Y, Z = Z)
    
}
