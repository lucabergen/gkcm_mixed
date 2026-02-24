
GKCM_wrapper <- function(data, job, instance, ...){
  
  gaussian_kernel <- function(x, y = NULL, sigma = 1) {
    
    if(!is.vector(x)) x <- as.vector(x)
    
    if (is.null(y)) y <- x
    
    exp(-outer(x, y, "-")^2 / (2 * sigma^2))
    
  }
  
  center_kernel <- function(K) {
    
    n <- nrow(K)
    H <- diag(n) - matrix(1/n, n, n)
    
    H %*% K %*% H
    
  }
  
  compute_residuals <- function(K, W) {
    
    n <- nrow(K)
    
    # Compute B = (I - W)
    B <- diag(n) - W
    R <- B %*% K %*% t(B)
    
    # Symmetrize
    (R + t(R)) / 2
    
  }
  
  GKCM <- function(X, Y, Z, eps = 1e-16, seed = NULL) {
    
    n <- nrow(Z)
    if (!is.null(seed)) {set.seed(seed)}
    
    X <- scale(X)
    Y <- scale(Y)
    Z <- scale(Z)
    
    K <- gaussian_kernel(X, sigma = 1)
    L <- gaussian_kernel(Y, sigma = 1)
    
    # Set honesty = F to avoid sample splitting
    DRF_X <- drf::drf(Z, X, num.trees = ncol(Z)*100, bandwidth = 1, 
                      honesty = F, response.scaling = F, 
                      min.node.size = 5, mtry = 40)
    
    DRF_Y <- drf::drf(Z, Y, num.trees = ncol(Z)*100, bandwidth = 1, 
                      honesty = F, response.scaling = F, 
                      min.node.size = 5, mtry = 40)
    
    W_X <- drf::get_sample_weights(DRF_X) |> as.matrix()
    W_Y <- drf::get_sample_weights(DRF_Y) |> as.matrix()
    
    K_XZ <- compute_residuals(K, W_X) |> center_kernel()
    L_YZ <- compute_residuals(L, W_Y) |> center_kernel()
    
    # Compute test statistic
    t_stat <- sum(K_XZ * L_YZ) / n
    
    G <- K_XZ * L_YZ
    H <- center_kernel(G)/(n - 1)
    
    ev <- eigen(H, symmetric = TRUE)$values
    ev <- ev[ev > eps]
    
    # Moment based approximation of the null  
    p_val <- tryCatch({
      1 - momentchi2::lpb4(ev, t_stat)
    }, error = function(e) {
      1 - momentchi2::hbe(ev, t_stat)
    })
    
    list(t_stat = t_stat, p_val = p_val) 
    
  }
  
  GKCM(
    X = instance$X,
    Y = instance$Y,
    Z = instance$Z
  )$p_val
  
}

PCM_wrapper <- function(data, job, instance, ...){
  
  comets::pcm(
    X = instance$X,
    Y = instance$Y,
    Z = instance$Z,
    reg_YonXZ = "tuned_rf",
    reg_YonZ = "tuned_rf",
    reg_YhatonZ = "tuned_rf",
    reg_VonXZ = "tuned_rf",
    reg_RonZ = "tuned_rf"
  )$p.value
  
}

wGCM_wrapper <- function(data, job, instance, ...){
  
  comets::wgcm(
    X = instance$X,
    Y = instance$Y,
    Z = instance$Z,
    reg_YonZ = "tuned_rf", 
    reg_XonZ = "tuned_rf", 
    reg_wfun = "tuned_rf"
  )$p.value[[1]]
  
}

GCM_wrapper <- function(data, job, instance, ...){
  
  comets::gcm(
    X = instance$X,
    Y = instance$Y,
    Z = instance$Z,
    reg_YonZ = "tuned_rf",
    reg_XonZ = "tuned_rf"
  )$p.value
  
}

KCIT_wrapper <- function(data, job, instance, ...){

  RCIT::KCIT(
    x = instance$X,
    y = instance$Y,
    z = instance$Z)
}

RCIT_wrapper <- function(data, job, instance, ...){

  RCIT::RCIT(
    x = instance$X,
    y = instance$Y,
    z = instance$Z)$p
}

RCoT_wrapper <- function(data, job, instance, ...){
  
  RCIT::RCoT(
    x = instance$X,
    y = instance$Y,
    z = instance$Z)$p
}

