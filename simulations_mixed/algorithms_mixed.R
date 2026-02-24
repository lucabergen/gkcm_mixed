
GKCM_wrapper <- function(data, job, instance, ...){
  
  gaussian_kernel <- function(x, y = NULL, sigma = 1) {
    
    if(!is.vector(x)) x <- as.vector(x)
    
    if (length(sigma) != 1) {
      stop("For gaussian_kernel, sigma must be a scalar.")
    }
    
    if (is.null(y)) y <- x
    
    K <- exp(-outer(x, y, "-")^2 / (2 * sigma^2))
    
    K
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
  
  GKCM_mixed <- function(X, Y, Z, eps = 1e-16, seed = NULL) {
    
    n <- nrow(Z)
    if (!is.null(seed)) {set.seed(seed)}
    
    Y <- scale(Y)
    Z <- scale(Z)
    
    # Continuous Y
    L <- gaussian_kernel(Y, sigma = 1)
    
    DRF <- drf::drf(Z, Y, num.trees = ncol(Z)*100, bandwidth = 1, 
                    honesty = F, response.scaling = F, 
                    min.node.size = 5, mtry = 40)
    
    # Use newdata = Z as additional argument to get in-sample (non-OOB) weights 
    W <- drf::get_sample_weights(DRF) |> as.matrix()
    
    L_YZ <- compute_residuals(L, W) |> center_kernel()
    
    
    # Discrete X
    k <- length(unique(X))
    
    rf <- ranger::ranger(formula  = X ~ ., data = data.frame(X = X, Z = Z), 
                         mtry = 5, num.trees = ncol(Z)*50, 
                         probability = TRUE, splitrule = "gini")
    
    pred <- predict(rf, data = data.frame(Z = Z), type = "response")
    P <- pred$predictions          # N x K matrix
    
    # Ensure label levels match columns of P
    classes <- colnames(P)
    X <- factor(X, levels = classes)
    
    # Generate one hot encondings
    OH <- matrix(0, nrow = n, ncol = k)
    colnames(OH) <- classes
    OH[cbind(seq_len(n), as.integer(X))] <- 1
    
    # Probability residuals
    R <- OH - P
    
    # Inner products
    K_XZ <- (R %*% t(R)) |> center_kernel()
    
    # Compute test statistic
    t_stat <- sum(K_XZ * L_YZ) / n
    
    G <- K_XZ * L_YZ
    H <- center_kernel(G)/(n - 1)
    
    ev <- eigen(H, symmetric = TRUE)$values
    
    # Moment based approximation of the null  
    p_val <- tryCatch({
      1 - momentchi2::lpb4(pmax(ev, eps), t_stat)
    }, error = function(e) {
      1 - momentchi2::hbe(pmax(ev, eps), t_stat)
    })
    
    list(t_stat = t_stat, p_val = p_val) 
    
  }
  
  GKCM_mixed(
    X = instance$X,
    Y = instance$Y,
    Z = instance$Z
  )$p_val
  
}

PCM_wrapper <- function(data, job, instance, ...){
  
  df <- data.frame(
    X = instance$X, Y = instance$Y, 
    setNames(as.data.frame(instance$Z), paste0("Z",1:7))
  )
  
  comets::comet(
    Y ~ X | Z1 + Z2 + Z3 + Z4 + Z5 + Z6 + Z7,
    data = df,
    test = "pcm", 
    reg_YonZ = "tuned_rf", reg_YhatonZ = "tuned_rf",
    reg_YonXZ = "tuned_rf", reg_VonXZ = "tuned_rf", 
    reg_RonZ = "tuned_rf")$p.value
  
}

wGCM_wrapper <- function(data, job, instance, ...){
  
  df <- data.frame(
    X = instance$X, Y = instance$Y, 
    setNames(as.data.frame(instance$Z), paste0("Z",1:7))
  )
  
  comets::comet(
    Y ~ X | Z1 + Z2 + Z3 + Z4 + Z5 + Z6 + Z7,
    data = df,
    test = "wgcm", 
    reg_YonZ = "tuned_rf", 
    reg_XonZ = "tuned_rf", 
    reg_wfun = "tuned_rf"
  )$p.value[[1]]
  
}

GCM_wrapper <- function(data, job, instance, ...){
  
  df <- data.frame(
    X = instance$X, Y = instance$Y, 
    setNames(as.data.frame(instance$Z), paste0("Z",1:7))
  )
  
  comets::comet(
    Y ~ X | Z1 + Z2 + Z3 + Z4 + Z5 + Z6 + Z7,
    data = df,
    test = "gcm",
    reg_YonZ = "tuned_rf",
    reg_XonZ = "tuned_rf"
  )$p.value
  
}
