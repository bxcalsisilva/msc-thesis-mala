# File: R/data_generation_power.R

#' Universal Power Logit Data Generator (Deterministic)
#' 
#' @param n Integer. Number of observations.
#' @param betas Numeric vector. True coefficients (including intercept).
#' @param lambda Numeric. Skewness parameter (lambda > 0). lambda=1 is standard logit.
#' @param rho Numeric. AR(1) correlation coefficient.
#' @param seed Integer. Seed for exact reproducibility of this function.
#' @param scale_X Boolean. If TRUE, standardizes covariates (Recommended).
#' 
#' @return A list containing X (Matrix), y (Vector), and metadata.

simulate_power_logit_data <- function(n, betas, lambda = 1, rho = 0, seed = NULL, scale_X = TRUE) {
  
  # 0. Ensure covariates have names
  if (is.null(names(betas))) {
    names(betas) <- c("(Intercept)", paste0("x", 1:(length(betas)-1)))
  }
  
  # 1. Seed Control (Monte Carlo Safety)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  p <- length(betas) - 1
  
  # 2. Generate Covariates
  if (rho == 0) {
    # --- INDEPENDENT CASE ---
    # Generate a matrix of p columns with standard normal noise N(0,1).
    X_raw <- matrix(rnorm(n * p), nrow = n, ncol = p)
    
  } else {
    Sigma <- matrix(0, nrow = p, ncol = p)
    
    # --- CORRELATED CASE (Auto-Regressive AR1 Structure) ---
    # A. Create the distance matrix
    exponent <- abs(outer(1:p, 1:p, "-"))
    
    # B. Apply correlation decay
    Sigma <- rho^exponent
    
    # C. Generate data using that "instruction" matrix (Sigma)
    X_raw <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  }
  
  # 3. Pre-processing (Scaling)
  # Crucial to do this BEFORE generating 'y'.
  if (scale_X) X_raw <- scale(X_raw)
  
  # 4. Final Design Matrix Construction
  X <- cbind(1, X_raw)
  colnames(X) <- names(betas)
  
  # 5. Response Generation (The Power Logit Process)
  eta <- as.numeric(X %*% betas)
  
  # --- POWER LOGIT LINK ---
  # Standard Logit Probability
  prob_logit <- 1 / (1 + exp(-eta))
  
  # Apply Power Parameter 'lambda'
  # if lambda < 1: Curve grows faster (Heavy left tail)
  # if lambda > 1: Curve grows slower (Heavy right tail)
  probs <- prob_logit^lambda
  
  y <- rbinom(n, size = 1, prob = probs)
  
  return(list(
    X = X,
    y = as.vector(y),
    p_true = probs,
    meta = list(seed = seed, rho = rho, n = n, lambda = lambda)
  ))
}

# ------------------------------------------------------------------------------
# Unit Testing / Sanity Checks
# Executes only if the script is called directly (not via source)
# ------------------------------------------------------------------------------
if (sys.nframe() == 0) {
  
  message("Running sanity checks for POWER logit generation...")
  
  # Test setup
  test_n <- 1000
  test_betas <- c("(Intercept)" = 0.0, "x1" = 0.5) # Centered at 0 to test skewness
  
  # 1. Test Standard Logit (lambda=1)
  sim_std <- simulate_power_logit_data(n = test_n, betas = test_betas, lambda = 1, seed = 42)
  mean_std <- mean(sim_std$y)
  
  # 2. Test Skewed Logit (lambda=0.5 -> Should increase probability mass)
  sim_skew <- simulate_power_logit_data(n = test_n, betas = test_betas, lambda = 0.5, seed = 42)
  mean_skew <- mean(sim_skew$y)
  
  # Validation output
  message(sprintf("Mean Y (lambda=1.0): %.3f (Expected ~0.5 for beta=0)", mean_std))
  message(sprintf("Mean Y (lambda=0.5): %.3f (Expected > Mean_std)", mean_skew))
  
  # Logic Verification
  if (mean_skew > mean_std) {
    message("Logic Check Passed: lambda < 1 increased probability mass.")
  } else {
    warning("Logic Check Failed: Expected higher mean for lambda < 1.")
  }
  
  # Structural Check
  stopifnot(
    length(sim_skew$y) == test_n,
    all(sim_skew$y %in% c(0, 1))
  )
  message("Structural checks passed.")
}
