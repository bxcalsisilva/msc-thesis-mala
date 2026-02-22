# File: R/data_generation.R

#' Universal Logit Data Generator (Deterministic)
#' 
#' @param n Integer. Number of observations.
#' @param betas Numeric vector. True coefficients (including intercept).
#' @param rho Numeric. AR(1) correlation coefficient.
#' @param seed Integer. Seed for exact reproducibility of this function.
#' @param scale_X Boolean. If TRUE, standardizes covariates (Recommended).
#' 
#' @return A list containing X (Matrix), y (Vector), and metadata.

simulate_logit_data <- function(n, betas, rho = 0, seed = NULL, scale_X = TRUE) {
  
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
    # There is no mathematical relationship between column 1 and column 2.
    X_raw <- matrix(rnorm(n * p), nrow = n, ncol = p)
    
  } else {
    Sigma <- matrix(0, nrow = p, ncol = p)
    
    # --- CORRELATED CASE (Auto-Regressive AR1 Structure) ---
    # Objective: Variables that are close (x1, x2) should be very similar,
    # while distant variables (x1, x50) should have little relationship.
    
    # A. Create the distance matrix (Vectorized vs Loops)
    # The 'outer' function takes indices 1:p and creates a cross-matrix by subtracting them.
    # Position [1,2] -> abs(1 - 2) = 1 (Distance 1)
    # Position [1,5] -> abs(1 - 5) = 4 (Distance 4)
    exponent <- abs(outer(1:p, 1:p, "-"))
    
    # B. Apply correlation decay
    # If rho=0.9:
    # Distance 0 (Diagonal): 0.9^0 = 1    (Total correlation with itself)
    # Distance 1 (Neighbors): 0.9^1 = 0.9  (High correlation)
    # Distance 10 (Distant): 0.9^10 = 0.34 (Low correlation)
    Sigma <- rho^exponent
    
    # C. Generate data using that "instruction" matrix (Sigma)
    X_raw <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  }
  
  # 3. Pre-processing (Scaling) - MOVED TO THE END OF X CREATION
  # It is crucial to do this BEFORE generating 'y'.
  # If we scale afterwards, the original 'betas' would no longer predict 'y' correctly.
  if (scale_X) X_raw <- scale(X_raw)
  
  # 4. Final Design Matrix Construction
  X <- cbind(1, X_raw)
  colnames(X) <- names(betas)
  
  # 5. Response Generation (The Logit Process)
  # Calculate TRUE probability using the processed/scaled X
  eta <- as.numeric(X %*% betas)
  probs <- 1 / (1 + exp(-eta))
  y <- rbinom(n, size = 1, prob = probs)
  
  return(list(
    X = X,
    y = as.vector(y),
    p_true = probs,
    meta = list(seed = seed, rho = rho, n = n)
  ))
}

# ------------------------------------------------------------------------------
# Unit Testing / Sanity Checks
# Executes only if the script is called directly (not via source)
# ------------------------------------------------------------------------------
if (sys.nframe() == 0) {
  
  message("Running sanity checks for data generation...")
  
  # Test setup
  test_n <- 100
  test_betas <- c("(Intercept)" = 0.5, "x1" = 1.2, "x2" = -0.8, "x3" = 0.6)
  
  # 1. Test base case (Independent)
  sim_indep <- simulate_logit_data(n = test_n, betas = test_betas, rho = 0, seed = 42)
  
  # Basic structural validations
  stopifnot(
    nrow(sim_indep$X) == test_n,
    ncol(sim_indep$X) == length(test_betas),
    all(sim_indep$X[, 1] == 1),     # Intercept must be 1
    all(sim_indep$y %in% c(0, 1))   # Response must be binary
  )
  
  # 2. Test correlation (AR1)
  # We use large n so sample correlation converges to theoretical correlation
  sim_corr <- simulate_logit_data(n = 500, betas = test_betas, rho = 0.8, seed = 123)
  
  # Check if sample correlation approaches 0.8
  matrix_x_no_int <- sim_corr$X[, -1]
  cor_obs <- cor(matrix_x_no_int)[1, 2]
  
  message(sprintf("Observed correlation (Target ~0.8): %.3f", cor_obs))
  
  if (abs(cor_obs - 0.8) > 0.05) {
    warning("Correlation structure differs more than expected.")
  } else {
    message("Checks passed.")
  }
}