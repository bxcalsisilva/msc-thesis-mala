# File: R/jags_driver.R
suppressPackageStartupMessages({
  library(runjags)
  library(coda)
  library(posterior) # The Gold Standard for MCMC diagnostics (used by Stan)
  library(mcmcse)    # For Multivariate ESS
})

#' JAGS Execution Engine (Parallel)
#' 
#' @param data_list List containing data (X, y, N, P, sigma).
#' @param params List containing MCMC settings.
#' @param base_seed Integer. Master seed.
#' 
#' @return List containing summary dataframe, MCMC object, detailed timing, and multivariate metrics.

run_jags_logit_parallel <- function(data_list, params, base_seed) {
  
  # 1. Setup Initial Values & Seeds
  chain_seeds <- base_seed + (1:params$n_chains)
  # 
  inits_list <- lapply(1:params$n_chains, function(i) {
    list(
      .RNG.name = "base::Wichmann-Hill",
      .RNG.seed = chain_seeds[i],
      beta = rep(0, data_list$P) 
    )
  })
  
  # 2. Execution (Parallel)
  sys_start <- Sys.time()
  
  out_runjags <- tryCatch({
      run.jags(
        model = "models/logit.jags", 
        monitor = c("beta"),
        data = data_list,
        inits = inits_list,
        n.chains = params$n_chains,
        adapt = params$n_adapt,
        burnin = params$n_burnin,
        sample = params$n_iter,
        method = "parallel",
        summarise = FALSE,
        silent.jags = TRUE
      )
  }, error = function(e) {
    message("Critical JAGS Error: ", e$message)
    return(NULL)
  })
  
  sys_end <- Sys.time()
  total_wall_time <- as.numeric(difftime(sys_end, sys_start, units = "secs"))
  
  if (is.null(out_runjags)) return(NULL)
  
  # 3. Output Processing & Metrics Calculation
  
  # A. Convert to 'coda' object (for legacy support and mcmcse)
  mcmc_coda <- as.mcmc.list(out_runjags)
  
  # B. Convert to 'posterior' draw object (for Robust/Stan-like metrics)
  # This ensures formulas for Rhat and ESS are identical to Stan's output
  draws <- as_draws(mcmc_coda)
  
  # --- Calculate Univariate Metrics (Per Parameter) ---
  # summarise_draws automatically calculates: 
  # mean, median, sd, mad, q5, q95, rhat, ess_bulk, ess_tail
  stats_summary <- tryCatch({
    
    summ <- summarise_draws(
      draws, 
      "mean", "sd",
      "rhat",
      "ess_bulk", "ess_tail",
      "mcse_mean", "mcse_sd",
      # Specific quantiles for 95% Coverage (2.5% and 97.5%)
      ~quantile(.x, probs = c(0.025, 0.975)) 
    )
    
    # Calculate IAT (Integrated Autocorrelation Time)
    # IAT = N_draws / ESS
    total_draws <- ndraws(draws)
    summ$iat <- total_draws / summ$ess_bulk
    
    # Rename columns to match your existing CSV structure / Preference
    colnames(summ)[colnames(summ) == "variable"]  <- "parametro"
    colnames(summ)[colnames(summ) == "2.5%"]      <- "lower_ci"
    colnames(summ)[colnames(summ) == "97.5%"]     <- "upper_ci"
    
    as.data.frame(summ)
    
  }, error = function(e) {
    message("Error calculating univariate stats: ", e$message)
    return(NULL)
  })
  
  # --- Calculate Multivariate Metrics (Global) ---
  # Multivariate ESS tells you if the joint distribution is explored well.
  # We extract the matrix of samples (Rows=Iter, Cols=Params)
  multi_metrics <- tryCatch({
    mat <- as.matrix(mcmc_coda)
    # Calculate MultiESS using 'mcmcse' package
    m_ess <- mcmcse::multiESS(mat)
    list(multivariate_ess = m_ess)
  }, error = function(e) {
    return(list(multivariate_ess = NA))
  })
  
  return(list(
    resumen = stats_summary,
    muestras = mcmc_coda,
    timetotal = total_wall_time,
    # Global metrics (scalars)
    global_metrics = multi_metrics
  ))
}