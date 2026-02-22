# File: R/jags_driver_power.R
suppressPackageStartupMessages({
  library(runjags)
  library(coda)
  library(posterior) 
  library(mcmcse)    
})

#' JAGS Execution Engine for POWER LOGIT (Parallel)
#' 
#' @param data_list List containing data (X, y, N, P, sigma).
#' @param params List containing MCMC settings.
#' @param base_seed Integer. Master seed.
#' 
#' @return List containing summary dataframe, MCMC object, detailed timing, and multivariate metrics.


  
run_jags_power_logit_parallel <- function(data_list, params, base_seed) {
  
  # 1. Setup Initial Values & Seeds
  set.seed(base_seed)
  chain_seeds <- sample.int(1e9, params$n_chains)
  
  inits_list <- lapply(1:params$n_chains, function(i) {
    list(
      .RNG.name = "base::Wichmann-Hill",
      .RNG.seed = chain_seeds[i],
      beta = rep(0, data_list$P),
      lambda = 1.0 # Paper Section 3: lambda=1 recupera el logit estándar. Buen punto de partida.
    )
  })
  
  # 2. Execution (Parallel)
  sys_start <- Sys.time()
  
  out_runjags <- tryCatch({
    run.jags(
      model = "models/power_logit.jags", # Updated model file
      monitor = c("beta", "lambda"),
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
  
  # --- Extract Internal Timings ---
  times <- out_runjags$timetaken
  time_warmup <- as.numeric(times["adapt"]) + as.numeric(times["burnin"])
  time_sample <- as.numeric(times["sample"])
  
  # 3. Output Processing & Metrics Calculation
  
  # A. Convert to 'coda' object
  mcmc_coda <- as.mcmc.list(out_runjags)
  
  # B. Convert to 'posterior' draw object
  draws <- as_draws(mcmc_coda)
  
  # --- Calculate Univariate Metrics (Per Parameter) ---
  stats_summary <- tryCatch({
    
    # Safe quantile functions
    q2.5  <- function(x) quantile(x, probs = 0.025, names = FALSE)
    q97.5 <- function(x) quantile(x, probs = 0.975, names = FALSE)
    
    summ <- summarise_draws(
      draws, 
      "mean", "sd",
      "rhat",
      "ess_bulk", "ess_tail",
      "mcse_mean", "mcse_sd",
      lower_ci = q2.5, 
      upper_ci = q97.5
    )
    
    # Calculate IAT
    total_draws <- ndraws(draws)
    summ$iat <- total_draws / summ$ess_bulk
    
    # Rename columns
    colnames(summ)[colnames(summ) == "variable"]  <- "parametro"
    
    as.data.frame(summ)
    
  }, error = function(e) {
    message("Error calculating univariate stats: ", e$message)
    return(NULL)
  })
  
  # --- Calculate Multivariate Metrics (Global) ---
  multi_metrics <- tryCatch({
    mat <- as.matrix(mcmc_coda)
    m_ess <- mcmcse::multiESS(mat)
    list(multivariate_ess = m_ess)
  }, error = function(e) {
    return(list(multivariate_ess = NA))
  })
  
  return(list(
    resumen = stats_summary,
    muestras = mcmc_coda,
    
    # Detailed timings
    times = list(
      total = total_wall_time,
      warmup = time_warmup,
      sample = time_sample
    ),
    
    # Global metrics
    global_metrics = multi_metrics
  ))
}