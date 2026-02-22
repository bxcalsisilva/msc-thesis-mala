# File: R/stan_driver_power.R
suppressPackageStartupMessages({
  library(rstan)
  library(posterior) 
  library(mcmcse)    
})

# Configuración recomendada
rstan_options(auto_write = TRUE)

#' Stan Execution Engine for POWER LOGIT (Parallel)
#' 
#' @param data_list List containing data (X, y, N, P, sigma).
#' @param params List containing MCMC settings (chains, iter, warmup).
#' @param model_obj Compiled Stan model object.
#' @param base_seed Integer. Master seed.
#' 
#' @return List containing summary dataframe, MCMC object, timings, and metrics.

run_stan_power_logit_parallel <- function(data_list, params, model_obj, base_seed) {
  
  # 1. Execution
  sys_start <- Sys.time()
  
  stan_fit <- tryCatch({
    sampling(
      object = model_obj,
      data = data_list,
      chains = params$n_chains,
      iter = params$n_iter,     
      warmup = params$n_warmup,
      cores = params$n_chains,  
      seed = base_seed,
      refresh = 0,              
      show_messages = FALSE
    )
  }, error = function(e) {
    message("Critical Stan Error: ", e$message)
    return(NULL)
  })
  
  sys_end <- Sys.time()
  total_wall_time <- as.numeric(difftime(sys_end, sys_start, units = "secs"))
  
  if (is.null(stan_fit)) return(NULL)
  
  # 2. Extract Timings
  all_times <- get_elapsed_time(stan_fit)
  time_warmup <- max(all_times[, "warmup"])
  time_sample <- max(all_times[, "sample"])
  
  # 3. Output Processing
  
  # A. Convert to 'posterior' draw object
  draws <- as_draws(stan_fit)
  
  # --- CAMBIO CLAVE: Extraemos 'beta' Y 'lambda' ---
  # Esto asegura que el resumen incluya el parámetro de asimetría
  draws_params <- subset_draws(draws, variable = c("beta", "lambda")) 
  
  # --- Calculate Univariate Metrics ---
  stats_summary <- tryCatch({
    
    # Funciones seguras para cuantiles
    q2.5  <- function(x) quantile(x, probs = 0.025, names = FALSE)
    q97.5 <- function(x) quantile(x, probs = 0.975, names = FALSE)
    
    summ <- summarise_draws(
      draws_params,   # Usamos el objeto con betas y lambda
      "mean", "sd", "rhat",
      "ess_bulk", "ess_tail",
      "mcse_mean", "mcse_sd",
      lower_ci = q2.5, 
      upper_ci = q97.5
    )
    
    # IAT
    total_draws <- ndraws(draws_params)
    summ$iat <- total_draws / summ$ess_bulk
    
    # Renombrar columnas
    colnames(summ)[colnames(summ) == "variable"]  <- "parametro"
    
    as.data.frame(summ)
    
  }, error = function(e) {
    message("Error calculating univariate stats: ", e$message)
    return(NULL)
  })
  
  # --- Calculate Multivariate Metrics ---
  # Incluimos lambda en la matriz para evaluar la mezcla conjunta
  multi_metrics <- tryCatch({
    mat <- as.matrix(stan_fit, pars = c("beta", "lambda"))
    list(multivariate_ess = mcmcse::multiESS(mat))
  }, error = function(e) return(list(multivariate_ess = NA)))
  
  return(list(
    resumen = stats_summary,
    muestras = stan_fit, 
    times = list(
      total = total_wall_time, 
      warmup = time_warmup, 
      sample = time_sample
    ),
    global_metrics = multi_metrics
  ))
}