# ==============================================================================
# R/mala_driver_shill.R
# MALA DRIVER MODULAR: SHILL BIDDING + REPORTING EXTENDIDO
# ==============================================================================

suppressPackageStartupMessages({
  library(parallel)
  library(RTMB)      
  library(coda)
  library(posterior) 
  library(mcmcse)    
})

# ==============================================================================
# 1. WORKERS (Lógica matemática por enlace)
# ==============================================================================

worker_mala_shill_standard <- function(seed, data_list, params) {
  set.seed(seed)
  f_model <- function(theta) {
    beta <- theta
    log_prior_beta <- sum(dnorm(beta, 0, sqrt(data_list$sigma2_beta), log = TRUE))
    # logit link
    log_lik <- sum(dbinom(data_list$y, size = 1, prob = plogis(data_list$X %*% beta), log = TRUE))
    return(log_prior_beta + log_lik)
  }
  run_adaptive_mala_loop(f_model, rep(0, data_list$p_beta), params)
}

worker_mala_shill_power <- function(seed, data_list, params) {
  set.seed(seed)
  f_model <- function(theta) {
    p <- data_list$p_beta
    beta <- theta[1:p]; lambda <- exp(theta[p + 1])
    log_priors <- sum(dnorm(beta, 0, sqrt(data_list$sigma2_beta), log = TRUE)) + 
      dnorm(theta[p+1], 0, sqrt(data_list$sigma2_delta), log = TRUE)
    # power logit link
    p_power <- plogis(data_list$X %*% beta) ^ lambda
    log_lik <- sum(dbinom(data_list$y, size = 1, prob = p_power, log = TRUE))
    return(log_priors + log_lik)
  }
  run_adaptive_mala_loop(f_model, c(rep(0, data_list$p_beta), 0), params)
}

worker_mala_shill_reversal <- function(seed, data_list, params) {
  set.seed(seed)
  f_model <- function(theta) {
    p <- data_list$p_beta
    beta <- theta[1:p]; lambda <- exp(theta[p + 1])
    log_priors <- sum(dnorm(beta, 0, sqrt(data_list$sigma2_beta), log = TRUE)) + 
      dnorm(theta[p+1], 0, sqrt(data_list$sigma2_delta), log = TRUE)
    # reversal power logit link
    p_rev <- 1 - (plogis(-(data_list$X %*% beta)) ^ lambda)
    log_lik <- sum(dbinom(data_list$y, size = 1, prob = p_rev, log = TRUE))
    return(log_priors + log_lik)
  }
  run_adaptive_mala_loop(f_model, c(rep(0, data_list$p_beta), 0), params)
}

# ==============================================================================
# 2. MOTOR MALA (Loop Adaptativo)
# ==============================================================================

run_adaptive_mala_loop <- function(f_model, theta_init, params) {
  # Crear cinta AD en el nodo local
  obj <- MakeADFun(f_model, theta_init, silent = TRUE)
  get_info <- function(th) list(lp=obj$fn(th), grad=obj$gr(th))
  
  n_total <- params$n_iter; n_warmup <- params$n_warmup; eps <- params$eps_init
  current_theta <- theta_init + rnorm(length(theta_init), 0, 0.05)
  info_curr <- get_info(current_theta)
  draws <- matrix(NA, n_total, length(theta_init))
  accepted <- logical(n_total)
  
  for (i in 1:n_total) {
    z <- rnorm(length(theta_init))
    mu_curr <- current_theta + eps * info_curr$grad
    prop_theta <- mu_curr + sqrt(2*eps)*z
    info_prop <- get_info(prop_theta)
    mu_prop <- prop_theta + eps * info_prop$grad
    
    log_ratio <- (info_prop$lp - info_curr$lp) + 
      (sum(dnorm(current_theta, mu_prop, sqrt(2*eps), log=TRUE)) - 
         sum(dnorm(prop_theta, mu_curr, sqrt(2*eps), log=TRUE)))
    
    acc <- FALSE
    if (is.finite(log_ratio) && log(runif(1)) < log_ratio) {
      current_theta <- prop_theta; info_curr <- info_prop; acc <- TRUE
    }
    draws[i,] <- current_theta
    accepted[i] <- acc
    
    if (i <= n_warmup) {
      eps <- exp(log(eps) + (i^-0.5)*(as.numeric(acc) - 0.574))
      eps <- min(max(eps, 1e-7), 2.0)
    }
  }
  return(list(
    draws = draws[(n_warmup+1):n_total,,drop=FALSE], 
    final_eps = eps, 
    accept_rate = mean(accepted[(n_warmup+1):n_total])
  ))
}

# ==============================================================================
# 3. DISPATCHER CON REPORTE EXTENDIDO
# ==============================================================================

run_mala_shill_parallel <- function(data_list, params, model_type = "power", n_cores = 4) {
  
  worker_name <- switch(model_type,
                        "standard" = "worker_mala_shill_standard",
                        "power"    = "worker_mala_shill_power",
                        "reversal" = "worker_mala_shill_reversal"
  )
  
  cl <- makeCluster(n_cores)
  clusterExport(cl, c(worker_name, "run_adaptive_mala_loop", "data_list", "params"), envir=environment())
  clusterEvalQ(cl, library(RTMB))
  
  seeds <- sample.int(1e6, n_cores)
  sys_start <- Sys.time()
  
  # Ejecución en paralelo
  res_list <- parLapply(cl, seeds, function(s) get(worker_name)(s, data_list, params))
  
  stopCluster(cl)
  sys_end <- Sys.time()
  
  # --------------------------------------------------------------------------
  # CÁLCULOS DE TIEMPOS Y MÉTRICAS GLOBALES
  # --------------------------------------------------------------------------
  time_total_val  <- as.numeric(difftime(sys_end, sys_start, units="secs"))
  time_warmup_val <- time_total_val * (params$n_warmup / params$n_iter)
  time_sample_val <- time_total_val * ((params$n_iter - params$n_warmup) / params$n_iter)
  
  final_eps_val   <- mean(sapply(res_list, `[[`, "final_eps"))
  accept_rate_val <- mean(sapply(res_list, `[[`, "accept_rate"))
  
  # --------------------------------------------------------------------------
  # PROCESAMIENTO DE CADENAS
  # --------------------------------------------------------------------------
  mcmc_list <- mcmc.list(lapply(res_list, function(r) {
    mat <- r$draws
    # Transformación lambda si aplica
    if(model_type != "standard") mat[,ncol(mat)] <- exp(mat[,ncol(mat)]) 
    
    # Nombres
    nms <- data_list$param_names
    if(!is.null(nms)) {
      if(model_type == "standard") colnames(mat) <- nms
      else colnames(mat) <- c(nms, "lambda")
    }
    mcmc(mat)
  }))
  
  # --------------------------------------------------------------------------
  # GENERACIÓN DE "RESUMEN" CON TODOS LOS CAMPOS
  # --------------------------------------------------------------------------
  draws_obj <- as_draws(mcmc_list)
  
  # 1. Métricas por parámetro (posterior)
  stats_summary <- tryCatch({
    # Funciones para intervalos de credibilidad
    q2.5  <- function(x) quantile(x, probs = 0.025, names = FALSE)
    q97.5 <- function(x) quantile(x, probs = 0.975, names = FALSE)
    
    summ <- summarise_draws(
      draws_obj, 
      "mean", "sd", "rhat",
      "ess_bulk", "ess_tail", 
      "mcse_mean", "mcse_sd",
      lower_ci = q2.5, 
      upper_ci = q97.5
    )
    as.data.frame(summ)
  }, error = function(e) { message("Error stats: ", e$message); return(NULL) })
  
  # Renombrar columna 'variable' a 'parametro'
  colnames(stats_summary)[colnames(stats_summary) == "variable"] <- "parametro"
  
  # 2. Agregar métricas derivadas por parámetro
  total_draws <- ndraws(draws_obj)
  stats_summary$iat <- total_draws / stats_summary$ess_bulk
  
  # 3. Agregar métricas de eficiencia computacional
  stats_summary$bulk_ess_per_sec <- stats_summary$ess_bulk / time_total_val
  stats_summary$tail_ess_per_sec <- stats_summary$ess_tail / time_total_val
  
  # 4. Agregar métricas globales (repetidas para cada fila)
  stats_summary$time_total   <- time_total_val
  stats_summary$time_warmup  <- time_warmup_val
  stats_summary$time_sample  <- time_sample_val
  stats_summary$final_eps    <- final_eps_val
  stats_summary$accept_rate  <- accept_rate_val
  
  # ESS Multivariado
  multi_ess_val <- tryCatch({
    mcmcse::multiESS(as.matrix(mcmc_list))
  }, error = function(e) NA)
  stats_summary$multi_ess <- multi_ess_val
  
  # --------------------------------------------------------------------------
  # ORDENAR COLUMNAS SEGÚN SOLICITUD
  # --------------------------------------------------------------------------
  cols_order <- c(
    "parametro", "mean", "sd", "rhat", "ess_bulk", "ess_tail", 
    "mcse_mean", "mcse_sd", "lower_ci", "upper_ci", "iat",
    "time_total", "time_warmup", "time_sample", 
    "final_eps", "accept_rate", "multi_ess", 
    "bulk_ess_per_sec", "tail_ess_per_sec"
  )
  
  # Intersección para evitar errores si falta alguna columna auxiliar
  cols_final <- intersect(cols_order, colnames(stats_summary))
  stats_summary <- stats_summary[, cols_final]
  
  return(list(
    resumen         = stats_summary,
    muestras        = mcmc_list,
    times           = list(total = time_total_val),
    global_metrics  = list(multi_ess = multi_ess_val),
    diagnostics     = list(final_eps = final_eps_val, accept_rate = accept_rate_val),
    model_type      = model_type
  ))
}