# File: R/mala_driver.R
suppressPackageStartupMessages({
  library(parallel)
  library(RTMB)      
  library(posterior) 
  library(mcmcse)    
  library(coda)
})

#' Función interna que corre una cadena MALA única con RTMB y Adaptación
worker_mala_rtmb <- function(seed, data_list, params) {
  set.seed(seed)
  
  # --- 1. Definir Modelo con RTMB ---
  f_model <- function(beta) {
    X <- data_list$X
    y <- data_list$y
    sigma2 <- data_list$sigma2
    
    log_prior <- sum(dnorm(beta, 0, sqrt(sigma2), log = TRUE))
    eta <- X %*% beta
    log_lik <- sum(dbinom(y, size = 1, prob = 1 / (1 + exp(-eta)), log = TRUE))
    
    return(log_prior + log_lik)
  }
  
  # Crear objeto AD
  p <- data_list$P
  beta_init <- rep(0, p)
  
  # MakeADFun debe ejecutarse DENTRO del worker
  obj <- MakeADFun(f_model, beta_init, silent = TRUE)
  
  get_info <- function(b) {
    val <- obj$fn(b)
    gr  <- obj$gr(b)
    return(list(lp = val, grad = gr))
  }
  
  # --- 2. Configuración MCMC ---
  n_total <- params$n_iter
  n_warmup <- params$n_warmup
  eps <- params$eps_init
  
  draws <- matrix(NA, nrow = n_total, ncol = p)
  accepted <- logical(n_total)
  
  current_beta <- beta_init + rnorm(p, 0, 0.1) 
  info_curr <- get_info(current_beta)
  
  target_accept <- 0.574 
  
  for (i in 1:n_total) {
    z <- rnorm(p)
    grad_curr <- info_curr$grad
    
    # Propuesta MALA
    mu_curr <- current_beta + eps * grad_curr
    prop_beta <- mu_curr + sqrt(2 * eps) * z
    
    info_prop <- get_info(prop_beta)
    
    # Ratio MH
    mu_prop <- prop_beta + eps * info_prop$grad
    log_q_fwd <- sum(dnorm(prop_beta, mean = mu_curr, sd = sqrt(2 * eps), log = TRUE))
    log_q_bwd <- sum(dnorm(current_beta, mean = mu_prop, sd = sqrt(2 * eps), log = TRUE))
    
    log_ratio <- info_prop$lp - info_curr$lp + log_q_bwd - log_q_fwd
    
    acc <- FALSE
    if (log(runif(1)) < log_ratio) {
      current_beta <- prop_beta
      info_curr <- info_prop
      acc <- TRUE
    }
    
    draws[i, ] <- current_beta
    accepted[i] <- acc
    
    # ADAPTACIÓN (Robbins-Monro)
    if (i <= n_warmup) {
      learn_rate <- i^(-0.5) 
      log_eps <- log(eps) + learn_rate * (as.numeric(acc) - target_accept)
      eps <- exp(log_eps)
      eps <- min(max(eps, 1e-6), 1.0)
    }
  }
  
  return(list(
    draws = draws[(n_warmup + 1):n_total, , drop = FALSE],
    final_eps = eps,
    accept_rate = mean(accepted[(n_warmup + 1):n_total]),
    warmup_time = NA 
  ))
}

#' MALA Execution Engine (Parallel)
run_mala_logit_parallel <- function(data_list, params, base_seed) {
  
  # 1. Configurar Cluster
  n_cores <- params$n_chains
  cl <- makeCluster(n_cores)
  
  # 1. Exportamos la función worker.
  # 2. Exportamos 'data_list' y 'params' para que los workers los vean.
  # 3. Usamos envir=environment() para que R busque estas variables dentro de 
  #    ESTA función (run_mala_logit_parallel), no en el entorno global.
  clusterExport(
    cl, 
    c("worker_mala_rtmb", "data_list", "params"), 
    envir = environment()
  )
  
  # Cargar librerías en los workers
  clusterEvalQ(cl, {library(RTMB)})
  
  # 2. Ejecución Paralela
  sys_start <- Sys.time()
  
  seeds <- base_seed + (1:n_cores)
  
  # Ahora el parLapply funcionará porque 'data_list' y 'params' existen en el worker
  res_list <- parLapply(cl, seeds, function(s) {
    worker_mala_rtmb(s, data_list, params)
  })
  
  stopCluster(cl)
  
  sys_end <- Sys.time()
  total_wall_time <- as.numeric(difftime(sys_end, sys_start, units = "secs"))
  
  # 3. Procesamiento de Resultados 
  mcmc_list <- mcmc.list(lapply(res_list, function(r) {
    mat <- r$draws
    if(!is.null(data_list$param_names)) colnames(mat) <- data_list$param_names
    mcmc(mat)
  }))
  
  avg_eps <- mean(sapply(res_list, `[[`, "final_eps"))
  avg_acc <- mean(sapply(res_list, `[[`, "accept_rate"))
  
  # Métricas con 'posterior'
  draws_obj <- as_draws(mcmc_list)
  
  stats_summary <- tryCatch({
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
    
    total_draws <- ndraws(draws_obj)
    summ$iat <- total_draws / summ$ess_bulk
    colnames(summ)[colnames(summ) == "variable"]  <- "parametro"
    
    as.data.frame(summ)
    
  }, error = function(e) {
    message("Error stats: ", e$message)
    return(NULL)
  })
  
  multi_metrics <- tryCatch({
    mat <- as.matrix(mcmc_list)
    list(multivariate_ess = mcmcse::multiESS(mat))
  }, error = function(e) return(list(multivariate_ess = NA)))
  
  return(list(
    resumen = stats_summary,
    muestras = mcmc_list,
    times = list(
      total = total_wall_time,
      warmup = total_wall_time * (params$n_warmup / params$n_iter),
      sample = total_wall_time * ((params$n_iter - params$n_warmup) / params$n_iter)
    ),
    global_metrics = multi_metrics,
    diagnostics = list(
      final_eps = avg_eps,
      accept_rate = avg_acc
    )
  ))
}