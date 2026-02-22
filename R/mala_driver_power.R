# File: R/mala_driver_power.R
suppressPackageStartupMessages({
  library(parallel)
  library(RTMB)      
  library(posterior) 
  library(mcmcse)    
  library(coda)
})

#' Función interna que corre una cadena MALA única con RTMB y Adaptación
worker_mala_power_rtmb <- function(seed, data_list, params) {
  set.seed(seed)
  
  # --- 1. Definir Modelo con RTMB ---
  # Theta vector: [beta_1, ..., beta_p, log_lambda]
  f_model <- function(theta) {
    # Desempaquetar parámetros
    p <- data_list$P
    beta <- theta[1:p]
    log_lambda <- theta[p + 1]
    lambda <- exp(log_lambda)
    
    X <- data_list$X
    y <- data_list$y
    sigma2 <- data_list$sigma2
    
    # A. Priors
    # Beta ~ Normal(0, sigma2)
    log_prior_beta <- sum(dnorm(beta, 0, sqrt(sigma2), log = TRUE))
    # Lambda ~ LogNormal(0, 1) => log(lambda) ~ Normal(0, 1)
    log_prior_lambda <- dnorm(log_lambda, 0, 1, log = TRUE)
    
    # B. Likelihood Power Logit
    # p_base = 1 / (1 + exp(-eta))
    eta <- X %*% beta
    
    # RTMB tiene la función plogis(x) = 1/(1+exp(-x))
    # Es numéricamente estable
    p_base <- plogis(eta)
    
    # Transformación Power: p = p_base ^ lambda
    p_power <- p_base ^ lambda
    
    # Bernoulli Log-Likelihood
    # Evitamos log(0) con un epsilon muy pequeño si fuera necesario, 
    # pero dbinom suele manejarlo bien.
    log_lik <- sum(dbinom(y, size = 1, prob = p_power, log = TRUE))
    
    return(log_prior_beta + log_prior_lambda + log_lik)
  }
  
  # Dimensiones
  p <- data_list$P
  n_params <- p + 1 # Betas + log_lambda
  
  # Inicialización
  # Betas en 0, Lambda en 1 (log_lambda = 0)
  theta_init <- rep(0, n_params)
  
  # MakeADFun DENTRO del worker
  obj <- MakeADFun(f_model, theta_init, silent = TRUE)
  
  get_info <- function(th) {
    val <- obj$fn(th)
    gr  <- obj$gr(th)
    return(list(lp = val, grad = gr))
  }
  
  # --- 2. Configuración MCMC ---
  n_total <- params$n_iter
  n_warmup <- params$n_warmup
  eps <- params$eps_init
  
  draws <- matrix(NA, nrow = n_total, ncol = n_params)
  accepted <- logical(n_total)
  
  # Estado inicial con ruido pequeño
  current_theta <- theta_init + rnorm(n_params, 0, 0.1) 
  info_curr <- get_info(current_theta)
  
  target_accept <- 0.574 
  
  for (i in 1:n_total) {
    z <- rnorm(n_params)
    grad_curr <- info_curr$grad
    
    # Propuesta MALA
    mu_curr <- current_theta + eps * grad_curr
    prop_theta <- mu_curr + sqrt(2 * eps) * z
    
    info_prop <- get_info(prop_theta)
    
    # Ratio MH
    mu_prop <- prop_theta + eps * info_prop$grad
    
    log_q_fwd <- sum(dnorm(prop_theta, mean = mu_curr, sd = sqrt(2 * eps), log = TRUE))
    log_q_bwd <- sum(dnorm(current_theta, mean = mu_prop, sd = sqrt(2 * eps), log = TRUE))
    
    log_ratio <- info_prop$lp - info_curr$lp + log_q_bwd - log_q_fwd
    
    acc <- FALSE
    if (log(runif(1)) < log_ratio) {
      current_theta <- prop_theta
      info_curr <- info_prop
      acc <- TRUE
    }
    
    draws[i, ] <- current_theta
    accepted[i] <- acc
    
    # ADAPTACIÓN (Robbins-Monro)
    if (i <= n_warmup) {
      learn_rate <- i^(-0.5) 
      log_eps <- log(eps) + learn_rate * (as.numeric(acc) - target_accept)
      eps <- exp(log_eps)
      eps <- min(max(eps, 1e-6), 1.0)
    }
  }
  
  # Retorno crudo (log_lambda incluido)
  return(list(
    draws = draws[(n_warmup + 1):n_total, , drop = FALSE],
    final_eps = eps,
    accept_rate = mean(accepted[(n_warmup + 1):n_total]),
    warmup_time = NA 
  ))
}

#' MALA Execution Engine for POWER LOGIT (Parallel)
run_mala_power_logit_parallel <- function(data_list, params, base_seed) {
  
  # 1. Configurar Cluster
  n_cores <- params$n_chains
  cl <- makeCluster(n_cores)
  
  # Exportar explícitamente variables y funciones
  # 1. Exportamos la función worker.
  # 2. Exportamos 'data_list' y 'params' para que los workers los vean.
  # 3. Usamos envir=environment() para que R busque estas variables dentro de 
  #    ESTA función (run_mala_logit_parallel), no en el entorno global.
  clusterExport(
    cl, 
    c("worker_mala_power_rtmb", "data_list", "params"), 
    envir = environment()
  )
  
  # Cargar librerías
  clusterEvalQ(cl, {library(RTMB)})
  
  # 2. Ejecución Paralela
  sys_start <- Sys.time()
  
  seeds <- base_seed + (1:n_cores)
  
  # Ahora el parLapply funcionará porque 'data_list' y 'params' existen en el worker
  res_list <- parLapply(cl, seeds, function(s) {
    worker_mala_power_rtmb(s, data_list, params)
  })
  
  stopCluster(cl)
  
  sys_end <- Sys.time()
  total_wall_time <- as.numeric(difftime(sys_end, sys_start, units = "secs"))
  
  # 3. Procesamiento y Transformación (Log_Lambda -> Lambda)
  
  mcmc_list <- mcmc.list(lapply(res_list, function(r) {
    mat <- r$draws
    
    # TRANSFORMACIÓN: Exponenciar la última columna (log_lambda -> lambda)
    k <- ncol(mat)
    mat[, k] <- exp(mat[, k])
    
    # Nombres
    if(!is.null(data_list$param_names)) {
      # Asumimos que param_names trae los nombres de betas. Agregamos lambda.
      colnames(mat) <- c(data_list$param_names, "lambda")
    }
    
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