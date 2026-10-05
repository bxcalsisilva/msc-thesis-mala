# File: R/stan_driver_shill.R
# Driver Stan modular para los 3 modelos de enlace:
#   - "standard"  -> Logit estándar         (solo beta)
#   - "power"     -> Power Logit            (beta + lambda)
#   - "reversal"  -> Reversal Power Logit   (beta + lambda)

suppressPackageStartupMessages({
  library(rstan)
  library(posterior)
  library(mcmcse)
})

# Configuración recomendada para Stan en C++
rstan_options(auto_write = TRUE)

# ==============================================================================
# FUNCIÓN GENÉRICA (Nueva) — Soporta los 3 modelos
# ==============================================================================

run_stan_shill_parallel <- function(data_list, params, model_obj, base_seed,
                                    model_type = "power") {
  
  cat(sprintf("\nIniciando muestreo Stan (HMC/NUTS) — Modelo: %s\n", toupper(model_type)))
  sys_start <- Sys.time()
  
  # --------------------------------------------------------------------------
  # 1. Ejecución del muestreo
  # --------------------------------------------------------------------------
  adapt_delta   <- if (!is.null(params$adapt_delta))   params$adapt_delta   else 0.9
  max_treedepth <- if (!is.null(params$max_treedepth)) params$max_treedepth else 12

  cat(sprintf("  adapt_delta=%.2f | max_treedepth=%d\n", adapt_delta, max_treedepth))

  stan_fit <- tryCatch({
    sampling(
      object        = model_obj,
      data          = data_list,
      chains        = params$n_chains,
      iter          = params$n_iter,
      warmup        = params$n_warmup,
      cores         = params$n_chains,
      seed          = base_seed,
      refresh       = max(1, params$n_iter / 10),
      show_messages = TRUE,
      control       = list(adapt_delta   = adapt_delta,
                           max_treedepth = max_treedepth)
    )
  }, error = function(e) {
    message("Error crítico en Stan: ", e$message)
    return(NULL)
  })

  sys_end <- Sys.time()
  total_wall_time <- as.numeric(difftime(sys_end, sys_start, units = "secs"))

  if (is.null(stan_fit)) return(NULL)

  # --------------------------------------------------------------------------
  # 2. Extracción de tiempos internos de Stan
  # --------------------------------------------------------------------------
  all_times   <- get_elapsed_time(stan_fit)
  time_warmup <- max(all_times[, "warmup"])
  time_sample <- max(all_times[, "sample"])
  
  # --------------------------------------------------------------------------
  # 3. Métricas del sampler (análogo a accept_rate y final_eps de MALA)
  # --------------------------------------------------------------------------
  sp <- get_sampler_params(stan_fit, inc_warmup = FALSE)
  
  accept_stat_val <- tryCatch(
    mean(sapply(sp, function(x) mean(x[, "accept_stat__"]))),
    error = function(e) NA
  )
  
  # Tamaño de paso final: promedio del último valor por cadena
  final_stepsize_val <- tryCatch(
    mean(sapply(sp, function(x) tail(x[, "stepsize__"], 1))),
    error = function(e) NA
  )
  
  # --------------------------------------------------------------------------
  # 4. Variables a extraer según el modelo
  # --------------------------------------------------------------------------
  vars_to_extract <- if (model_type == "standard") "beta" else c("beta", "lambda")
  
  draws        <- as_draws(stan_fit)
  draws_params <- subset_draws(draws, variable = vars_to_extract)
  
  # --------------------------------------------------------------------------
  # 5. Resumen estadístico (mismas columnas que mala_driver_shill.R)
  # --------------------------------------------------------------------------
  stats_summary <- tryCatch({
    q2.5  <- function(x) quantile(x, probs = 0.025, names = FALSE)
    q97.5 <- function(x) quantile(x, probs = 0.975, names = FALSE)
    
    summ <- summarise_draws(
      draws_params,
      "mean", "sd",
      lower_ci = q2.5,
      upper_ci = q97.5,
      "rhat",
      "ess_bulk",
      "ess_tail",
      "mcse_mean",
      "mcse_sd"
    )
    
    # Renombrar 'variable' → 'parametro'
    colnames(summ)[colnames(summ) == "variable"] <- "parametro"
    
    # Reemplazar "beta[1]", "beta[2]", ... por los nombres reales
    if (!is.null(data_list$param_names)) {
      nombres_stan <- paste0("beta[", seq_along(data_list$param_names), "]")
      summ$parametro <- as.character(summ$parametro)
      for (i in seq_along(nombres_stan)) {
        summ$parametro[summ$parametro == nombres_stan[i]] <- data_list$param_names[i]
      }
    }
    
    as.data.frame(summ)
  }, error = function(e) {
    message("Error calculando estadísticos: ", e$message)
    return(NULL)
  })
  
  if (is.null(stats_summary)) {
    return(list(resumen = NULL, muestras = stan_fit,
                times = list(total = total_wall_time, warmup = time_warmup, sample = time_sample)))
  }
  
  # --------------------------------------------------------------------------
  # 6. Métricas derivadas (IAT, ESS/seg)
  # --------------------------------------------------------------------------
  total_draws <- ndraws(draws_params)
  
  stats_summary$iat              <- total_draws / stats_summary$ess_bulk
  stats_summary$bulk_ess_per_sec <- stats_summary$ess_bulk / total_wall_time
  stats_summary$tail_ess_per_sec <- stats_summary$ess_tail / total_wall_time
  
  # --------------------------------------------------------------------------
  # 7. Métricas globales (repetidas por fila para exportación CSV plana)
  # --------------------------------------------------------------------------
  stats_summary$time_total      <- total_wall_time
  stats_summary$time_warmup     <- time_warmup
  stats_summary$time_sample     <- time_sample
  stats_summary$accept_rate     <- accept_stat_val
  stats_summary$final_stepsize  <- final_stepsize_val
  
  # ESS Multivariado
  multi_ess_val <- tryCatch({
    draws_mat <- as.matrix(stan_fit, pars = vars_to_extract)
    mcmcse::multiESS(draws_mat)
  }, error = function(e) NA)
  
  stats_summary$multi_ess <- multi_ess_val

  # --------------------------------------------------------------------------
  # 8. Ordenar columnas (mismo orden que mala_driver_shill.R)
  # --------------------------------------------------------------------------
  cols_order <- c(
    "parametro", "mean", "sd", "rhat", "ess_bulk", "ess_tail",
    "mcse_mean", "mcse_sd", "lower_ci", "upper_ci", "iat",
    "time_total", "time_warmup", "time_sample",
    "final_stepsize", "accept_rate", "multi_ess",
    "bulk_ess_per_sec", "tail_ess_per_sec"
  )
  cols_final    <- intersect(cols_order, colnames(stats_summary))
  stats_summary <- stats_summary[, cols_final]
  
  # --------------------------------------------------------------------------
  # 9. Retorno (misma estructura que mala_driver_shill.R)
  # --------------------------------------------------------------------------
  return(list(
    resumen        = stats_summary,
    muestras       = stan_fit,
    times          = list(total = total_wall_time, warmup = time_warmup, sample = time_sample),
    global_metrics = list(multi_ess = multi_ess_val),
    diagnostics    = list(accept_rate = accept_stat_val, final_stepsize = final_stepsize_val),
    model_type     = model_type
  ))
  
}

# ==============================================================================
# FUNCIÓN LEGADO — Mantiene compatibilidad con run_stan_shill.R existente
# ==============================================================================

run_stan_power_shill_parallel <- function(data_list, params, model_obj, base_seed) {
  
  cat("\nIniciando muestreo con Stan (HMC/NUTS)...\n")
  sys_start <- Sys.time()
  
  stan_fit <- tryCatch({
    sampling(
      object        = model_obj,
      data          = data_list,
      chains        = params$n_chains,
      iter          = params$n_iter,
      warmup        = params$n_warmup,
      cores         = params$n_chains,
      seed          = base_seed,
      refresh       = max(1, params$n_iter / 10),
      show_messages = TRUE
    )
  }, error = function(e) {
    message("Error crítico en Stan: ", e$message)
    return(NULL)
  })
  
  sys_end         <- Sys.time()
  total_wall_time <- as.numeric(difftime(sys_end, sys_start, units = "secs"))
  
  if (is.null(stan_fit)) return(NULL)
  
  all_times   <- get_elapsed_time(stan_fit)
  time_warmup <- max(all_times[, "warmup"])
  time_sample <- max(all_times[, "sample"])
  
  draws        <- as_draws(stan_fit)
  draws_params <- subset_draws(draws, variable = c("beta", "lambda"))
  
  stats_summary <- tryCatch({
    q2.5  <- function(x) quantile(x, probs = 0.025, names = FALSE)
    q97.5 <- function(x) quantile(x, probs = 0.975, names = FALSE)
    
    summ <- summarise_draws(
      draws_params,
      "mean", "sd",
      lower_ci = q2.5,
      upper_ci = q97.5,
      "rhat",
      "ess_bulk"
    )
    
    colnames(summ)[colnames(summ) == "variable"] <- "parametro"
    
    if (!is.null(data_list$param_names)) {
      nombres_stan <- paste0("beta[", seq_along(data_list$param_names), "]")
      summ$parametro <- as.character(summ$parametro)
      for (i in seq_along(nombres_stan)) {
        summ$parametro[summ$parametro == nombres_stan[i]] <- data_list$param_names[i]
      }
    }
    as.data.frame(summ)
  }, error = function(e) {
    message("Error calculando estadísticos: ", e$message)
    return(NULL)
  })
  
  return(list(
    resumen  = stats_summary,
    muestras = stan_fit,
    times    = list(total = total_wall_time, warmup = time_warmup, sample = time_sample)
  ))
}
