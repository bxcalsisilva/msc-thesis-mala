# scripts/run_stan_thesis_pipeline.R
# Pipeline Stan completo — Shill Bidding Application
# Modelos: Logit Estándar | Power Logit | Reversal Power Logit
#
# Outputs (todos archivos planos, legibles sin R):
#   output/tables/
#     stan_shill_summary_{model}.csv          — estimaciones por parámetro
#     stan_convergence_comparison.csv         — Rhat, ESS, tiempo por modelo
#     stan_loo_comparison.csv                 — selección de modelo (LOO-ELPD)
#     stan_power_logit_estimates.csv          — estimaciones Power Logit para tesis
#     stan_lambda_summary.csv                 — resumen posterior de lambda por modelo
#   output/loo/
#     loo_diagnostics_{model}.csv             — Pareto-k + ELPD por modelo
#     loo_influential_obs_{model}.csv         — observaciones influyentes (k > 0.5)
#   output/
#     stan_results_summary.txt                — resumen interpretado completo
#   output/plots/  (RDS para scripts de ggplot)
#     stan_lambda_trace_data.rds
#     stan_estimates_caterpillar_data.rds
#
# Diseño de memoria: el fit Stan se libera DENTRO del loop tras extraer
# todo lo necesario. Solo un modelo ocupa RAM a la vez.
# ==============================================================================
rm(list = ls())

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(rstan)
  library(posterior)
  library(loo)
})

source(here("R", "stan_driver_shill.R"))

rstan_options(auto_write = TRUE)

# ==============================================================================
# 1. DATOS
# ==============================================================================
cat("Cargando datos Shill Bidding...\n")
raw_data <- read.csv(here("data", "Data_shillB.csv"))

X_raw    <- raw_data[, 4:12]
y        <- raw_data$Class
X_scaled <- scale(X_raw)
X_matrix <- cbind(Intercepto = 1, X_scaled)

stan_data <- list(
  N           = nrow(X_matrix),
  P           = ncol(X_matrix),
  X           = X_matrix,
  y           = as.integer(y),
  sigma       = 10,
  param_names = colnames(X_matrix)
)

cat(sprintf("N = %d observaciones | P = %d predictores\n\n", stan_data$N, stan_data$P))

# ==============================================================================
# 2. CONFIGURACIÓN MCMC
# ==============================================================================
PARAMS <- list(
  n_chains = 4,
  n_warmup = 1000,
  n_iter   = 4000   # 1500 post-warmup por cadena = 6000 totales
)

BASE_SEED     <- 12345
MODELS_TO_RUN <- c("standard", "power", "reversal")

stan_files <- list(
  standard = here("models", "logit_shill.stan"),
  power    = here("models", "power_logit_shill.stan"),
  reversal = here("models", "reversal_power_logit_shill.stan")
)

for (d in c("output/tables", "output/plots", "output/loo")) {
  dir.create(here(d), showWarnings = FALSE, recursive = TRUE)
}

# ==============================================================================
# 3. COMPILACIÓN (una sola vez antes del loop)
# ==============================================================================
cat("Compilando modelos Stan (se cachean tras la primera compilación)...\n")
compiled_models <- list()
for (model in MODELS_TO_RUN) {
  cat(sprintf("  -> %s\n", basename(stan_files[[model]])))
  compiled_models[[model]] <- rstan::stan_model(file = stan_files[[model]])
}
cat("Compilación completa.\n\n")

# ==============================================================================
# 4. LOOP PRINCIPAL
# ==============================================================================
RESULTS_LIST    <- list()  # resúmenes de parámetros (sin el fit)
LOO_LIST        <- list()  # objetos loo para loo_compare()
LAMBDA_TRACES   <- list()  # trazas completas de lambda (para RDS de plots)
PARETO_SUMMARY  <- list()  # diagnósticos Pareto-k por modelo (para resumen final)

for (model in MODELS_TO_RUN) {

  cat(strrep("=", 60), "\n", sep = "")
  cat(sprintf("MODELO: %s\n", toupper(model)))
  cat(strrep("=", 60), "\n\n", sep = "")

  # --- 4a. Muestreo ---
  res <- run_stan_shill_parallel(
    data_list  = stan_data,
    params     = PARAMS,
    model_obj  = compiled_models[[model]],
    base_seed  = BASE_SEED,
    model_type = model
  )

  if (is.null(res)) {
    cat(sprintf("Muestreo fallido para %s. Continuando.\n\n", model))
    next
  }

  # --- 4b. CSV de estimaciones de parámetros ---
  if (!is.null(res$resumen)) {
    write.csv(
      res$resumen,
      here("output", "tables", sprintf("stan_shill_summary_%s.csv", model)),
      row.names = FALSE
    )
  }

  ess_lambda <- if (model == "standard") "N/A" else
    sprintf("%.0f", tail(res$resumen$ess_bulk, 1))

  cat(sprintf(
    "Tiempo: %.2f s | Aceptación: %.1f%% | StepSize: %.5f | ESS(lambda): %s\n\n",
    res$times$total,
    res$diagnostics$accept_rate * 100,
    res$diagnostics$final_stepsize,
    ess_lambda
  ))

  # --- 4c. LOO-ELPD (con el fit todavía en memoria) ---
  cat("Calculando LOO-ELPD...\n")

  loo_obj <- tryCatch({
    log_lik <- extract_log_lik(
      res$muestras,
      parameter_name = "log_lik",
      merge_chains   = FALSE
    )
    rel_eff <- relative_eff(exp(log_lik), cores = 1)
    loo_out <- loo(log_lik, r_eff = rel_eff, cores = 1)
    rm(log_lik, rel_eff)
    loo_out
  }, error = function(e) {
    message("LOO falló para '", model, "': ", e$message)
    message("Verifica que 'generated quantities { vector[N] log_lik; }' ",
            "existe en el .stan correspondiente.")
    NULL
  })

  if (!is.null(loo_obj)) {

    print(loo_obj)
    k_vals <- pareto_k_values(loo_obj)

    cat(sprintf("\nDiagnóstico Pareto-k (%s):\n", toupper(model)))
    cat(sprintf("  Bueno      k < 0.5  : %4d  (%.1f%%)\n",
                sum(k_vals < 0.5),                         100 * mean(k_vals < 0.5)))
    cat(sprintf("  Aceptable  0.5-0.7  : %4d  (%.1f%%)\n",
                sum(k_vals >= 0.5 & k_vals < 0.7),         100 * mean(k_vals >= 0.5 & k_vals < 0.7)))
    cat(sprintf("  Malo       0.7-1.0  : %4d  (%.1f%%)\n",
                sum(k_vals >= 0.7 & k_vals < 1.0),         100 * mean(k_vals >= 0.7 & k_vals < 1.0)))
    cat(sprintf("  Muy malo   k >= 1.0 : %4d  (%.1f%%)\n",
                sum(k_vals >= 1.0),                        100 * mean(k_vals >= 1.0)))

    pct_bad <- mean(k_vals > 0.7) * 100
    if (pct_bad > 10)
      cat(sprintf("  ADVERTENCIA: %.1f%% de k malos — LOO puede ser poco fiable.\n", pct_bad))

    # Guardar diagnósticos Pareto-k + ELPD por modelo
    write.csv(
      data.frame(
        modelo      = toupper(model),
        categoria   = c("Bueno (k<0.5)", "Aceptable (0.5-0.7)",
                        "Malo (0.7-1.0)", "Muy malo (k>=1.0)"),
        n_obs       = c(sum(k_vals < 0.5),
                        sum(k_vals >= 0.5 & k_vals < 0.7),
                        sum(k_vals >= 0.7 & k_vals < 1.0),
                        sum(k_vals >= 1.0)),
        pct         = round(c(mean(k_vals < 0.5),
                              mean(k_vals >= 0.5 & k_vals < 0.7),
                              mean(k_vals >= 0.7 & k_vals < 1.0),
                              mean(k_vals >= 1.0)) * 100, 1),
        elpd_loo    = round(loo_obj$estimates["elpd_loo", "Estimate"], 2),
        se_elpd_loo = round(loo_obj$estimates["elpd_loo", "SE"], 2),
        p_loo       = round(loo_obj$estimates["p_loo",    "Estimate"], 2),
        looic       = round(loo_obj$estimates["looic",    "Estimate"], 2)
      ),
      here("output", "loo", sprintf("loo_diagnostics_%s.csv", model)),
      row.names = FALSE
    )

    # Guardar observaciones influyentes (k > 0.5) — útil para diagnóstico
    idx_influential <- which(k_vals > 0.5)
    if (length(idx_influential) > 0) {
      influential_df <- data.frame(
        obs_idx = idx_influential,
        k_value = round(k_vals[idx_influential], 4),
        clase   = stan_data$y[idx_influential]
      )
      influential_df <- influential_df[order(-influential_df$k_value), ]
      write.csv(
        influential_df,
        here("output", "loo", sprintf("loo_influential_obs_%s.csv", model)),
        row.names = FALSE
      )
      cat(sprintf("  %d observaciones influyentes guardadas en loo_influential_obs_%s.csv\n",
                  nrow(influential_df), model))
    }

    # Guardar objeto LOO y acumular para comparación final
    saveRDS(loo_obj, here("output", "loo", sprintf("loo_%s.rds", model)))
    LOO_LIST[[model]] <- loo_obj

    # Guardar diagnóstico para el resumen final
    PARETO_SUMMARY[[model]] <- list(
      elpd     = round(loo_obj$estimates["elpd_loo", "Estimate"], 2),
      se_elpd  = round(loo_obj$estimates["elpd_loo", "SE"], 2),
      p_loo    = round(loo_obj$estimates["p_loo",    "Estimate"], 2),
      looic    = round(loo_obj$estimates["looic",    "Estimate"], 2),
      pct_good = round(mean(k_vals < 0.5) * 100, 1),
      pct_ok   = round(mean(k_vals >= 0.5 & k_vals < 0.7) * 100, 1),
      pct_bad  = round(mean(k_vals >= 0.7) * 100, 1)
    )

    cat(sprintf("\n  -> loo_diagnostics_%s.csv | loo_%s.rds\n\n", model, model))
  }

  # --- 4d. Trazas de lambda (antes de liberar el fit) ---
  if (model != "standard") {
    draws_df <- as_draws_df(as_draws(res$muestras))
    LAMBDA_TRACES[[model]] <- data.frame(
      Iteracion = draws_df$.iteration,
      Cadena    = as.factor(draws_df$.chain),
      Lambda    = as.numeric(draws_df$lambda),
      Modelo    = toupper(model)
    )
  }

  # --- 4e. Liberar el fit ---
  res$muestras <- NULL
  RESULTS_LIST[[model]] <- res
  gc()

  cat(sprintf("Modelo %s completado. RAM liberada.\n\n", toupper(model)))
}

# ==============================================================================
# 5. TABLA DE CONVERGENCIA
# ==============================================================================
convergence_table <- do.call(rbind, lapply(MODELS_TO_RUN, function(m) {
  if (is.null(RESULTS_LIST[[m]])) return(NULL)
  res  <- RESULTS_LIST[[m]]
  summ <- res$resumen
  data.frame(
    Modelo      = toupper(m),
    Tiempo_s    = round(res$times$total, 2),
    Accept_Rate = round(res$diagnostics$accept_rate, 3),
    Step_Size   = round(res$diagnostics$final_stepsize, 5),
    Rhat_Max    = round(max(summ$rhat,     na.rm = TRUE), 4),
    ESS_Min     = round(min(summ$ess_bulk, na.rm = TRUE), 0),
    MultiESS    = round(ifelse(is.null(res$global_metrics$multi_ess), NA,
                               res$global_metrics$multi_ess), 0)
  )
}))

write.csv(convergence_table,
          here("output", "tables", "stan_convergence_comparison.csv"),
          row.names = FALSE)

# ==============================================================================
# 6. TABLA LOO-ELPD
# ==============================================================================
loo_comparison_df <- NULL

if (length(LOO_LIST) >= 2) {
  comparison        <- loo_compare(LOO_LIST)
  loo_comparison_df        <- as.data.frame(comparison)
  loo_comparison_df$model  <- rownames(loo_comparison_df)
  numeric_cols             <- setdiff(names(loo_comparison_df), "model")
  loo_comparison_df[, numeric_cols] <- round(loo_comparison_df[, numeric_cols], 2)
  loo_comparison_df        <- loo_comparison_df[, c("model", numeric_cols)]

  write.csv(loo_comparison_df,
            here("output", "tables", "stan_loo_comparison.csv"),
            row.names = FALSE)
}

# ==============================================================================
# 7. ESTIMACIONES POWER LOGIT
# ==============================================================================
pl_estimates <- NULL

if (!is.null(RESULTS_LIST[["power"]])) {
  pl_estimates <- RESULTS_LIST[["power"]]$resumen %>%
    select(parametro, mean, sd, lower_ci, upper_ci, rhat, ess_bulk)

  write.csv(pl_estimates,
            here("output", "tables", "stan_power_logit_estimates.csv"),
            row.names = FALSE)

  saveRDS(pl_estimates, here("output", "plots", "stan_estimates_caterpillar_data.rds"))
}

# ==============================================================================
# 8. RESUMEN POSTERIOR DE LAMBDA
# ==============================================================================
if (length(LAMBDA_TRACES) > 0) {

  # RDS para scripts de ggplot (traza completa)
  df_lambda_full <- do.call(rbind, LAMBDA_TRACES)
  saveRDS(df_lambda_full, here("output", "plots", "stan_lambda_trace_data.rds"))

  # CSV con estadísticas resumen de lambda por modelo
  lambda_summary <- do.call(rbind, lapply(names(LAMBDA_TRACES), function(m) {
    vals <- LAMBDA_TRACES[[m]]$Lambda
    data.frame(
      modelo  = toupper(m),
      mean    = round(mean(vals), 4),
      sd      = round(sd(vals), 4),
      q2.5    = round(quantile(vals, 0.025), 4),
      q50     = round(quantile(vals, 0.500), 4),
      q97.5   = round(quantile(vals, 0.975), 4),
      p_gt1   = round(mean(vals > 1) * 100, 1),   # % de muestras con lambda > 1
      p_lt1   = round(mean(vals < 1) * 100, 1)    # % de muestras con lambda < 1
    )
  }))

  write.csv(lambda_summary,
            here("output", "tables", "stan_lambda_summary.csv"),
            row.names = FALSE)
}

# ==============================================================================
# 9. RESUMEN INTERPRETADO EN TEXTO PLANO
#    Archivo único que contiene todos los resultados con interpretación.
#    Legible sin R, importable como contexto para Claude.
# ==============================================================================

summary_lines <- c(
  strrep("=", 70),
  "STAN SHILL BIDDING — RESUMEN COMPLETO DE RESULTADOS",
  strrep("=", 70),
  sprintf("Generado: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",

  # --- Datos ---
  strrep("-", 70),
  "DATOS",
  strrep("-", 70),
  sprintf("Dataset : Shill Bidding (UCI ML Repository)"),
  sprintf("N obs   : %d", stan_data$N),
  sprintf("P pred  : %d (incluye intercepto)", stan_data$P),
  sprintf("Balance : %.1f%% positivos (shill bidders)",
          100 * mean(stan_data$y == 1)),
  "",

  # --- MCMC ---
  strrep("-", 70),
  "CONFIGURACIÓN MCMC (Stan HMC/NUTS)",
  strrep("-", 70),
  sprintf("Cadenas             : %d", PARAMS$n_chains),
  sprintf("Warmup por cadena   : %d", PARAMS$n_warmup),
  sprintf("Iteraciones totales : %d", PARAMS$n_iter),
  sprintf("Post-warmup totales : %d",
          (PARAMS$n_iter - PARAMS$n_warmup) * PARAMS$n_chains),
  sprintf("Semilla maestra     : %d", BASE_SEED),
  ""
)

# --- Convergencia ---
if (!is.null(convergence_table)) {
  summary_lines <- c(summary_lines,
    strrep("-", 70),
    "CONVERGENCIA Y EFICIENCIA",
    strrep("-", 70),
    "Criterio: Rhat < 1.01 (estricto) o < 1.05 (aceptable) | ESS > 400",
    ""
  )
  for (i in seq_len(nrow(convergence_table))) {
    r <- convergence_table[i, ]
    rhat_ok <- ifelse(r$Rhat_Max < 1.01, "EXCELENTE",
               ifelse(r$Rhat_Max < 1.05, "ACEPTABLE", "PROBLEMATICO"))
    ess_ok   <- ifelse(r$ESS_Min > 400, "OK", "BAJO")
    summary_lines <- c(summary_lines,
      sprintf("  %s:", r$Modelo),
      sprintf("    Tiempo       = %.2f s", r$Tiempo_s),
      sprintf("    Accept Rate  = %.1f%%", r$Accept_Rate * 100),
      sprintf("    Step Size    = %.5f", r$Step_Size),
      sprintf("    Rhat máximo  = %.4f  [%s]", r$Rhat_Max, rhat_ok),
      sprintf("    ESS mínimo   = %d     [%s]", r$ESS_Min, ess_ok),
      sprintf("    ESS multiv.  = %s",
              ifelse(is.na(r$MultiESS), "N/A", as.character(r$MultiESS))),
      ""
    )
  }
}

# --- LOO-ELPD ---
if (!is.null(loo_comparison_df)) {
  best_model <- loo_comparison_df$model[1]

  summary_lines <- c(summary_lines,
    strrep("-", 70),
    "SELECCIÓN DE MODELO — LOO-ELPD",
    strrep("-", 70),
    "Método: Leave-One-Out Cross-Validation (PSIS-LOO)",
    "Mayor ELPD = mejor capacidad predictiva fuera de muestra.",
    "Diferencia significativa si |elpd_diff| > 2 × se_diff  (|z| > 2)",
    "",
    sprintf("  %-12s  %8s  %8s  %8s  %10s  %8s",
            "Modelo", "ELPD", "SE", "LOOIC", "elpd_diff", "se_diff"),
    sprintf("  %s", strrep("-", 60))
  )

  for (i in seq_len(nrow(loo_comparison_df))) {
    r     <- loo_comparison_df[i, ]
    z_val <- if (i == 1) "—" else
      sprintf("%.2f", abs(r$elpd_diff) / r$se_diff)
    sig   <- if (i == 1) "(referencia)" else
      if (abs(r$elpd_diff) > 2 * r$se_diff) "* significativo" else "  no significativo"

    summary_lines <- c(summary_lines,
      sprintf("  %-12s  %8.2f  %8.2f  %8.2f  %10.2f  %8.2f  z=%s  %s",
              toupper(r$model),
              r$elpd_loo, r$se_elpd_loo, r$looic,
              r$elpd_diff, r$se_diff, z_val, sig)
    )
  }

  summary_lines <- c(summary_lines,
    "",
    sprintf("  MEJOR MODELO: %s", toupper(best_model)),
    ""
  )
}

# --- Pareto-k ---
if (length(PARETO_SUMMARY) > 0) {
  summary_lines <- c(summary_lines,
    strrep("-", 70),
    "DIAGNÓSTICO PARETO-K (Fiabilidad del LOO)",
    strrep("-", 70),
    "k < 0.5: bueno | 0.5-0.7: aceptable | >0.7: estimados poco fiables",
    ""
  )
  for (m in names(PARETO_SUMMARY)) {
    ps <- PARETO_SUMMARY[[m]]
    fiab <- if (ps$pct_bad < 5) "LOO FIABLE" else
      if (ps$pct_bad < 10) "LOO ACEPTABLE" else "LOO POCO FIABLE"
    summary_lines <- c(summary_lines,
      sprintf("  %s  [%s]:", toupper(m), fiab),
      sprintf("    Bueno (k<0.5)    : %.1f%%", ps$pct_good),
      sprintf("    Aceptable (0.5-0.7): %.1f%%", ps$pct_ok),
      sprintf("    Malo/Muy malo (>0.7): %.1f%%", ps$pct_bad),
      sprintf("    ELPD = %.2f (SE = %.2f) | p_loo = %.2f | LOOIC = %.2f",
              ps$elpd, ps$se_elpd, ps$p_loo, ps$looic),
      ""
    )
  }
}

# --- Lambda ---
if (length(LAMBDA_TRACES) > 0 && !is.null(lambda_summary)) {
  summary_lines <- c(summary_lines,
    strrep("-", 70),
    "PARÁMETRO DE ASIMETRÍA LAMBDA",
    strrep("-", 70),
    "Lambda = 1 → logit estándar (simétrico)",
    "Lambda < 1 → la probabilidad crece más rápido (asimetría positiva)",
    "Lambda > 1 → la probabilidad crece más lento (asimetría negativa)",
    ""
  )
  for (i in seq_len(nrow(lambda_summary))) {
    r <- lambda_summary[i, ]
    evidencia <- if (r$q2.5 > 1 | r$q97.5 < 1)
      "ASIMETRIA SIGNIFICATIVA (IC95% no contiene 1)"
    else
      "asimetría no significativa (IC95% contiene 1)"
    summary_lines <- c(summary_lines,
      sprintf("  %s:", r$modelo),
      sprintf("    Media  = %.4f | SD = %.4f", r$mean, r$sd),
      sprintf("    IC 95%% = [%.4f, %.4f]", r$q2.5, r$q97.5),
      sprintf("    P(lambda > 1) = %.1f%% | P(lambda < 1) = %.1f%%",
              r$p_gt1, r$p_lt1),
      sprintf("    → %s", evidencia),
      ""
    )
  }
}

# --- Estimaciones Power Logit ---
if (!is.null(pl_estimates)) {
  summary_lines <- c(summary_lines,
    strrep("-", 70),
    "ESTIMACIONES POSTERIORES — POWER LOGIT",
    strrep("-", 70),
    sprintf("  %-30s  %8s  %8s  %10s  %10s",
            "Parámetro", "Media", "SD", "IC2.5%", "IC97.5%"),
    sprintf("  %s", strrep("-", 72))
  )
  for (i in seq_len(nrow(pl_estimates))) {
    r <- pl_estimates[i, ]
    summary_lines <- c(summary_lines,
      sprintf("  %-30s  %8.4f  %8.4f  %10.4f  %10.4f",
              r$parametro, r$mean, r$sd, r$lower_ci, r$upper_ci)
    )
  }
  summary_lines <- c(summary_lines, "")
}

# --- Archivos generados ---
summary_lines <- c(summary_lines,
  strrep("-", 70),
  "ARCHIVOS GENERADOS",
  strrep("-", 70),
  "output/tables/",
  "  stan_shill_summary_standard.csv      — estimaciones logit estándar",
  "  stan_shill_summary_power.csv         — estimaciones power logit",
  "  stan_shill_summary_reversal.csv      — estimaciones reversal power logit",
  "  stan_convergence_comparison.csv      — Rhat, ESS, tiempo por modelo",
  "  stan_loo_comparison.csv              — comparación LOO-ELPD entre modelos",
  "  stan_power_logit_estimates.csv       — tabla de tesis: estimaciones PL",
  "  stan_lambda_summary.csv              — resumen posterior de lambda",
  "output/loo/",
  "  loo_diagnostics_{model}.csv          — Pareto-k + ELPD por modelo",
  "  loo_influential_obs_{model}.csv      — observaciones con k > 0.5",
  "output/",
  "  stan_results_summary.txt             — este archivo",
  "",
  strrep("=", 70),
  "FIN DEL RESUMEN",
  strrep("=", 70)
)

writeLines(summary_lines, here("output", "stan_results_summary.txt"))

# ==============================================================================
# 10. CONSOLA FINAL
# ==============================================================================
cat("\n")
cat(strrep("=", 60), "\n", sep = "")
cat("PIPELINE COMPLETADO\n")
cat(strrep("=", 60), "\n")
cat("\nArchivos planos generados:\n")
cat("  output/tables/stan_convergence_comparison.csv\n")
cat("  output/tables/stan_loo_comparison.csv\n")
cat("  output/tables/stan_power_logit_estimates.csv\n")
cat("  output/tables/stan_lambda_summary.csv\n")
cat("  output/loo/loo_diagnostics_{model}.csv  (x3)\n")
cat("  output/loo/loo_influential_obs_{model}.csv  (si aplica)\n")
cat("  output/stan_results_summary.txt  <-- resumen interpretado completo\n")
