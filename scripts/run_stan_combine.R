# scripts/run_stan_combine.R
# Combina los resultados de los tres modelos Stan y genera todos los outputs.
#
# PRE-REQUISITO: haber corrido run_stan_single_model.R para cada modelo:
#   output/loo/loo_standard.rds
#   output/loo/loo_power.rds
#   output/loo/loo_reversal.rds
#   output/tables/stan_shill_summary_{standard,power,reversal}.csv
#   output/tables/stan_lambda_trace_{power,reversal}.csv
#
# Outputs generados:
#   output/tables/stan_loo_comparison.csv
#   output/tables/stan_lambda_summary.csv
#   output/tables/stan_convergence_comparison.csv
#   output/tables/stan_power_logit_estimates.csv
#   output/stan_results_summary.txt
# ==============================================================================
rm(list = ls())

suppressPackageStartupMessages({
  library(here)
  library(loo)
  library(dplyr)
  library(bridgesampling)
})

MODELS <- c("standard", "power", "reversal")

# ==============================================================================
# 1. VERIFICAR QUE EXISTEN LOS ARCHIVOS
# ==============================================================================
missing <- c()
for (m in MODELS) {
  rds <- here("output", "loo", sprintf("loo_%s.rds", m))
  csv <- here("output", "tables", sprintf("stan_shill_summary_%s.csv", m))
  if (!file.exists(rds)) missing <- c(missing, rds)
  if (!file.exists(csv)) missing <- c(missing, csv)
}

if (length(missing) > 0) {
  cat("Faltan los siguientes archivos:\n")
  cat(paste(" -", missing, collapse = "\n"), "\n")
  stop("Corre run_stan_single_model.R para cada modelo antes de combinar.")
}

cat("Todos los archivos encontrados. Combinando resultados...\n\n")

# ==============================================================================
# 2. CARGAR LOO OBJECTS Y COMPARAR
# ==============================================================================
loo_list <- lapply(setNames(MODELS, MODELS), function(m) {
  readRDS(here("output", "loo", sprintf("loo_%s.rds", m)))
})

comparison        <- loo_compare(loo_list)
loo_comparison_df <- as.data.frame(comparison)
loo_comparison_df$model <- rownames(loo_comparison_df)
num_cols <- setdiff(names(loo_comparison_df), "model")
loo_comparison_df[, num_cols] <- round(loo_comparison_df[, num_cols], 2)
loo_comparison_df <- loo_comparison_df[, c("model", num_cols)]

write.csv(loo_comparison_df,
          here("output", "tables", "stan_loo_comparison.csv"),
          row.names = FALSE)
cat("-> stan_loo_comparison.csv\n")

# ==============================================================================
# 2b. BRIDGE SAMPLING — BAYES FACTORS
# ==============================================================================
bridge_files     <- sapply(MODELS, function(m)
  here("output", "loo", sprintf("bridge_%s.rds", m)))
bridge_available <- all(file.exists(bridge_files))

logml    <- NULL
bf_table <- NULL

if (bridge_available) {
  bridge_list <- lapply(setNames(MODELS, MODELS), function(m)
    readRDS(here("output", "loo", sprintf("bridge_%s.rds", m))))

  logml        <- sapply(bridge_list, function(b) b$logml)
  names(logml) <- MODELS

  interp_bf <- function(bf) {
    bfabs <- max(bf, 1 / bf)
    if      (bfabs > 100) "Decisiva"
    else if (bfabs > 30)  "Muy fuerte"
    else if (bfabs > 10)  "Fuerte"
    else if (bfabs > 3)   "Moderada"
    else                  "Debil"
  }

  comparaciones <- list(
    list(m1 = "power",    m2 = "standard", label = "Power vs Standard"),
    list(m1 = "reversal", m2 = "standard", label = "Reversal vs Standard"),
    list(m1 = "power",    m2 = "reversal", label = "Power vs Reversal")
  )

  bf_rows <- lapply(comparaciones, function(x) {
    log_bf <- logml[x$m1] - logml[x$m2]
    bf_val <- exp(log_bf)
    data.frame(
      comparacion = x$label,
      log_BF      = round(log_bf, 3),
      BF          = round(bf_val, 2),
      favorece    = ifelse(bf_val >= 1, toupper(x$m1), toupper(x$m2)),
      evidencia   = interp_bf(bf_val),
      stringsAsFactors = FALSE
    )
  })
  bf_table <- do.call(rbind, bf_rows)

  write.csv(bf_table,
            here("output", "tables", "stan_bayes_factors.csv"),
            row.names = FALSE)
  cat("-> stan_bayes_factors.csv\n")
  cat("\nBAYES FACTORS:\n")
  print(bf_table)
  cat("\n")

} else {
  cat("Bridge sampling opcional: corre run_stan_single_model.R para generarlo.\n\n")
}

# ==============================================================================
# 3. CARGAR RESUMENES DE PARAMETROS
# ==============================================================================
summaries <- lapply(setNames(MODELS, MODELS), function(m) {
  read.csv(here("output", "tables", sprintf("stan_shill_summary_%s.csv", m)))
})

# Tabla de convergencia
convergence_table <- do.call(rbind, lapply(MODELS, function(m) {
  s <- summaries[[m]]
  data.frame(
    Modelo      = toupper(m),
    Rhat_Max    = round(max(s$rhat,     na.rm = TRUE), 4),
    ESS_Min     = round(min(s$ess_bulk, na.rm = TRUE), 0),
    Accept_Rate = round(unique(s$accept_rate)[1], 3),
    Step_Size   = round(unique(s$final_stepsize)[1], 5),
    Tiempo_s    = round(unique(s$time_total)[1], 2)
  )
}))

write.csv(convergence_table,
          here("output", "tables", "stan_convergence_comparison.csv"),
          row.names = FALSE)
cat("-> stan_convergence_comparison.csv\n")

# Estimaciones Power Logit para tabla de tesis
pl_estimates <- summaries[["power"]] %>%
  select(parametro, mean, sd, lower_ci, upper_ci, rhat, ess_bulk)

write.csv(pl_estimates,
          here("output", "tables", "stan_power_logit_estimates.csv"),
          row.names = FALSE)
cat("-> stan_power_logit_estimates.csv\n")

# ==============================================================================
# 4. RESUMEN DE LAMBDA
# ==============================================================================
lambda_summary <- do.call(rbind, lapply(c("power", "reversal"), function(m) {
  trace_file <- here("output", "tables", sprintf("stan_lambda_trace_%s.csv", m))
  if (!file.exists(trace_file)) {
    cat(sprintf("Advertencia: no se encontro %s\n", basename(trace_file)))
    return(NULL)
  }
  vals <- read.csv(trace_file)$lambda
  data.frame(
    modelo = toupper(m),
    mean   = round(mean(vals), 4),
    sd     = round(sd(vals), 4),
    q2.5   = round(quantile(vals, 0.025), 4),
    q50    = round(quantile(vals, 0.500), 4),
    q97.5  = round(quantile(vals, 0.975), 4),
    p_gt1  = round(mean(vals > 1) * 100, 1),
    p_lt1  = round(mean(vals < 1) * 100, 1)
  )
}))

if (!is.null(lambda_summary)) {
  write.csv(lambda_summary,
            here("output", "tables", "stan_lambda_summary.csv"),
            row.names = FALSE)
  cat("-> stan_lambda_summary.csv\n")
}

# ==============================================================================
# 5. DIAGNOSTICOS PARETO-K POR MODELO
# ==============================================================================
pareto_summary <- lapply(setNames(MODELS, MODELS), function(m) {
  diag_file <- here("output", "loo", sprintf("loo_diagnostics_%s.csv", m))
  if (!file.exists(diag_file)) return(NULL)
  d <- read.csv(diag_file)
  list(
    elpd     = d$elpd_loo[1],
    se_elpd  = d$se_elpd_loo[1],
    p_loo    = d$p_loo[1],
    looic    = d$looic[1],
    pct_good = d$pct[d$categoria == "Bueno (k<0.5)"],
    pct_ok   = d$pct[d$categoria == "Aceptable (0.5-0.7)"],
    pct_bad  = 100 - d$pct[d$categoria == "Bueno (k<0.5)"] -
                     d$pct[d$categoria == "Aceptable (0.5-0.7)"]
  )
})

# ==============================================================================
# 6. RESUMEN EN TEXTO PLANO
# ==============================================================================
best_model <- loo_comparison_df$model[1]

lines <- c(
  strrep("=", 70),
  "STAN SHILL BIDDING — RESUMEN COMPLETO DE RESULTADOS",
  strrep("=", 70),
  sprintf("Generado: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",

  strrep("-", 70),
  "DATOS",
  strrep("-", 70),
  "Dataset : Shill Bidding (UCI ML Repository)",
  "N obs   : 6321 | P pred : 10 (incluye intercepto)",
  "Balance : 10.7% positivos (shill bidders)",
  "",

  strrep("-", 70),
  "CONFIGURACION MCMC (Stan HMC/NUTS)",
  strrep("-", 70),
  "Modelos corridos en sesiones R independientes (evita acumulacion de RAM).",
  "  Standard : 2 cadenas | warmup 1000 | iter 5000 |  8000 draws totales",
  "  Power    : 2 cadenas | warmup 1000 | iter 5000 |  8000 draws totales",
  "  Reversal : 2 cadenas | warmup 1000 | iter 5000 |  8000 draws totales",
  "(Power y Reversal usan 2 cadenas para mantener log_lik < 0.8 GB por modelo)",
  ""
)

# Convergencia
lines <- c(lines,
  strrep("-", 70),
  "CONVERGENCIA Y EFICIENCIA",
  strrep("-", 70),
  "Criterio: Rhat < 1.01 (excelente) | ESS > 400",
  ""
)
for (i in seq_len(nrow(convergence_table))) {
  r <- convergence_table[i, ]
  rhat_label <- ifelse(r$Rhat_Max < 1.01, "EXCELENTE",
                ifelse(r$Rhat_Max < 1.05, "ACEPTABLE", "PROBLEMATICO"))
  lines <- c(lines,
    sprintf("  %s:", r$Modelo),
    sprintf("    Tiempo      = %.2f s", r$Tiempo_s),
    sprintf("    Accept Rate = %.1f%%", r$Accept_Rate * 100),
    sprintf("    Step Size   = %.5f", r$Step_Size),
    sprintf("    Rhat max    = %.4f  [%s]", r$Rhat_Max, rhat_label),
    sprintf("    ESS min     = %.0f", r$ESS_Min),
    ""
  )
}

# LOO
lines <- c(lines,
  strrep("-", 70),
  "SELECCION DE MODELO — LOO-ELPD (PSIS-LOO)",
  strrep("-", 70),
  "Mayor ELPD = mejor predictivo. Diferencia significativa si |elpd_diff| > 2 x se_diff",
  "",
  sprintf("  %-12s  %8s  %6s  %8s  %10s  %8s",
          "Modelo", "ELPD", "SE", "LOOIC", "elpd_diff", "se_diff"),
  sprintf("  %s", strrep("-", 58))
)
for (i in seq_len(nrow(loo_comparison_df))) {
  r     <- loo_comparison_df[i, ]
  z_val <- if (i == 1) "—" else sprintf("%.2f", abs(r$elpd_diff) / r$se_diff)
  sig   <- if (i == 1) "(referencia)" else
    if (abs(r$elpd_diff) > 2 * r$se_diff) "* sig." else "  n.s."
  lines <- c(lines,
    sprintf("  %-12s  %8.2f  %6.2f  %8.2f  %10.2f  %8.2f  z=%s  %s",
            toupper(r$model), r$elpd_loo, r$se_elpd_loo, r$looic,
            r$elpd_diff, r$se_diff, z_val, sig)
  )
}
lines <- c(lines, "", sprintf("  MEJOR MODELO (LOO): %s", toupper(best_model)), "")

# Bridge Sampling
if (!is.null(logml) && !is.null(bf_table)) {
  lines <- c(lines,
    strrep("-", 70),
    "BAYES FACTORS (Bridge Sampling)",
    strrep("-", 70),
    "BF > 1: evidencia a favor del primer modelo.",
    "Escala Jeffreys: >3 moderada | >10 fuerte | >30 muy fuerte | >100 decisiva",
    "",
    "  Log-verosimilitud marginal:",
    sprintf("    Standard : %.4f", logml["standard"]),
    sprintf("    Power    : %.4f", logml["power"]),
    sprintf("    Reversal : %.4f", logml["reversal"]),
    ""
  )
  for (i in seq_len(nrow(bf_table))) {
    r <- bf_table[i, ]
    lines <- c(lines,
      sprintf("  %-24s  log_BF=%8.3f  BF=%9.2f  [%-12s]  -> %s",
              r$comparacion, r$log_BF, r$BF, r$evidencia, r$favorece)
    )
  }
  bf_power_rev <- bf_table[bf_table$comparacion == "Power vs Reversal", ]
  if (nrow(bf_power_rev) > 0) {
    lines <- c(lines, "",
      sprintf("  MODELO RECOMENDADO (BF): %s", bf_power_rev$favorece))
  }
  lines <- c(lines, "")
}

# Pareto-k
lines <- c(lines,
  strrep("-", 70),
  "DIAGNOSTICO PARETO-K",
  strrep("-", 70),
  "k < 0.5: bueno | 0.5-0.7: aceptable | > 0.7: LOO poco fiable",
  ""
)
for (m in MODELS) {
  if (is.null(pareto_summary[[m]])) next
  ps    <- pareto_summary[[m]]
  fiab  <- if (ps$pct_bad < 5)  "LOO FIABLE" else
            if (ps$pct_bad < 10) "LOO ACEPTABLE" else "LOO POCO FIABLE"
  lines <- c(lines,
    sprintf("  %s  [%s]:", toupper(m), fiab),
    sprintf("    Bueno (k<0.5)    : %.1f%%", ps$pct_good),
    sprintf("    Aceptable (0.5-0.7): %.1f%%", ps$pct_ok),
    sprintf("    Malo/Muy malo (>0.7): %.1f%%", ps$pct_bad),
    sprintf("    ELPD = %.2f (SE = %.2f) | p_loo = %.2f | LOOIC = %.2f",
            ps$elpd, ps$se_elpd, ps$p_loo, ps$looic),
    ""
  )
}

# Lambda
if (!is.null(lambda_summary) && nrow(lambda_summary) > 0) {
  lines <- c(lines,
    strrep("-", 70),
    "PARAMETRO DE ASIMETRIA LAMBDA",
    strrep("-", 70),
    "Lambda != 1 indica asimetria significativa respecto al logit estandar.",
    ""
  )
  for (i in seq_len(nrow(lambda_summary))) {
    r   <- lambda_summary[i, ]
    ev  <- if (r$q2.5 > 1 | r$q97.5 < 1)
      "ASIMETRIA SIGNIFICATIVA (IC95% no contiene 1)"
    else "asimetria no significativa (IC95% contiene 1)"
    lines <- c(lines,
      sprintf("  %s:", r$modelo),
      sprintf("    Media = %.4f | SD = %.4f | IC95%% = [%.4f, %.4f]",
              r$mean, r$sd, r$q2.5, r$q97.5),
      sprintf("    P(lambda>1) = %.1f%% | P(lambda<1) = %.1f%%",
              r$p_gt1, r$p_lt1),
      sprintf("    -> %s", ev),
      ""
    )
  }
}

# Estimaciones Power Logit
lines <- c(lines,
  strrep("-", 70),
  "ESTIMACIONES POSTERIORES — POWER LOGIT",
  strrep("-", 70),
  sprintf("  %-28s  %8s  %8s  %10s  %10s",
          "Parametro", "Media", "SD", "IC2.5%", "IC97.5%"),
  sprintf("  %s", strrep("-", 70))
)
for (i in seq_len(nrow(pl_estimates))) {
  r <- pl_estimates[i, ]
  lines <- c(lines,
    sprintf("  %-28s  %8.4f  %8.4f  %10.4f  %10.4f",
            r$parametro, r$mean, r$sd, r$lower_ci, r$upper_ci)
  )
}
lines <- c(lines, "")

lines <- c(lines,
  strrep("=", 70),
  "FIN DEL RESUMEN",
  strrep("=", 70)
)

writeLines(lines, here("output", "stan_results_summary.txt"))
cat("-> stan_results_summary.txt\n")

cat("\n")
cat(strrep("=", 60), "\n", sep = "")
cat("COMBINACION COMPLETA\n")
cat(strrep("=", 60), "\n")
cat(sprintf("Mejor modelo (LOO): %s\n", toupper(best_model)))
print(comparison)
