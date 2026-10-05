# scripts/run_mala_combine.R
# Combina los resultados de los tres modelos MALA y genera el resumen.
# Espejo de run_stan_combine.R, pero con las metricas de MALA.
#
# Bayes Factor: NO se recomputa con MALA (la verosimilitud marginal es propiedad
# del modelo, no del sampler). Si existen los bridge_*.rds de Stan, se incluyen
# como REFERENCIA claramente etiquetada (independiente del sampler).
#
# PRE-REQUISITO: run_mala_single_model.R para cada modelo:
#   output/loo/loo_mala_{standard,power,reversal}.rds
#   output/tables/mala_shill_summary_{standard,power,reversal}.csv
#   output/tables/mala_lambda_trace_{power,reversal}.csv
#
# Outputs:
#   output/tables/mala_loo_comparison.csv
#   output/tables/mala_lambda_summary.csv
#   output/tables/mala_convergence_comparison.csv
#   output/tables/mala_power_logit_estimates.csv
#   output/mala_results_summary.txt
# ==============================================================================
rm(list = ls())

suppressPackageStartupMessages({
  library(here)
  library(loo)
  library(dplyr)
})

MODELS <- c("standard", "power", "reversal")

# ==============================================================================
# 1. VERIFICAR ARCHIVOS
# ==============================================================================
missing <- c()
for (m in MODELS) {
  rds <- here("output", "loo", sprintf("loo_mala_%s.rds", m))
  csv <- here("output", "tables", sprintf("mala_shill_summary_%s.csv", m))
  if (!file.exists(rds)) missing <- c(missing, rds)
  if (!file.exists(csv)) missing <- c(missing, csv)
}
if (length(missing) > 0) {
  cat("Faltan los siguientes archivos:\n")
  cat(paste(" -", missing, collapse = "\n"), "\n")
  stop("Corre run_mala_single_model.R para cada modelo antes de combinar.")
}
cat("Todos los archivos encontrados. Combinando resultados MALA...\n\n")

# ==============================================================================
# 2. LOO + COMPARACION
# ==============================================================================
loo_list <- lapply(setNames(MODELS, MODELS), function(m)
  readRDS(here("output", "loo", sprintf("loo_mala_%s.rds", m))))

comparison        <- loo_compare(loo_list)
loo_comparison_df <- as.data.frame(comparison)
loo_comparison_df$model <- rownames(loo_comparison_df)
num_cols <- setdiff(names(loo_comparison_df), "model")
loo_comparison_df[, num_cols] <- round(loo_comparison_df[, num_cols], 2)
loo_comparison_df <- loo_comparison_df[, c("model", num_cols)]

write.csv(loo_comparison_df, here("output", "tables", "mala_loo_comparison.csv"),
          row.names = FALSE)
cat("-> mala_loo_comparison.csv\n")

# ==============================================================================
# 2b. BAYES FACTOR — REFERENCIA desde Stan (independiente del sampler)
# ==============================================================================
bridge_files     <- sapply(MODELS, function(m) here("output","loo",sprintf("bridge_%s.rds", m)))
bridge_available <- all(file.exists(bridge_files))
logml <- NULL; bf_table <- NULL

if (bridge_available) {
  logml        <- sapply(setNames(MODELS, MODELS),
                         function(m) readRDS(here("output","loo",sprintf("bridge_%s.rds", m)))$logml)
  interp_bf <- function(bf) {
    a <- max(bf, 1/bf)
    if (a>100) "Decisiva" else if (a>30) "Muy fuerte" else if (a>10) "Fuerte" else if (a>3) "Moderada" else "Debil"
  }
  comps <- list(c("power","standard"), c("reversal","standard"), c("power","reversal"))
  bf_table <- do.call(rbind, lapply(comps, function(x) {
    lbf <- logml[x[1]] - logml[x[2]]; bfv <- exp(lbf)
    data.frame(comparacion = sprintf("%s vs %s", tools::toTitleCase(x[1]), tools::toTitleCase(x[2])),
               log_BF = round(lbf,3), BF = round(bfv,2),
               favorece = ifelse(bfv>=1, toupper(x[1]), toupper(x[2])),
               evidencia = interp_bf(bfv), stringsAsFactors = FALSE)
  }))
  cat("-> BF (referencia de Stan) cargado\n")
} else {
  cat("(BF de referencia no disponible: faltan bridge_*.rds de Stan)\n")
}

# ==============================================================================
# 3. RESUMENES DE PARAMETROS + CONVERGENCIA
# ==============================================================================
summaries <- lapply(setNames(MODELS, MODELS), function(m)
  read.csv(here("output", "tables", sprintf("mala_shill_summary_%s.csv", m))))

convergence_table <- do.call(rbind, lapply(MODELS, function(m) {
  s <- summaries[[m]]
  data.frame(
    Modelo      = toupper(m),
    Rhat_Max    = round(max(s$rhat,     na.rm = TRUE), 4),
    ESS_Min     = round(min(s$ess_bulk, na.rm = TRUE), 0),
    Accept_Rate = round(unique(s$accept_rate)[1], 3),
    Eps_Final   = round(unique(s$final_eps)[1], 6),
    Tiempo_s    = round(unique(s$time_total)[1], 2)
  )
}))
write.csv(convergence_table, here("output","tables","mala_convergence_comparison.csv"),
          row.names = FALSE)
cat("-> mala_convergence_comparison.csv\n")

pl_estimates <- summaries[["power"]] %>%
  select(parametro, mean, sd, lower_ci, upper_ci, rhat, ess_bulk)
write.csv(pl_estimates, here("output","tables","mala_power_logit_estimates.csv"),
          row.names = FALSE)
cat("-> mala_power_logit_estimates.csv\n")

# ==============================================================================
# 4. RESUMEN DE LAMBDA
# ==============================================================================
lambda_summary <- do.call(rbind, lapply(c("power", "reversal"), function(m) {
  tf <- here("output","tables",sprintf("mala_lambda_trace_%s.csv", m))
  if (!file.exists(tf)) { cat(sprintf("Advertencia: falta %s\n", basename(tf))); return(NULL) }
  v <- read.csv(tf)$lambda
  data.frame(modelo = toupper(m), mean = round(mean(v),4), sd = round(sd(v),4),
             q2.5 = round(quantile(v,.025),4), q50 = round(quantile(v,.5),4),
             q97.5 = round(quantile(v,.975),4),
             p_gt1 = round(mean(v>1)*100,1), p_lt1 = round(mean(v<1)*100,1))
}))
if (!is.null(lambda_summary)) {
  write.csv(lambda_summary, here("output","tables","mala_lambda_summary.csv"), row.names = FALSE)
  cat("-> mala_lambda_summary.csv\n")
}

# ==============================================================================
# 5. PARETO-K POR MODELO
# ==============================================================================
pareto_summary <- lapply(setNames(MODELS, MODELS), function(m) {
  df <- here("output","loo",sprintf("loo_diagnostics_mala_%s.csv", m))
  if (!file.exists(df)) return(NULL)
  d <- read.csv(df)
  list(elpd = d$elpd_loo[1], se_elpd = d$se_elpd_loo[1], p_loo = d$p_loo[1], looic = d$looic[1],
       pct_good = d$pct[d$categoria=="Bueno (k<0.5)"],
       pct_ok   = d$pct[d$categoria=="Aceptable (0.5-0.7)"],
       pct_bad  = 100 - d$pct[d$categoria=="Bueno (k<0.5)"] - d$pct[d$categoria=="Aceptable (0.5-0.7)"])
})

# ==============================================================================
# 6. RESUMEN EN TEXTO
# ==============================================================================
best_model <- loo_comparison_df$model[1]

lines <- c(
  strrep("=", 70),
  "MALA SHILL BIDDING — RESUMEN COMPLETO DE RESULTADOS",
  strrep("=", 70),
  sprintf("Generado: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "Sampler: MALA (RTMB) | Prior log_lambda: Uniforme(-2,2) [rechazo de frontera]",
  "Covariables: 4 (seleccion Fase 1) + intercepto",
  "",
  strrep("-", 70), "DATOS", strrep("-", 70),
  "Dataset : Shill Bidding (UCI ML Repository)",
  "N obs   : 6321 | P pred : 5 (incluye intercepto)",
  "Balance : 10.7% positivos (shill bidders)",
  "",
  strrep("-", 70), "CONVERGENCIA Y EFICIENCIA", strrep("-", 70),
  "Criterio: Rhat < 1.01 (excelente) | ESS > 400", ""
)
for (i in seq_len(nrow(convergence_table))) {
  r <- convergence_table[i, ]
  lab <- ifelse(r$Rhat_Max < 1.01, "EXCELENTE", ifelse(r$Rhat_Max < 1.05, "ACEPTABLE", "PROBLEMATICO"))
  lines <- c(lines,
    sprintf("  %s:", r$Modelo),
    sprintf("    Tiempo      = %.2f s", r$Tiempo_s),
    sprintf("    Accept Rate = %.1f%%", r$Accept_Rate * 100),
    sprintf("    Eps final   = %.6f", r$Eps_Final),
    sprintf("    Rhat max    = %.4f  [%s]", r$Rhat_Max, lab),
    sprintf("    ESS min     = %.0f", r$ESS_Min), "")
}

lines <- c(lines,
  strrep("-", 70), "SELECCION DE MODELO — LOO-ELPD (PSIS-LOO sobre cadenas MALA)", strrep("-", 70),
  "Mayor ELPD = mejor predictivo. Diferencia significativa si |elpd_diff| > 2 x se_diff", "",
  sprintf("  %-12s  %8s  %6s  %8s  %10s  %8s", "Modelo","ELPD","SE","LOOIC","elpd_diff","se_diff"),
  sprintf("  %s", strrep("-", 58)))
for (i in seq_len(nrow(loo_comparison_df))) {
  r <- loo_comparison_df[i, ]
  z <- if (i==1) "—" else sprintf("%.2f", abs(r$elpd_diff)/r$se_diff)
  sig <- if (i==1) "(referencia)" else if (abs(r$elpd_diff) > 2*r$se_diff) "* sig." else "  n.s."
  lines <- c(lines, sprintf("  %-12s  %8.2f  %6.2f  %8.2f  %10.2f  %8.2f  z=%s  %s",
                            toupper(r$model), r$elpd_loo, r$se_elpd_loo, r$looic, r$elpd_diff, r$se_diff, z, sig))
}
lines <- c(lines, "", sprintf("  MEJOR MODELO (LOO, MALA): %s", toupper(best_model)), "")

# BF de referencia
if (!is.null(bf_table)) {
  lines <- c(lines,
    strrep("-", 70), "BAYES FACTORS — REFERENCIA (bridge sampling de Stan)", strrep("-", 70),
    "La verosimilitud marginal es propiedad del modelo, no del sampler;",
    "se reporta una sola vez (calculada con Stan) y aplica tambien a MALA.", "",
    "  Log-verosimilitud marginal (Stan):",
    sprintf("    Standard : %.4f", logml["standard"]),
    sprintf("    Power    : %.4f", logml["power"]),
    sprintf("    Reversal : %.4f", logml["reversal"]), "")
  for (i in seq_len(nrow(bf_table))) {
    r <- bf_table[i, ]
    lines <- c(lines, sprintf("  %-24s  log_BF=%8.3f  BF=%9.2f  [%-12s]  -> %s",
                              r$comparacion, r$log_BF, r$BF, r$evidencia, r$favorece))
  }
  lines <- c(lines, "")
}

# Pareto-k
lines <- c(lines, strrep("-", 70), "DIAGNOSTICO PARETO-K", strrep("-", 70),
  "k < 0.5: bueno | 0.5-0.7: aceptable | > 0.7: LOO poco fiable", "")
for (m in MODELS) {
  if (is.null(pareto_summary[[m]])) next
  ps <- pareto_summary[[m]]
  fiab <- if (ps$pct_bad < 5) "LOO FIABLE" else if (ps$pct_bad < 10) "LOO ACEPTABLE" else "LOO POCO FIABLE"
  lines <- c(lines,
    sprintf("  %s  [%s]:", toupper(m), fiab),
    sprintf("    Bueno (k<0.5)    : %.1f%%", ps$pct_good),
    sprintf("    Aceptable (0.5-0.7): %.1f%%", ps$pct_ok),
    sprintf("    Malo/Muy malo (>0.7): %.1f%%", ps$pct_bad),
    sprintf("    ELPD = %.2f (SE = %.2f) | p_loo = %.2f | LOOIC = %.2f",
            ps$elpd, ps$se_elpd, ps$p_loo, ps$looic), "")
}

# Lambda
if (!is.null(lambda_summary) && nrow(lambda_summary) > 0) {
  lines <- c(lines, strrep("-", 70), "PARAMETRO DE ASIMETRIA LAMBDA", strrep("-", 70),
    "Lambda != 1 indica asimetria significativa respecto al logit estandar.", "")
  for (i in seq_len(nrow(lambda_summary))) {
    r <- lambda_summary[i, ]
    ev <- if (r$q2.5 > 1 | r$q97.5 < 1) "ASIMETRIA SIGNIFICATIVA (IC95% no contiene 1)" else "asimetria no significativa (IC95% contiene 1)"
    lines <- c(lines,
      sprintf("  %s:", r$modelo),
      sprintf("    Media = %.4f | SD = %.4f | IC95%% = [%.4f, %.4f]", r$mean, r$sd, r$q2.5, r$q97.5),
      sprintf("    P(lambda>1) = %.1f%% | P(lambda<1) = %.1f%%", r$p_gt1, r$p_lt1),
      sprintf("    -> %s", ev), "")
  }
}

# Estimaciones Power
lines <- c(lines, strrep("-", 70), "ESTIMACIONES POSTERIORES — POWER LOGIT (MALA)", strrep("-", 70),
  sprintf("  %-28s  %8s  %8s  %10s  %10s", "Parametro","Media","SD","IC2.5%","IC97.5%"),
  sprintf("  %s", strrep("-", 70)))
for (i in seq_len(nrow(pl_estimates))) {
  r <- pl_estimates[i, ]
  lines <- c(lines, sprintf("  %-28s  %8.4f  %8.4f  %10.4f  %10.4f",
                            r$parametro, r$mean, r$sd, r$lower_ci, r$upper_ci))
}
lines <- c(lines, "", strrep("=", 70), "FIN DEL RESUMEN (MALA)", strrep("=", 70))

writeLines(lines, here("output", "mala_results_summary.txt"))
cat("-> mala_results_summary.txt\n")

cat("\n", strrep("=", 60), "\nCOMBINACION MALA COMPLETA\n", strrep("=", 60), "\n", sep = "")
cat(sprintf("Mejor modelo (LOO, MALA): %s\n", toupper(best_model)))
print(comparison)
