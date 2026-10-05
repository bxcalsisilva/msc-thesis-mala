# scripts/run_stan_fase1_seleccion.R
# FASE 1 — Seleccion de covariables (corroboracion LOO).
#
# Compara, via LOO-CV, el modelo Standard COMPLETO (9 covariables) contra el
# REDUCIDO (4 covariables: las que importan segun el IC 95% en los 3 enlaces).
# Ambos se ajustan con la MISMA configuracion congelada para que la comparacion
# sea apples-to-apples.
#
# Objetivo: validar que el reducido NO es peor que el completo (y, idealmente,
# es mejor por menor sobreajuste). Si empatan -> se elige el reducido por
# parsimonia. La X reducida queda como la matriz de diseno FINAL para el resto
# de la tesis (fases 2 y 3).
#
# Salidas:
#   output/tables/stan_seleccion_loo_comparison.csv
#   output/tables/stan_seleccion_standard_full.csv
#   output/tables/stan_seleccion_standard_reducido.csv
#   output/stan_seleccion_covariables.txt
# ==============================================================================
rm(list = ls())

suppressPackageStartupMessages({
  library(here)
  library(rstan)
  library(posterior)
  library(loo)
})

source(here("R", "stan_driver_shill.R"))
rstan_options(auto_write = TRUE)

# ------------------------------------------------------------------------------
# Configuracion congelada (identica al resto del pipeline)
# ------------------------------------------------------------------------------
PARAMS <- list(n_chains = 2, n_warmup = 1000, n_iter = 5000,
               adapt_delta = 0.90, max_treedepth = 12)
BASE_SEED <- 12345

# Covariables que SOBREVIVEN la seleccion por IC 95% (en los 3 enlaces)
COVARIABLES_REDUCIDAS <- c("Bidder_Tendency", "Successive_Outbidding",
                           "Winning_Ratio", "Auction_Duration")

for (d in c("output/tables")) dir.create(here(d), showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 1. DATOS — construir X completa (9 cov) y X reducida (4 cov)
# ==============================================================================
cat("Cargando datos...\n")
raw_data <- read.csv(here("data", "Data_shillB.csv"))
y        <- as.integer(raw_data$Class)

X_full_raw    <- raw_data[, 4:12]                       # 9 covariables
X_reduced_raw <- raw_data[, COVARIABLES_REDUCIDAS]      # 4 covariables

make_stan_data <- function(X_raw, y) {
  X_scaled <- scale(X_raw)
  X_matrix <- cbind(Intercepto = 1, X_scaled)
  list(
    N = nrow(X_matrix), P = ncol(X_matrix),
    X = X_matrix, y = y, sigma = 10,
    param_names = colnames(X_matrix)
  )
}

data_full    <- make_stan_data(X_full_raw,    y)
data_reduced <- make_stan_data(X_reduced_raw, y)

cat(sprintf("Completo : P=%d (%s)\n", data_full$P,
            paste(colnames(data_full$X)[-1], collapse = ", ")))
cat(sprintf("Reducido : P=%d (%s)\n\n", data_reduced$P,
            paste(COVARIABLES_REDUCIDAS, collapse = ", ")))

# ==============================================================================
# 2. COMPILAR (el mismo logit_shill.stan sirve para cualquier P)
# ==============================================================================
cat("Compilando logit_shill.stan...\n")
compiled <- rstan::stan_model(file = here("models", "logit_shill.stan"))

# ==============================================================================
# 3. Funcion: ajusta Standard + calcula LOO para una X dada
# ==============================================================================
fit_and_loo <- function(stan_data, etiqueta) {
  cat(sprintf("\n--- Ajustando modelo %s (P=%d) ---\n", etiqueta, stan_data$P))
  res <- run_stan_shill_parallel(
    data_list = stan_data, params = PARAMS,
    model_obj = compiled, base_seed = BASE_SEED, model_type = "standard"
  )
  if (is.null(res)) stop(sprintf("Muestreo fallido (%s).", etiqueta))

  ll      <- extract_log_lik(res$muestras, "log_lik", merge_chains = FALSE)
  r_eff   <- relative_eff(exp(ll), cores = 1)
  loo_obj <- loo(ll, r_eff = r_eff, cores = 1)

  k <- pareto_k_values(loo_obj)
  cat(sprintf("  ELPD = %.2f | p_loo = %.2f | k>0.7: %.1f%% | Rhat max=%.4f | ESS min=%.0f\n",
              loo_obj$estimates["elpd_loo", "Estimate"],
              loo_obj$estimates["p_loo",    "Estimate"],
              100 * mean(k > 0.7),
              max(res$resumen$rhat, na.rm = TRUE),
              min(res$resumen$ess_bulk, na.rm = TRUE)))

  res$muestras <- NULL; gc()
  list(loo = loo_obj, resumen = res$resumen, k = k)
}

ajuste_full    <- fit_and_loo(data_full,    "COMPLETO")
ajuste_reduced <- fit_and_loo(data_reduced, "REDUCIDO")

# ==============================================================================
# 4. COMPARACION LOO (el corazon de la fase)
# ==============================================================================
loo_list <- list(completo = ajuste_full$loo, reducido = ajuste_reduced$loo)
comp     <- loo_compare(loo_list)

cat("\n========== COMPARACION LOO (completo vs reducido) ==========\n")
print(comp)

# loo_compare ordena por mejor ELPD; la fila 2 trae elpd_diff y se_diff vs la 1
comp_df         <- as.data.frame(comp)
comp_df$modelo  <- rownames(comp_df)
mejor           <- comp_df$modelo[1]
elpd_diff_2     <- comp_df$elpd_diff[2]
se_diff_2       <- comp_df$se_diff[2]
z_val           <- if (se_diff_2 > 0) abs(elpd_diff_2) / se_diff_2 else NA

veredicto <- if (is.na(z_val)) {
  "modelos identicos"
} else if (z_val < 2) {
  sprintf("EMPATE estadistico (z=%.2f < 2) -> elegir REDUCIDO por parsimonia", z_val)
} else if (mejor == "reducido") {
  sprintf("REDUCIDO significativamente MEJOR (z=%.2f) -> elegir REDUCIDO", z_val)
} else {
  sprintf("COMPLETO significativamente mejor (z=%.2f) -> revisar la seleccion", z_val)
}

cat(sprintf("\nMejor (nominal): %s | elpd_diff=%.2f | se_diff=%.2f | z=%.2f\n",
            toupper(mejor), elpd_diff_2, se_diff_2,
            ifelse(is.na(z_val), 0, z_val)))
cat(sprintf("VEREDICTO: %s\n", veredicto))

# Comparacion de complejidad efectiva (p_loo): deberia bajar en el reducido
ploo_full <- ajuste_full$loo$estimates["p_loo", "Estimate"]
ploo_red  <- ajuste_reduced$loo$estimates["p_loo", "Estimate"]
cat(sprintf("p_loo: completo=%.2f -> reducido=%.2f (menor = menos sobreajuste)\n",
            ploo_full, ploo_red))

# ==============================================================================
# 5. GUARDAR RESULTADOS
# ==============================================================================
out_df <- data.frame(
  modelo      = c("completo", "reducido"),
  n_cov       = c(data_full$P - 1, data_reduced$P - 1),
  elpd_loo    = c(ajuste_full$loo$estimates["elpd_loo", "Estimate"],
                  ajuste_reduced$loo$estimates["elpd_loo", "Estimate"]),
  se_elpd     = c(ajuste_full$loo$estimates["elpd_loo", "SE"],
                  ajuste_reduced$loo$estimates["elpd_loo", "SE"]),
  p_loo       = c(ploo_full, ploo_red),
  looic       = c(ajuste_full$loo$estimates["looic", "Estimate"],
                  ajuste_reduced$loo$estimates["looic", "Estimate"]),
  pct_k_malos = c(100 * mean(ajuste_full$k    > 0.7),
                  100 * mean(ajuste_reduced$k > 0.7))
)
out_df[, -c(1,2)] <- round(out_df[, -c(1,2)], 2)
write.csv(out_df, here("output", "tables", "stan_seleccion_loo_comparison.csv"),
          row.names = FALSE)

write.csv(ajuste_full$resumen,
          here("output", "tables", "stan_seleccion_standard_full.csv"), row.names = FALSE)
write.csv(ajuste_reduced$resumen,
          here("output", "tables", "stan_seleccion_standard_reducido.csv"), row.names = FALSE)

# Reporte de texto
lines <- c(
  strrep("=", 70),
  "FASE 1 — SELECCION DE COVARIABLES (corroboracion LOO)",
  strrep("=", 70),
  sprintf("Generado: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  sprintf("Modelo COMPLETO : %d covariables (%s)",
          data_full$P - 1, paste(colnames(data_full$X)[-1], collapse = ", ")),
  sprintf("Modelo REDUCIDO : %d covariables (%s)",
          data_reduced$P - 1, paste(COVARIABLES_REDUCIDAS, collapse = ", ")),
  "",
  "Criterio: |elpd_diff| > 2*se_diff (z>2) => diferencia significativa.",
  "",
  sprintf("  %-10s  %5s  %9s  %7s  %7s  %9s",
          "Modelo", "nCov", "ELPD", "SE", "p_loo", "k>0.7(%)"),
  sprintf("  %s", strrep("-", 55))
)
for (i in seq_len(nrow(out_df))) {
  r <- out_df[i, ]
  lines <- c(lines, sprintf("  %-10s  %5d  %9.2f  %7.2f  %7.2f  %9.1f",
                            r$modelo, r$n_cov, r$elpd_loo, r$se_elpd, r$p_loo, r$pct_k_malos))
}
lines <- c(lines, "",
  sprintf("Mejor (nominal): %s | elpd_diff=%.2f | se_diff=%.2f | z=%.2f",
          toupper(mejor), elpd_diff_2, se_diff_2, ifelse(is.na(z_val), 0, z_val)),
  sprintf("VEREDICTO: %s", veredicto),
  "",
  strrep("=", 70))
writeLines(lines, here("output", "stan_seleccion_covariables.txt"))

cat("\n-> stan_seleccion_loo_comparison.csv\n")
cat("-> stan_seleccion_standard_full.csv / _reducido.csv\n")
cat("-> stan_seleccion_covariables.txt\n")
cat("\nFASE 1 completada. La X reducida es la matriz de diseno FINAL.\n")
