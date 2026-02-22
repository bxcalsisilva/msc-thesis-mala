# scripts/shill_bidding_app.R
rm(list = ls())

library(here)
library(dplyr)
library(posterior)

source(here("R", "mala_driver_shill.R"))

# --- Preparación de datos ---
cat("Cargando datos shill bidding...\n")
raw_data <- read.csv(here("data", "Data_shillB.csv"))

X_raw    <- raw_data[, 4:12]
y        <- raw_data$Class
X_scaled <- scale(X_raw)
X_matrix <- cbind(Intercepto = 1, X_scaled)

shill_data <- list(
  X            = X_matrix,
  y            = y,
  p_beta       = ncol(X_matrix),
  sigma2_beta  = 100, 
  sigma2_delta = 1,   
  param_names  = colnames(X_matrix)
)

# --- Configuración MCMC ---
PARAMS <- list(
  n_iter   = 110000,
  n_warmup = 10000,
  eps_init = 0.0005
)

MODELS_TO_RUN <- c("standard", "power", "reversal")
RESULTS_LIST  <- list()

# --- Ejecución ---
for (model in MODELS_TO_RUN) {
  cat(sprintf("\nModelo: %s\n", toupper(model)))
  
  res <- run_mala_shill_parallel(
    data_list  = shill_data, 
    params     = PARAMS, 
    model_type = model, 
    n_cores    = 4
  )
  
  RESULTS_LIST[[model]] <- res
  
  # Guardar resumen parcial
  write.csv(res$resumen, 
            here("output", "tables", sprintf("shill_summary_%s.csv", model)), 
            row.names = FALSE)
  
  cat(sprintf(" Time: %.2f s | Acc: %.1f%% | ESS Bulk(L): %.0f\n", 
              res$times$total, 
              res$diagnostics$accept_rate * 100,
              ifelse(model=="standard", NA, tail(res$resumen$ess_bulk, 1))))
}

# --- Generación de tablas para tesis ---
cat("\nGenerando tablas finales...\n")

# Tabla 1: Comparativa de convergencia
comparison_table <- do.call(rbind, lapply(MODELS_TO_RUN, function(m) {
  res  <- RESULTS_LIST[[m]]
  summ <- res$resumen
  data.frame(
    Modelo   = toupper(m),
    Tiempo   = res$times$total,
    Acc_Rate = res$diagnostics$accept_rate,
    Rhat_Max = max(summ$rhat, na.rm=TRUE),
    ESS_Min  = min(summ$ess_bulk, na.rm=TRUE),
    MultiESS = ifelse(is.null(res$global_metrics$multi_ess), NA, res$global_metrics$multi_ess)
  )
}))

write.table(comparison_table, here("output", "tables", "thesis_convergence_comparison.csv"), 
            row.names = FALSE, sep = ";")

# Tabla 2: Estimaciones Power Logit
pl_estimates <- RESULTS_LIST[["power"]]$resumen %>%
  select(parametro, mean, sd, lower_ci, upper_ci, rhat, ess_bulk)

write.table(pl_estimates, here("output", "tables", "thesis_power_logit_estimates.csv"), 
            row.names = FALSE, sep = ";")

# --- Exportación para gráficos ---
cat("Exportando rds para visualización...\n")

get_lambda <- function(m_name) {
  if (m_name == "standard") return(NULL)
  mcmc <- RESULTS_LIST[[m_name]]$muestras
  do.call(rbind, lapply(1:length(mcmc), function(ch) {
    data.frame(
      Iteracion = 1:nrow(mcmc[[ch]]),
      Cadena    = as.factor(ch),
      Lambda    = as.numeric(mcmc[[ch]][, "lambda"]),
      Modelo    = toupper(m_name)
    )
  }))
}

df_lambda <- rbind(get_lambda("power"), get_lambda("reversal"))

saveRDS(df_lambda, here("output", "plots", "lambda_trace_data.rds"))
saveRDS(pl_estimates, here("output", "plots", "estimates_caterpillar_data.rds"))