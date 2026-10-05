# scripts/run_shill_generic.R
rm(list=ls())
suppressPackageStartupMessages({
  library(here)
  library(coda)
})

source(here("R", "mala_driver_shill.R"))

# --- Configuración del Modelo ---
MODEL_TO_RUN <- "standard"  

# Preparación de matriz de diseño y respuesta
data <- read.csv(here("data", "Data_shillB.csv"))
X    <- scale(data[, c("Bidder_Tendency", "Successive_Outbidding", "Winning_Ratio", "Auction_Duration")])
X    <- cbind(Intercepto = 1, X)
y    <- data$Class

shill_data <- list(
  X            = X, 
  y            = y, 
  p_beta       = ncol(X), 
  sigma2_beta  = 100, 
  sigma2_delta = 1, 
  param_names  = colnames(X)
)

params <- list(
  n_iter   = 110000, 
  n_warmup = 10000, 
  eps_init = 0.001
)

# --- Ejecución ---
res <- run_mala_shill_parallel(
  shill_data, 
  params, 
  model_type = MODEL_TO_RUN, 
  n_cores    = 4
)

# --- Diagnósticos de salida ---
cat(sprintf("\nResultados MALA: %s\n", toupper(res$model_type)))
cat(sprintf("Tiempo: %.2f s | Aceptación: %.2f%% | Final eps: %.6f\n", 
            res$times$total, 
            res$diagnostics$accept_rate * 100, 
            res$diagnostics$final_eps))

if(!is.null(res$global_metrics$multivariate_ess)) {
  cat(sprintf("Multi-ESS: %.0f\n", res$global_metrics$multivariate_ess))
}

# Tabla de parámetros
df_summary <- res$resumen
df_summary[,-1] <- round(df_summary[,-1], 3)
print(head(df_summary, 15), row.names = FALSE)

