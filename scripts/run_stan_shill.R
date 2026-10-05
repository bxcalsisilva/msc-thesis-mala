# File: scripts/run_stan_shill.R
rm(list = ls())

library(here)
library(dplyr)
library(rstan)

# Cargar el driver
source(here("R", "stan_driver_shill.R"))

# ==============================================================================
# 1. PREPARACIÓN DE DATOS (LA NUEVA ESTRUCTURA)
# ==============================================================================
cat("Cargando datos Shill Bidding...\n")
raw_data <- read.csv(here("data", "Data_shillB.csv"))

# Seleccionamos SOLO las 4 variables significativas
X_raw <- raw_data[, c("Bidder_Tendency", "Successive_Outbidding", "Winning_Ratio", "Auction_Duration")]
y     <- raw_data$Class

# Estandarizamos (Obligatorio para buena convergencia en HMC)
X_scaled <- scale(X_raw)

# Añadimos Intercepto
X_matrix <- cbind(Intercepto = 1, X_scaled)

# Generamos la lista de datos para Stan
stan_data <- list(
  N = nrow(X_matrix),
  P = ncol(X_matrix),
  X = X_matrix,
  y = as.vector(y),
  sigma = 10, # Equivalente a varianza=100
  param_names = colnames(X_matrix)
)

# ==============================================================================
# 2. CONFIGURACIÓN DEL MCMC
# ==============================================================================
PARAMS <- list(
  n_chains = 4,
  n_warmup = 1000,
  n_iter   = 5000 
)

# ==============================================================================
# 3. COMPILACIÓN Y EJECUCIÓN
# ==============================================================================
cat("Compilando el modelo Stan en C++ (Esto tomará un minuto)...\n")
compiled_model <- rstan::stan_model(file = here("models/power_logit_shill.stan"))

# set.seed(42) # Semilla maestra

res_stan <- run_stan_power_shill_parallel(
  data_list = stan_data, 
  params    = PARAMS, 
  model_obj = compiled_model, 
  base_seed = 12345
)

# ==============================================================================
# 4. EXPORTACIÓN DE RESULTADOS
# ==============================================================================
if (!is.null(res_stan$resumen)) {
  
  cat("RESULTADOS STAN (POWER LOGIT) \n")
  print(res_stan$resumen)
  
  cat(sprintf("\n Tiempo Total: %.2f segundos \n", res_stan$times$total))
  
  # Guardar CSV
  dir.create(here("output", "tables"), showWarnings = FALSE, recursive = TRUE)
  write.csv(res_stan$resumen, 
            here("output", "tables", "stan_shill_power_summary.csv"), 
            row.names = FALSE)
  
  # Guardar objeto Muestras (opcional para Traceplots)
  dir.create(here("output", "chains"), showWarnings = FALSE, recursive = TRUE)
  saveRDS(res_stan$muestras, 
          file = here("output", "chains", "stan_shill_power_fit.rds"))
}
