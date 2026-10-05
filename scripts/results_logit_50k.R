library(stringr)
library(tidyverse)
library(here)

# 1. Definir la ruta CORRECTA relativa al RProject
# -----------------------------------------------------------------------------
# Como es un RProject, la raíz es la base. Apuntamos a la carpeta "output"
path_files <- "output" 

# Verificamos si existe la carpeta para evitar errores
if (!dir.exists(path_files)) {
  stop("¡Cuidado! No encuentro la carpeta 'output'. Verifica que estés en la raíz del RProject.")
}

files <- list.files(path = path_files, 
                    pattern = "_results\\.csv$", 
                    full.names = TRUE)

# 2. Función para leer y extraer metadatos (Idéntica a la anterior)
# -----------------------------------------------------------------------------
read_and_process <- function(filepath) {
  
  fname <- basename(filepath)
  
  # --- Extracción de Metadatos ---
  algorithm <- str_extract(fname, "^[a-z]+") 
  model_type <- ifelse(str_detect(fname, "power"), "Power Logit", "Logit Estándar")
  
  n_val <- as.numeric(str_extract(fname, "(?<=_n)[0-9]+"))
  p_val <- as.numeric(str_extract(fname, "(?<=_p)[0-9]+"))
  
  rho_raw <- as.numeric(str_extract(fname, "(?<=_rho)[0-9]+"))
  rho_val <- rho_raw / 10
  
  lambda_val <- NA
  if (model_type == "Power Logit") {
    s_raw <- as.numeric(str_extract(fname, "(?<=_s)[0-9]+"))
    lambda_val <- s_raw / 10
  }
  
  # D. DETECCIÓN DE ITERACIONES (NUEVO)
  # Lógica:
  # 1. Si el nombre tiene "_n50k_", son 50,000.
  # 2. Si es Power Logit (por defecto dijimos que eran 50k en la metodología), también 50,000.
  # 3. Si no, asumimos el estándar de 10,000.
  
  is_high_iter <- str_detect(fname, "_n50k_")
  
  iterations_count <- case_when(
    is_high_iter ~ 50000,
    model_type == "Power Logit" ~ 50000, # Asumiendo tu regla para Power Logit
    TRUE ~ 10000
  )
  
  # --- Lectura del CSV ---
  df <- read_csv(filepath, show_col_types = FALSE) %>%
    mutate(
      algorithm = toupper(algorithm),
      model     = model_type,
      n         = n_val,
      p         = p_val,
      rho       = rho_val,
      lambda_sim = lambda_val,
      iterations = iterations_count, # Nueva columna
      source_file = fname
    ) %>%
    relocate(algorithm, model, n, p, rho, lambda_sim, iterations, source_file)
  
  return(df)
}

# 3. Ejecución y Limpieza de Duplicados
# -----------------------------------------------------------------------------
# 3. Ejecución y Limpieza de Duplicados (CORREGIDO)
# -----------------------------------------------------------------------------
if (length(files) > 0) {
  cat("Procesando", length(files), "archivos...\n")
  
  raw_master_df <- map_dfr(files, read_and_process)
  
  # IMPORTANTE: Verificamos el nombre exacto de la columna de parámetros
  # Asumo que se llama 'parametro' basándome en tu prompt anterior. 
  # Si se llama 'parameter', cámbialo abajo.
  col_param <- if("parametro" %in% names(raw_master_df)) "parametro" else "parameter"
  
  cat("Usando columna de parámetro:", col_param, "\n")
  
  # LÓGICA DE LIMPIEZA CORREGIDA
  master_df <- raw_master_df %>%
    # AHORA SÍ: Agrupamos también por el parámetro para no perder datos
    group_by(algorithm, model, n, p, rho, lambda_sim, iteration, !!sym(col_param)) %>%
    
    # Prioridad: Más iteraciones primero
    arrange(desc(iterations)) %>% 
    
    # Nos quedamos con la mejor versión de ESE parámetro en ESA iteración
    slice(1) %>% 
    ungroup()
  
  cat("Carga completa. Total de filas:", nrow(master_df), "\n")
  
  # Verificación rápida: 
  # Para JAGS P=10, deberíamos tener 10 filas por iteración.
  check_filas <- master_df %>% 
    filter(algorithm == "JAGS", p == 10, iteration == 1) %>% 
    nrow()
  
  cat("Verificación de integridad: Para JAGS P=10, iteración 1, hay", check_filas, "filas (deberían ser 10).\n")
  
} else {
  cat("No se encontraron archivos.\n")
}


# 4. RESULTADOS: MODELO LOGIT ESTÁNDAR
# ==============================================================================

# 1. Filtrar datos del modelo Logit
logit_df <- master_df %>% filter(model == 'Logit Estándar')

# 2. Generación de Tabla de Convergencia (R-hat y MCSE)
# -----------------------------------------------------------------------------
# Paso A: Agregación por Réplica (Nivel 1)
# Calculamos el "peor escenario" de cada simulación individual (max R-hat)
logit_replica_stats <- logit_df %>%
  group_by(algorithm, n, p, rho, iterations, iteration) %>%
  summarise(
    # ¿Cuál fue el parámetro que peor convergió en esta réplica específica?
    max_rhat_replica = max(rhat, na.rm = TRUE),
    
    # Calidad de la estimación: Ratio MCSE / SD (promedio de todos los parámetros)
    # Si es bajo (<0.1), tenemos suficientes muestras.
    avg_mcse_ratio_replica = mean(mcse_mean / sd, na.rm = TRUE),
    
    .groups = "drop"
  )

# Paso B: Resumen Global por Escenario (Nivel 2)
# Resumimos las 100 réplicas para la tabla final
logit_convergence_table <- logit_replica_stats %>%
  group_by(algorithm, n, p, rho, iterations) %>%
  summarise(
    # Promedio del R-hat máximo (Robustez media)
    Avg_Max_Rhat = mean(max_rhat_replica, na.rm = TRUE),
    
    # El peor caso absoluto observado en las 100 réplicas (Estabilidad extrema)
    Worst_Case_Rhat = max(max_rhat_replica, na.rm = TRUE),
    
    # Tasa de Convergencia: % de réplicas donde el peor R-hat fue aceptable (< 1.05)
    Convergence_Rate_Pct = mean(max_rhat_replica < 1.05, na.rm = TRUE) * 100,
    
    # Calidad promedio de Monte Carlo
    Avg_MCSE_Ratio = mean(avg_mcse_ratio_replica, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  # Ordenar para presentación lógica
  arrange(n, p, rho, iterations, algorithm)

# 3. Visualización y Exportación
# -----------------------------------------------------------------------------
print(logit_convergence_table)

# Guardar la tabla lista para tu tesis
write_delim(logit_convergence_table, here("output", "tables", "logit_convergence_summary.csv"), delim=";")

cat(">>> Tabla de convergencia Logit generada: output/tables/logit_convergence_summary.csv\n")

# Gráfico 1: Boxplot de R-hat Máximo
# -----------------------------------------------------------------------------
plot_rhat <- logit_df %>%
  group_by(algorithm, n, p, rho, iteration) %>%
  summarise(max_rhat = max(rhat, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = algorithm, y = max_rhat, fill = algorithm)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1) +
  geom_hline(yintercept = 1.05, linetype = "dashed", color = "red") + # Umbral
  facet_grid(p ~ rho, labeller = label_both) + # Facetas por P y Rho
  scale_fill_brewer(palette = "Set2") +
  coord_cartesian(ylim = c(0.9995, 1.0605)) +
  scale_y_continuous(breaks = seq(1.00, 1.10, by = 0.02)) +
  labs(
    x = "Algoritmo",
    y = expression("Máximo " * hat(R)),
    # title = "Diagnóstico de Convergencia: Distribución del R-hat Máximo",
    title = expression("Diagnóstico de Convergencia: Distribución del " * hat(R) * " Máximo"),
    subtitle = "Línea roja discontinua indica el umbral crítico de 1.05"
  ) +
  theme_bw() +
  theme(legend.position = "none")

print(plot_rhat)
ggsave("output/plots/logit_rhat_boxplot.pdf", plot_rhat, width = 8, height = 6)

# ==============================================================================
# 4.5.2 PRECISIÓN EN LA RECUPERACIÓN DE PARÁMETROS
# ==============================================================================

# 1. Cálculo de Métricas por Réplica (Nivel 1)
# -----------------------------------------------------------------------------
# Primero colapsamos los parámetros de cada réplica en un solo número promedio
accuracy_replica <- logit_df %>%
  group_by(algorithm, n, p, rho, iterations, iteration) %>%
  summarise(
    # RMSE por réplica: Raíz del promedio de los errores al cuadrado de los P parámetros
    rmse_val = sqrt(mean(squared_error, na.rm = TRUE)),
    
    # Sesgo Promedio (Raw Bias): Para ver si está centrado en 0
    bias_val = mean(mean - beta_true, na.rm = TRUE),
    
    # Cobertura Promedio: % de parámetros que cayeron dentro del intervalo
    coverage_val = mean(coverage, na.rm = TRUE),
    
    .groups = "drop"
  )

# 2. Tabla Resumen para LaTeX (Nivel 2)
# -----------------------------------------------------------------------------
accuracy_table <- accuracy_replica %>%
  group_by(algorithm, n, p, rho, iterations) %>%
  summarise(
    Avg_RMSE = mean(rmse_val),
    SD_RMSE  = sd(rmse_val),     # Para ver la estabilidad del error
    Avg_Bias = mean(bias_val),   # Debería ser cercano a 0
    Avg_Coverage = mean(coverage_val), # Objetivo: 0.95
    .groups = "drop"
  ) %>%
  arrange(n, p, rho, algorithm)

print(accuracy_table)
# Exportar tabla
write_delim(accuracy_table, here("output", "tables", "logit_accuracy_summary.csv"), delim = ";")

# 3. Gráficos
# -----------------------------------------------------------------------------

# --- Gráfico A: Boxplot de RMSE (Precisión) ---
plot_rmse <- ggplot(accuracy_replica, aes(x = algorithm, y = rmse_val, fill = algorithm)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  facet_grid(p ~ rho, scales = "free_y", labeller = label_both) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x = "Algoritmo", 
    y = "RMSE Promedio",
    title = "Precisión: Raíz del Error Cuadrático Medio (RMSE)",
    subtitle = "Menor valor indica mayor exactitud. Nótese la paridad entre métodos."
  ) +
  theme_bw() +
  theme(legend.position = "none")

print(plot_rmse)
ggsave("output/plots/logit_rmse_boxplot.pdf", plot_rmse, width = 8, height = 6)

# --- Gráfico B: Boxplot de Sesgo (Exactitud) ---
plot_bias <- ggplot(accuracy_replica, aes(x = algorithm, y = bias_val, fill = algorithm)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + # Referencia del 0
  facet_grid(p ~ rho, scales = "free_y", labeller = label_both) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x = "Algoritmo", 
    y = "Sesgo Promedio (Bias)",
    title = "Exactitud: Distribución del Sesgo",
    subtitle = "La línea roja en 0 indica estimación insesgada."
  ) +
  theme_bw() +
  theme(legend.position = "none")

print(plot_bias)
ggsave("output/plots/logit_bias_boxplot.pdf", plot_bias, width = 8, height = 6)

cat(">>> Gráficos y tablas de precisión generados.\n")

# ==============================================================================
# 4.5.3 ANÁLISIS DE EFICIENCIA COMPUTACIONAL
# ==============================================================================

# 1. Preparación de Datos
# -----------------------------------------------------------------------------
# Aseguramos que las métricas existan. 
# Nota: Usamos 'ess_bulk' (Tamaño de Muestra Efectivo del cuerpo de la distribución)
efficiency_df <- logit_df %>%
  mutate(
    # Calcular ESS/seg si no existe
    efficiency_metric = ess_bulk / time_total,
    # Costo por muestra independiente (opcional, inverso de eficiencia)
    seconds_per_sample = time_total / ess_bulk
  )

# 2. Tabla Resumen (Nivel 2)
# -----------------------------------------------------------------------------
efficiency_table <- efficiency_df %>%
  group_by(algorithm, n, p, rho, iterations) %>%
  summarise(
    # Tiempos
    Avg_Time_Total = mean(time_total, na.rm = TRUE),
    SD_Time_Total  = sd(time_total, na.rm = TRUE),
    
    # ESS (Calidad de Mezcla)
    Avg_ESS_Bulk   = mean(ess_bulk, na.rm = TRUE),
    
    # Eficiencia (Velocidad de Mezcla)
    Avg_ESS_per_Sec = mean(efficiency_metric, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  arrange(n, p, rho, algorithm)

# Exportar tabla
write_delim(efficiency_table, here("output", "tables", "logit_efficiency_summary.csv"), delim = ";")

# 3. Gráficos
# -----------------------------------------------------------------------------

# --- Gráfico A: Eficiencia (ESS / Segundo) ---
# Este es el gráfico PRINCIPAL. Muestra "cuánta información gano por segundo".
plot_efficiency <- ggplot(efficiency_df, aes(x = algorithm, y = efficiency_metric, fill = algorithm)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  facet_grid(p ~ rho, scales = "free_y", labeller = label_both) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x = "Algoritmo", 
    y = "ESS (Bulk) por Segundo",
    title = "Eficiencia Computacional: Muestras Efectivas por Segundo",
    subtitle = "Mayor valor es mejor. Stan suele dominar por su alta eficiencia de mezcla."
  ) +
  theme_bw() +
  theme(legend.position = "none")

ggsave("output/plots/logit_efficiency_boxplot.pdf", plot_efficiency, width = 8, height = 6)


# --- Gráfico B: Tiempo Total de Ejecución ---
# Este gráfico da contexto. Muestra que MALA es rápido (o comparable) incluso con 50k.
plot_time <- ggplot(efficiency_df, aes(x = algorithm, y = time_total, fill = algorithm)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  facet_grid(p ~ rho, scales = "free_y", labeller = label_both) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x = "Algoritmo", 
    y = "Tiempo Total de Ejecución (segundos)",
    title = "Costo Computacional Bruto",
    subtitle = "MALA 50k compite en tiempo total con Stan 10k gracias a iteraciones baratas."
  ) +
  theme_bw() +
  theme(legend.position = "none")

ggsave("output/plots/logit_time_boxplot.pdf", plot_time, width = 8, height = 6)

cat(">>> Tabla y gráficos de eficiencia generados.\n")

# ==============================================================================
# 5. RESULTADOS: MODELO POWER LOGIT
# ==============================================================================

# 1. Filtrar datos del modelo Power Logit
pl_df <- master_df %>% filter(model == 'Power Logit')

# 2. Generación de Tabla de Convergencia (R-hat y MCSE)
# -----------------------------------------------------------------------------
# Paso A: Agregación por Réplica (Nivel 1)
pl_replica_stats <- pl_df %>%
  # IMPORTANTE: Añadimos 'lambda_sim' al grupo por si tienes varios escenarios de asimetría
  group_by(algorithm, n, p, rho, lambda_sim, iterations, iteration) %>%
  summarise(
    max_rhat_replica = max(rhat, na.rm = TRUE),
    avg_mcse_ratio_replica = mean(mcse_mean / sd, na.rm = TRUE),
    .groups = "drop"
  )

# Paso B: Resumen Global por Escenario (Nivel 2)
pl_convergence_table <- pl_replica_stats %>%
  group_by(algorithm, n, p, rho, lambda_sim, iterations) %>%
  summarise(
    Avg_Max_Rhat = mean(max_rhat_replica, na.rm = TRUE),
    Worst_Case_Rhat = max(max_rhat_replica, na.rm = TRUE),
    Convergence_Rate_Pct = mean(max_rhat_replica < 1.05, na.rm = TRUE) * 100,
    Avg_MCSE_Ratio = mean(avg_mcse_ratio_replica, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Ordenamos incluyendo lambda
  arrange(n, p, rho, lambda_sim, iterations, algorithm)

# 3. Visualización y Exportación de la Tabla
# -----------------------------------------------------------------------------
print(pl_convergence_table)

write_delim(pl_convergence_table, here("output", "tables", "power_logit_convergence_summary.csv"), delim=";")

cat(">>> Tabla de convergencia Power Logit generada: output/tables/power_logit_convergence_summary.csv\n")

# Gráfico 1: Boxplot de R-hat Máximo (Power Logit)
# -----------------------------------------------------------------------------
plot_rhat_pl <- pl_df %>%
  group_by(algorithm, n, p, rho, lambda_sim, iteration) %>%
  summarise(max_rhat = max(rhat, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = algorithm, y = max_rhat, fill = algorithm)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1) +
  geom_hline(yintercept = 1.05, linetype = "dashed", color = "red") + 
  
  # Facetas: Aquí mantenemos P vs Rho. 
  # NOTA: Si tienes muchos lambdas distintos mezclados, considera usar: facet_grid(p ~ lambda_sim)
  facet_grid(p ~ rho, labeller = label_both) + 
  
  scale_fill_brewer(palette = "Set2") +
  
  # Ajustamos un poco el zoom porque Power Logit suele tener valores un poco más altos
  # Si ves que cortas muchos puntos, sube el límite superior (ej. a 1.10)
  # coord_cartesian(ylim = c(0.9995, 1.08)) + 
  # scale_y_continuous(breaks = seq(1.00, 1.10, by = 0.02)) +
  
  labs(
    x = "Algoritmo",
    y = expression("Máximo " * hat(R)),
    title = "Diagnóstico de Convergencia (Power Logit): R-hat Máximo",
    subtitle = "Línea roja discontinua indica el umbral crítico de 1.05"
  ) +
  theme_bw() +
  theme(legend.position = "none")

print(plot_rhat_pl)
ggsave("output/plots/power_logit_rhat_boxplot.pdf", plot_rhat_pl, width = 8, height = 6)

# ==============================================================================
# 5. RESULTADOS: MODELO POWER LOGIT
# ==============================================================================

# 0. Preparación y Separación de Datos
# -----------------------------------------------------------------------------
# Filtramos el modelo
pl_df <- master_df %>% filter(model == 'Power Logit')

# SEPARACIÓN CLAVE: Creamos dos dataframes para analizar por separado
# 1. Betas: Coeficientes de regresión
pl_betas <- pl_df %>% filter(str_detect(parametro, "beta"))

# 2. Lambda: Parámetro de forma (Asimetría)
pl_lambda <- pl_df %>% filter(str_detect(parametro, "lambda"))

cat("Datos filtrados. Betas:", nrow(pl_betas), "| Lambdas:", nrow(pl_lambda), "\n")


# ==============================================================================
# 4.6.2 PRECISIÓN EN LOS COEFICIENTES DE REGRESIÓN (BETAS)
# ==============================================================================

# 1. Cálculo de Métricas por Réplica (Solo Betas)
# -----------------------------------------------------------------------------
accuracy_replica_beta <- pl_betas %>%
  # Agrupamos también por lambda_sim para no mezclar escenarios
  group_by(algorithm, n, p, rho, lambda_sim, iterations, iteration) %>%
  summarise(
    rmse_val = sqrt(mean(squared_error, na.rm = TRUE)),
    bias_val = mean(mean - beta_true, na.rm = TRUE),
    coverage_val = mean(coverage, na.rm = TRUE),
    .groups = "drop"
  )

# 2. Tabla Resumen Betas (Nivel 2)
# -----------------------------------------------------------------------------
accuracy_table_beta <- accuracy_replica_beta %>%
  group_by(algorithm, n, p, rho, lambda_sim, iterations) %>%
  summarise(
    Avg_RMSE = mean(rmse_val),
    SD_RMSE  = sd(rmse_val),
    Avg_Bias = mean(bias_val),
    Avg_Coverage = mean(coverage_val),
    .groups = "drop"
  ) %>%
  arrange(n, p, rho, lambda_sim, algorithm)

# Exportar tabla
write_delim(accuracy_table_beta, here("output", "tables", "power_logit_beta_accuracy.csv"), delim = ";")

# 3. Gráficos Betas
# -----------------------------------------------------------------------------
# --- Gráfico RMSE Betas ---
plot_rmse_beta <- ggplot(accuracy_replica_beta, aes(x = algorithm, y = rmse_val, fill = algorithm)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  facet_grid(p ~ rho, scales = "free_y", labeller = label_both) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x = "Algoritmo", 
    y = "RMSE Promedio (Betas)",
    title = "Precisión en Betas (Power Logit)",
    subtitle = "Comparación de error en coeficientes de regresión."
  ) +
  theme_bw() + theme(legend.position = "none")

ggsave("output/plots/power_logit_beta_rmse.pdf", plot_rmse_beta, width = 8, height = 6)

# --- Gráfico Sesgo Betas ---
plot_bias_beta <- ggplot(accuracy_replica_beta, aes(x = algorithm, y = bias_val, fill = algorithm)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_grid(p ~ rho, scales = "free_y", labeller = label_both) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x = "Algoritmo", 
    y = "Sesgo Promedio (Betas)",
    title = "Exactitud en Betas (Power Logit)",
    subtitle = "Distribución del sesgo para los coeficientes."
  ) +
  theme_bw() + theme(legend.position = "none")

ggsave("output/plots/power_logit_beta_bias.pdf", plot_bias_beta, width = 8, height = 6)


# ==============================================================================
# 4.6.3 RECUPERACIÓN DEL PARÁMETRO DE FORMA (LAMBDA)
# ==============================================================================

# 1. Cálculo de Métricas por Réplica (Solo Lambda)
# -----------------------------------------------------------------------------
accuracy_replica_lambda <- pl_lambda %>%
  group_by(algorithm, n, p, rho, lambda_sim, iterations, iteration) %>%
  summarise(
    # Para un escalar, el RMSE es la raíz del error cuadrado (abs diff si es 1 param)
    rmse_val = sqrt(mean((mean - lambda_sim)^2, na.rm = TRUE)), 
    # Sesgo directo: Estimado - Verdadero (lambda_sim)
    bias_val = mean(mean - lambda_sim, na.rm = TRUE),
    coverage_val = mean(coverage, na.rm = TRUE),
    .groups = "drop"
  )

# 2. Tabla Resumen Lambda
# -----------------------------------------------------------------------------
accuracy_table_lambda <- accuracy_replica_lambda %>%
  group_by(algorithm, n, p, rho, lambda_sim, iterations) %>%
  summarise(
    Avg_RMSE_Lambda = mean(rmse_val),
    SD_RMSE_Lambda  = sd(rmse_val),
    Avg_Bias_Lambda = mean(bias_val),
    Avg_Cov_Lambda  = mean(coverage_val),
    .groups = "drop"
  ) %>%
  arrange(n, p, rho, lambda_sim, algorithm)

# Exportar tabla
write_delim(accuracy_table_lambda, here("output", "tables", "power_logit_lambda_accuracy.csv"), delim = ";")

# 3. Gráficos Lambda
# -----------------------------------------------------------------------------
# --- Gráfico Recuperación Lambda (Scatterplot vs Valor Real o Boxplot de Sesgo) ---
# Dado que lambda suele ser fijo (ej. 0.8), un Boxplot del Bias es lo más informativo
plot_bias_lambda <- ggplot(accuracy_replica_lambda, aes(x = algorithm, y = bias_val, fill = algorithm)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_grid(p ~ rho, scales = "free_y", labeller = label_both) +
  scale_fill_brewer(palette = "Oranges") + # Color distinto para diferenciar
  labs(
    x = "Algoritmo", 
    y = expression("Sesgo en " * lambda),
    title = expression("Exactitud del Parámetro de Forma (" * lambda * ")"),
    subtitle = "Capacidad de recuperar la asimetría del modelo."
  ) +
  theme_bw() + theme(legend.position = "none")

ggsave("output/plots/power_logit_lambda_bias.pdf", plot_bias_lambda, width = 8, height = 6)


# ==============================================================================
# 4.6.4 ANÁLISIS DE EFICIENCIA COMPUTACIONAL
# ==============================================================================
# Usamos pl_df completo (todo el modelo) para tiempos, 
# pero es útil ver el ESS promedio global.

# 1. Preparación de Datos
efficiency_df_pl <- pl_df %>%
  mutate(
    efficiency_metric = ess_bulk / time_total
  )

# 2. Tabla Resumen
efficiency_table_pl <- efficiency_df_pl %>%
  group_by(algorithm, n, p, rho, lambda_sim, iterations) %>%
  summarise(
    Avg_Time_Total = mean(time_total, na.rm = TRUE),
    Avg_ESS_Bulk   = mean(ess_bulk, na.rm = TRUE),
    Avg_ESS_per_Sec = mean(efficiency_metric, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(n, p, rho, algorithm)

write_delim(efficiency_table_pl, here("output", "tables", "power_logit_efficiency_summary.csv"), delim = ";")

# 3. Gráficos
# --- Gráfico Eficiencia (ESS/seg) ---
plot_efficiency_pl <- ggplot(efficiency_df_pl, aes(x = algorithm, y = efficiency_metric, fill = algorithm)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  facet_grid(p ~ rho, scales = "free_y", labeller = label_both) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x = "Algoritmo", 
    y = "ESS (Bulk) por Segundo",
    title = "Eficiencia Computacional (Power Logit)",
    subtitle = "Muestras efectivas generadas por segundo de cómputo."
  ) +
  theme_bw() + theme(legend.position = "none")

ggsave("output/plots/power_logit_efficiency_boxplot.pdf", plot_efficiency_pl, width = 8, height = 6)

# --- Gráfico Tiempo Total ---
plot_time_pl <- ggplot(efficiency_df_pl, aes(x = algorithm, y = time_total, fill = algorithm)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  facet_grid(p ~ rho, scales = "free_y", labeller = label_both) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x = "Algoritmo", 
    y = "Tiempo Total (segundos)",
    title = "Costo Computacional (Power Logit)",
    subtitle = "Tiempo de pared para completar las cadenas."
  ) +
  theme_bw() + theme(legend.position = "none")

ggsave("output/plots/power_logit_time_boxplot.pdf", plot_time_pl, width = 8, height = 6)

cat(">>> Generación completa de resultados Power Logit (Betas, Lambda, Eficiencia).\n")

# ==============================================================================
# 4.6.5 DINÁMICA DE SINTONIZACIÓN (POWER LOGIT)
# ==============================================================================

# 1. Preparación de Datos
# -----------------------------------------------------------------------------
# Filtramos solo Power Logit y solo MALA
tuning_df_pl <- master_df %>%
  filter(model == "Power Logit", algorithm == "MALA") %>%
  select(n, p, rho, lambda_sim, iterations, iteration, 
         # Asegúrate de que estos nombres coincidan con tu dataframe original
         final_step_size = matches("final_eps"), 
         final_accept_rate = matches("accept|tasa|rate"))

# 2. Tabla Resumen (Nivel 2)
# -----------------------------------------------------------------------------
tuning_table_pl <- tuning_df_pl %>%
  group_by(n, p, rho, lambda_sim, iterations) %>%
  summarise(
    # Paso de Aprendizaje (h)
    Avg_Step_Size = mean(final_step_size, na.rm = TRUE),
    SD_Step_Size  = sd(final_step_size, na.rm = TRUE),
    
    # Tasa de Aceptación
    Avg_Accept_Rate = mean(final_accept_rate, na.rm = TRUE),
    SD_Accept_Rate  = sd(final_accept_rate, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  arrange(n, p, rho)

# Exportar tabla
write_delim(tuning_table_pl, here("output", "tables", "power_logit_tuning_summary.csv"), delim = ";")

# 3. Gráficos
# -----------------------------------------------------------------------------

# --- Gráfico A: Dinámica del Paso (Step Size) ---
plot_step_pl <- ggplot(tuning_df_pl, aes(x = as.factor(rho), y = final_step_size, fill = as.factor(p))) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_brewer(palette = "Blues", name = "Dimensión (P)") +
  labs(
    x = "Correlación (rho)",
    y = "Tamaño de Paso Promedio (h)",
    title = "Adaptación del Tamaño de Paso (MALA - Power Logit)",
    subtitle = "Reducción de h ante la asimetría y dimensionalidad."
  ) +
  theme_bw()

ggsave("output/plots/power_logit_step_size.pdf", plot_step_pl, width = 7, height = 5)

# --- Gráfico B: Estabilidad de la Tasa de Aceptación ---
plot_accept_pl <- ggplot(tuning_df_pl, aes(x = as.factor(rho), y = final_accept_rate, fill = as.factor(p))) +
  geom_boxplot(alpha = 0.7) +
  geom_hline(yintercept = 0.574, linetype = "dashed", color = "red") + # Óptimo teórico
  scale_fill_brewer(palette = "Greens", name = "Dimensión (P)") +
  labs(
    x = "Correlación (rho)",
    y = "Tasa de Aceptación Final",
    title = "Tasas de Aceptación (MALA - Power Logit)",
    subtitle = "Línea roja: Óptimo teórico asintótico (0.574)."
  ) +
  theme_bw()

ggsave("output/plots/power_logit_accept_rate.pdf", plot_accept_pl, width = 7, height = 5)

cat(">>> Generados resultados de sintonización para Power Logit.\n")
