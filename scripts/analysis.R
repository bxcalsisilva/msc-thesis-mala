# ==============================================================================
# MASTER SCRIPT: GENERACIÓN DE RESULTADOS DEL CAPÍTULO DE SIMULACIÓN
# ==============================================================================
# Autor: [Tu Nombre]
# Descripción: Procesa los CSVs de JAGS, Stan y MALA para Logit y Power Logit.
# Genera Tablas (LaTeX/HTML) y Gráficos (ggplot2) para la tesis.
# ==============================================================================

# 1. SETUP Y LIBRERÍAS
# ==============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(kableExtra) # Para tablas bonitas
  library(patchwork)  # Para unir gráficos
  library(ggsci)      # Colores académicos (Nature/Lancet)
})

# Configuración global de gráficos
theme_set(theme_bw(base_size = 12))
scale_fill_custom <- scale_fill_jama()
scale_color_custom <- scale_color_jama()

# Crear directorios de salida si no existen
dir.create(here("output", "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("output", "tables"), recursive = TRUE, showWarnings = FALSE)

# 2. CARGA Y PROCESAMIENTO DE DATOS
# ==============================================================================
# Función para leer y etiquetar archivos basado en su nombre
read_sim_file <- function(filename) {
  # Parsear nombre: ej. "jags_n500_p10_rho0_results.csv"
  # O "stan_power_n1000_p10_rho0_s8_results.csv"
  
  parts <- str_split(filename, "_")[[1]]
  
  algorithm <- parts[1] # jags, stan, mala
  
  # Detectar si es Power Logit o Logit Estandar
  if (parts[2] == "power") {
    model_type <- "Power Logit"
    n_val <- parse_number(parts[3])
    p_val <- parse_number(parts[4])
    rho_val <- parse_number(parts[5]) / 10 # rho9 -> 0.9, rho0 -> 0.0
    # rho7 -> 0.7
    if(parts[5] == "rho7") rho_val <- 0.7 
    if(parts[5] == "rho0") rho_val <- 0.0
    
  } else {
    model_type <- "Logit Estándar"
    n_val <- parse_number(parts[2])
    p_val <- parse_number(parts[3])
    # Ajuste manual para rho (asumiendo formato rho0, rho7, rho9)
    rho_str <- parts[4]
    if(str_detect(rho_str, "rho0")) rho_val <- 0.0
    if(str_detect(rho_str, "rho7")) rho_val <- 0.7
    if(str_detect(rho_str, "rho9")) rho_val <- 0.9
  }
  
  # Leer CSV
  df <- read_csv(here("output", filename), show_col_types = FALSE) %>%
    mutate(
      Algorithm = str_to_upper(algorithm),
      Model = model_type,
      N = n_val,
      P = p_val,
      Rho = rho_val
    )
  
  return(df)
}

# Listar todos los archivos CSV en output
files <- list.files(here("output"), pattern = "*_results.csv", full.names = FALSE)

# Filtrar solo los archivos relevantes (evitar archivos basura)
files <- files[!str_detect(files, "chains")] 

cat(">>> Procesando", length(files), "archivos de simulación...\n")

# Cargar y Unir todo en un Gran DataFrame
master_df <- map_dfr(files, read_sim_file)

# Limpieza final y Factores para orden en gráficos
master_df <- master_df %>%
  mutate(
    Algorithm = factor(Algorithm, levels = c("JAGS", "STAN", "MALA")),
    Rho_Label = paste0("Correlación: ", Rho),
    Scenario_Label = paste0("N=", N, ", P=", P)
  )

# ==============================================================================
# 3. RESULTADOS: LOGIT ESTÁNDAR (Validación y Baseline)
# ==============================================================================
# Objetivo: Mostrar que MALA funciona igual que JAGS/Stan pero comparar tiempos.

df_logit <- master_df %>% filter(Model == "Logit Estándar")

# --- TABLA 1: Resumen de Eficiencia (Logit) ---
table_eff_logit <- df_logit %>%
  group_by(Algorithm, N, P, Rho) %>%
  summarise(
    Time_Avg = mean(time_total),
    ESS_Sec_Avg = mean(bulk_ess_per_sec),
    Bias_Avg = mean(abs(bias)), # Bias absoluto promedio
    Coverage_Avg = mean(coverage),
    .groups = "drop"
  )

# Guardar CSV para Tabla LaTeX
write_csv(table_eff_logit, here("output", "tables", "table_logit_summary.csv"))

# --- GRÁFICO 1: Comparación de Eficiencia Computacional (ESS/Seg) ---
p1 <- ggplot(df_logit, aes(x = Algorithm, y = bulk_ess_per_sec, fill = Algorithm)) +
  geom_boxplot(alpha = 0.7) +
  facet_grid(P ~ Rho_Label, scales = "free_y") +
  scale_y_log10() + # Escala Logarítmica es clave porque Stan es muy rápido en ESS
  labs(
    title = "Eficiencia de Muestreo: Logit Estándar",
    subtitle = "Escala Logarítmica. Comparativa entre Nivel de Correlación y Dimensión (P)",
    y = "ESS Bulk por Segundo (log10)",
    x = NULL
  ) +
  scale_fill_custom +
  theme(legend.position = "none")

ggsave(here("output", "plots", "fig01_logit_efficiency.png"), p1, width = 10, height = 6)


# --- GRÁFICO 2: Precisión (Bias) ---
# Solo mostramos que todos están cerca de 0
p2 <- ggplot(df_logit, aes(x = Algorithm, y = bias, color = Algorithm)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  stat_summary(fun = mean, geom = "point", size = 3, position = position_dodge(width = 0.5)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, position = position_dodge(width = 0.5)) +
  facet_wrap(~Scenario_Label) +
  labs(
    title = "Validación de Exactitud: Sesgo Promedio (Bias)",
    subtitle = "Todos los algoritmos convergen a la solución correcta (Bias ~ 0)",
    y = "Sesgo Medio de Parámetros",
    x = NULL
  ) +
  scale_color_custom

ggsave(here("output", "plots", "fig02_logit_bias.png"), p2, width = 10, height = 5)


# ==============================================================================
# 4. RESULTADOS: POWER LOGIT (El núcleo de la tesis)
# ==============================================================================
# Objetivo: Stan vs MALA. Impacto de Rho. Recuperación de Lambda.

df_power <- master_df %>% filter(Model == "Power Logit")

# --- GRÁFICO 3: Recuperación del Parámetro Lambda (Boxplot) ---
# Filtramos solo el parámetro Lambda
df_lambda <- df_power %>% filter(parametro == "lambda")

p3 <- ggplot(df_lambda, aes(x = as.factor(Rho), y = mean, fill = Algorithm)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red", size = 1) + # Valor Real
  facet_wrap(~ N, labeller = label_both) +
  labs(
    title = "Recuperación del Parámetro Lambda (Valor Real = 0.8)",
    subtitle = "Comparación Stan vs MALA bajo distintas correlaciones",
    y = "Lambda Estimado (Posterior Mean)",
    x = "Nivel de Correlación (Rho)"
  ) +
  scale_fill_manual(values = c("#DF8F44FF", "#00A1D5FF")) # Colores manuales Stan/Mala

ggsave(here("output", "plots", "fig03_power_lambda_recovery.png"), p3, width = 9, height = 6)


# --- GRÁFICO 4: El "Killer Chart" - Degradación por Correlación ---
# Muestra cómo el ESS/sec cae drásticamente en MALA cuando Rho sube
df_power_summary <- df_power %>%
  group_by(Algorithm, Rho, N) %>%
  summarise(mean_eff = mean(bulk_ess_per_sec), .groups = "drop")

p4 <- ggplot(df_power_summary, aes(x = Rho, y = mean_eff, color = Algorithm, group = Algorithm)) +
  geom_line(size = 1.2) +
  geom_point(size = 4) +
  facet_wrap(~N, labeller = label_both) +
  scale_y_log10() +
  labs(
    title = "Impacto de la Correlación en la Eficiencia",
    subtitle = "Modelo Power Logit: Note la caída de desempeño en MALA vs la robustez de Stan",
    y = "Eficiencia (ESS Bulk / Seg) - Escala Log",
    x = "Correlación (Rho)"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("#DF8F44FF", "#00A1D5FF"))

ggsave(here("output", "plots", "fig04_power_correlation_impact.png"), p4, width = 8, height = 5)


# --- GRÁFICO 5: Cobertura de Intervalos (Coverage) ---
# Verificamos si los intervalos de credibilidad son válidos (debería ser ~0.95)
df_cov <- df_power %>%
  group_by(Algorithm, Rho) %>%
  summarise(Coverage = mean(coverage), .groups = "drop")

p5 <- ggplot(df_cov, aes(x = as.factor(Rho), y = Coverage, fill = Algorithm)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
  coord_cartesian(ylim = c(0.80, 1.0)) + # Zoom a la parte importante
  labs(
    title = "Validación de Cobertura (Nominal = 0.95)",
    subtitle = "Proporción de veces que el valor real cae dentro del intervalo de credibilidad",
    y = "Cobertura Promedio",
    x = "Correlación (Rho)"
  ) +
  scale_fill_manual(values = c("#DF8F44FF", "#00A1D5FF"))

ggsave(here("output", "plots", "fig05_power_coverage.png"), p5, width = 8, height = 5)

# 5. GENERACIÓN DE TABLA FINAL COMPARATIVA (LaTeX)
# ==============================================================================
# Generamos una tabla resumen final para el Power Logit
final_table <- df_power %>%
  filter(parametro == "lambda") %>%
  group_by(Algorithm, N, Rho) %>%
  summarise(
    Lambda_Mean = mean(mean),
    Lambda_SD = sd(mean), # Dispersión de las estimaciones
    Bias = mean(bias),
    RMSE = sqrt(mean(squared_error)),
    ESS_Sec = mean(bulk_ess_per_sec),
    .groups = "drop"
  ) %>%
  arrange(N, Rho, Algorithm)

# Exportar para copiar en tesis
write_csv(final_table, here("output", "tables", "table_power_final_comparison.csv"))

cat(">>> Script finalizado con éxito. Gráficos en output/plots y Tablas en output/tables.\n")

# ==============================================================================
# 6. GENERACIÓN DE ARCHIVO MAESTRO DE PROMEDIOS (SOLICITUD EXTRA)
# ==============================================================================
# Objetivo: Generar un CSV grande con el promedio de todas las estadísticas
# desagregado por cada parámetro individual de cada configuración.

detailed_summary_df <- master_df %>%
  # Agrupar por todas las variables que definen el escenario y el parámetro específico
  group_by(Model, Algorithm, N, P, Rho, parametro) %>%
  # Calcular el promedio de TODAS las columnas numéricas (mean, sd, rhat, ess, bias, time, etc.)
  summarise(
    across(where(is.numeric), \(x) mean(x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  # Ordenar para que sea fácil de leer: primero por modelo, luego escenario, luego parámetro
  arrange(Model, N, P, Rho, Algorithm, parametro)

# Guardar el archivo grande en la carpeta de tablas
write.csv(detailed_summary_df, 
          file = here("output", "tables", "master_simulation_averages.csv"),
          fileEncoding = "latin1", # <--- Fuerza el encoding Latin
          row.names = FALSE)       # Importante para no guardar índice

cat(">>> Archivo maestro de promedios generado: output/tables/master_simulation_averages.csv\n")

