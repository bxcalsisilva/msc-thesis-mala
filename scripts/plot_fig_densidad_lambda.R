# scripts/plot_fig_densidad_lambda.R
# Figura: Densidad posterior del parametro de forma lambda (Power Logit).
#   Densidad rellena + linea de referencia en lambda = 1 + media e IC 95%.
#
# Lee:  output/tables/mala_lambda_trace_power.csv  (columna 'lambda')
# Crea: output/plots/fig3_densidad_lambda.pdf
#
# NOTA: con lambda ~ Uniforme(-2,2) la masa se concentra en [e^-2, e^2] ≈
# [0.135, 7.39]; la densidad ya NO presenta la cola larga hasta ~50 que aparecia
# bajo el prior normal. El eje X se ajusta automaticamente a los datos.
# ==============================================================================
suppressPackageStartupMessages({
  library(here)
  library(ggplot2)
})

ARCHIVO <- here("output", "tables", "mala_lambda_trace_power.csv")
SALIDA  <- here("output", "plots", "fig3_densidad_lambda.pdf")
dir.create(dirname(SALIDA), showWarnings = FALSE, recursive = TRUE)

COL_VERDE <- "#1B7837"

# ------------------------------------------------------------------------------
# 1. Datos y resumen
# ------------------------------------------------------------------------------
lam   <- read.csv(ARCHIVO)$lambda
media <- mean(lam)
ic    <- quantile(lam, c(0.025, 0.975))

# Banda inferior donde se dibuja el punto + IC (a 1/12 de la altura maxima)
dens   <- density(lam)
y_pico <- max(dens$y)
y_pos  <- y_pico / 12

# ------------------------------------------------------------------------------
# 2. Grafico
# ------------------------------------------------------------------------------
fig <- ggplot(data.frame(lambda = lam), aes(x = lambda)) +
  geom_density(fill = COL_VERDE, color = COL_VERDE, alpha = 0.45, linewidth = 0.6) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
  annotate("text", x = 1, y = y_pico * 0.6, label = "Logit estándar (λ=1)",
           angle = 90, vjust = -0.4, hjust = 0.5, color = "red", size = 3) +
  # Media e intervalo de credibilidad
  annotate("segment", x = ic[1], xend = ic[2], y = y_pos, yend = y_pos,
           color = "black", linewidth = 0.6) +
  annotate("segment", x = ic[1], xend = ic[1], y = y_pos * 0.6, yend = y_pos * 1.4,
           color = "black", linewidth = 0.6) +
  annotate("segment", x = ic[2], xend = ic[2], y = y_pos * 0.6, yend = y_pos * 1.4,
           color = "black", linewidth = 0.6) +
  annotate("point", x = media, y = y_pos, color = "black", size = 2.6) +
  labs(title = "Evidencia de Asimetría en el Fraude",
       subtitle = "Distribución Posterior del parámetro Lambda (Power Logit)",
       x = expression(lambda~"(Parámetro de Forma)"),
       y = "Densidad") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(color = COL_VERDE))

ggsave(SALIDA, fig, width = 7, height = 4.2, device = cairo_pdf)
cat(sprintf("-> %s\n", SALIDA))
