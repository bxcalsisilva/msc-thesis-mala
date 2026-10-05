# scripts/plot_fig_traza_lambda.R
# Figura: Trazas MCMC de log(lambda) para los modelos asimetricos.
#   Panel A: Power Logit (PL)   Panel B: Reversal Power Logit (RPL)
#
# Lee:  output/tables/mala_lambda_trace_power.csv
#       output/tables/mala_lambda_trace_reversal.csv
#       (columnas: iteracion, cadena, lambda, modelo)
# Crea: output/plots/fig1_traza_lambda_comparativa.pdf
#
# NOTA: con la configuracion final (440k iter) AMBOS modelos convergen. El RPL
# mezcla mas lento (mayor autocorrelacion) pero ya NO diverge, por lo que el
# panel B se etiqueta como convergencia con mezcla lenta, no como fallo.
# ==============================================================================
suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(ggplot2)
})

ARCH_PL  <- here("output", "tables", "mala_lambda_trace_power.csv")
ARCH_RPL <- here("output", "tables", "mala_lambda_trace_reversal.csv")
SALIDA   <- here("output", "plots", "fig1_traza_lambda_comparativa.pdf")
dir.create(dirname(SALIDA), showWarnings = FALSE, recursive = TRUE)

# Adelgazamiento solo para graficar (evita PDFs pesados con millones de puntos)
THIN_PLOT <- 100

# ------------------------------------------------------------------------------
# 1. Datos
# ------------------------------------------------------------------------------
leer_traza <- function(archivo, etiqueta) {
  d <- read.csv(archivo)
  d <- d[seq(1, nrow(d), by = THIN_PLOT), ]      # thinning por filas
  d$loglambda <- log(d$lambda)
  d$panel     <- etiqueta
  d$cadena    <- factor(d$cadena)
  d
}

pl  <- leer_traza(ARCH_PL,  "A) Modelo Power Logit (PL) — Convergencia robusta")
rpl <- leer_traza(ARCH_RPL, "B) Modelo Reversal Power Logit (RPL) — Convergencia con mezcla lenta")

datos <- rbind(pl, rpl)

# ------------------------------------------------------------------------------
# 2. Grafico (dos paneles, una serie por cadena)
# ------------------------------------------------------------------------------
fig <- ggplot(datos, aes(x = iteracion, y = loglambda, color = cadena)) +
  geom_line(linewidth = 0.25, alpha = 0.7) +
  facet_wrap(~ panel, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c("#2C7FB8", "#41B6C4", "#7FCDBB", "#253494")) +
  labs(title = "Diagnóstico de Mezcla: Parámetro de Asimetría",
       subtitle = "Comparación de las trazas de log(lambda) entre modelos asimétricos",
       x = "Iteraciones (Post-Warmup)",
       y = expression(log(lambda)),
       color = "Cadena") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(color = "#2C7FB8"),
        legend.position = "none",
        strip.text = element_text(size = 9, face = "bold"))

ggsave(SALIDA, fig, width = 9, height = 4, device = cairo_pdf)
cat(sprintf("-> %s\n", SALIDA))
