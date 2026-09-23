# funcion de Correlograma y grafico de este

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

correlograma <- function(datos, m = NULL) {
  
  valores <- datos$valores
  
  n <- length(valores) # T observaciones
  
  # Rezagos por defectos
  if (is.null(m)) {
    m <- min(floor(n / 4), 24)
  }
  
  # calculo manual del ACF
  media <- mean(valores)
  
  rezagos <- 1:m
  
  acf_mano <- function(valores, rezagos) {
    
    denominador <- sum((valores[1:n]-media)^2)
    
    map_dbl(rezagos, function(h) {
      numerador <- sum((valores[(h+1):n]-media) * (valores[1:(n-h)]-media))
      numerador/denominador
    }) 
  }
  
  tabla_acf <- tibble(
    h = rezagos,
    acf = acf_mano(valores, rezagos)
  ) 
  
  # calculo de pacf
  pacf_obj  <- pacf(valores, lag.max = m, plot = FALSE)
  pacf_vals <- as.numeric(pacf_obj$acf)
  datos_pacf <- tibble(h = rezagos, pacf = pacf_vals)
  
  # bandas de confianza
  ci     <- 0.95
  banda_rw <- qnorm((1 + ci) / 2) / sqrt(n)          # ruido blanco
  
  # Banda secuencial (rezago h): ± 1.96/sqrt(n) * sqrt(1 + 2*sum(r_k^2, k<h))
  banda_seq <- sapply(rezagos, function(h) {
    if (h == 1) {
      banda_rw
    } else {
      banda_rw * sqrt(1 + 2 * sum(tabla_acf$acf[1:(h - 1)]^2))
    }
  })
  
  # Tabla de bandas en formato largo para geom_line
  datos_bandas <- tibble(
    h            = rezagos,
    `Ruido blanco` = banda_rw,          # constante
    `Secuencial`   = banda_seq
  ) |>
    pivot_longer(
      cols      = c("Ruido blanco", "Secuencial"),
      names_to  = "tipo_banda",
      values_to = "valor_banda"
    )
  
  #  Gráfico ACF 
  p_acf <- ggplot() +
    
    geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.6) +
    
    # Barras de la ACF
    geom_segment(
      data = tabla_acf,
      aes(x = h, xend = h, y = 0, yend = acf),
      colour    = "#1F3B57",
      linewidth = 1.4
    ) +
    
    # Banda superior
    geom_line(
      data = datos_bandas,
      aes(x = h, y = valor_banda,
          colour   = tipo_banda,
          linetype = tipo_banda),
      linewidth = 0.7
    ) +
    
    # Banda inferior
    geom_line(
      data = datos_bandas,
      aes(x = h, y = -valor_banda,
          colour   = tipo_banda,
          linetype = tipo_banda),
      linewidth = 0.7
    ) +
    
    scale_colour_manual(
      name   = "Bandas de confianza",
      values = c("Ruido blanco" = "#A50026", "Secuencial" = "#1B6B30")
    ) +
    scale_linetype_manual(
      name   = "Bandas de confianza",
      values = c("Ruido blanco" = "dashed", "Secuencial" = "solid")
    ) +
    scale_x_continuous(breaks = rezagos) +
    
    labs(
      title = "ACF muestral",
      x     = "Rezago  h",
      y     = expression(hat(rho)(h))
    ) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "top")
  
  # Gráfico PACF 
  p_pacf <- ggplot() +
    
    geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.6) +
    
    # Barras de la PACF
    geom_segment(
      data = datos_pacf,
      aes(x = h, xend = h, y = 0, yend = pacf),
      colour    = "#4D4D4D",
      linewidth = 1.4
    ) +
    
    # (solo ruido blanco para PACF)
    geom_hline(
      yintercept =  banda_rw,
      colour     = "#A50026",
      linetype   = "dashed",
      linewidth  = 0.7
    ) +
    
    # Banda inferior
    geom_hline(
      yintercept = -banda_rw,
      colour     = "#A50026",
      linetype   = "dashed",
      linewidth  = 0.7
    ) +
    
    scale_x_continuous(breaks = rezagos) +
    
    labs(
      title = "PACF muestral",
      x     = "Rezago  h",
      y     = expression(hat(phi)(h))
    ) +
    theme_minimal(base_size = 13)
  
  #  Panel combinado ACF arriba + PACF abajo
  panel <- p_acf / p_pacf
  print(panel)
 
  
  acf_stats <- acf(valores, lag.max = m, plot = FALSE)$acf[-1]   # quita rezago 0
  acf_stats  <- as.numeric(acf_stats)
  
  max_dif <- max(abs(tabla_acf$acf - acf_stats))
  
  cat("\n── Verificación ACF ──────────────────────────────────────\n")
  cat(sprintf("  Máxima diferencia absoluta |r_h manual - r_h stats|: %.2e\n",
              max_dif))
  cat(sprintf("  ¿Diferencia < 1e-12? %s\n",
              ifelse(max_dif < 1e-12, "SÍ ✓", "NO ✗")))
  cat("──────────────────────────────────────────────────────────\n\n")
  
  # Devolver invisiblemente los valores calculados
  invisible(list(
    acf    = tabla_acf$acf,
    pacf   = datos_pacf,
    bandas = datos_bandas,
    n      = n,
    m      = m
  ))
   
 
}

