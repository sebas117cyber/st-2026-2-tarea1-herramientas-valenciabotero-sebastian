# paquetes permitidps


validar_errores <- function(modelo, serie, m = NULL, p = 0) {
  
  # Extraer componentes del modelo 
  errores <- modelo$residuales
  errores_limpios <- errores[!is.na(errores)]   # quitar NA del calentamiento
  T_error <- length(errores_limpios)
  
  # Gráfico de errores en el tiempo 
  
  # Recuperar fechas de la serie original alineadas con errores no-NA
  fechas_orig <- serie$fecha                          # columna fecha de leer_serie
  idx_validos <- which(!is.na(errores))               # posiciones sin NA
  fechas_err  <- fechas_orig[idx_validos]
  
  serie_errores <- data.frame(
    t      = seq_along(errores_limpios),
    fecha  = fechas_err,
    valores = errores_limpios
  )
  
  # Asignar atributos que espera Graficar_serie()
  attr(serie_errores, "unidad") <- "Valor error"
  
  
  Graficar_serie(serie_errores, titulo = paste("Errores en el tiempo", modelo$metodo))
  
  # Correlograma de errores
  serie_err_corr <- data.frame(
    t      = seq_along(errores_limpios),
    fecha  = fechas_err,
    valores = errores_limpios
  )
  
  # m por defecto: como en correlograma(), floor(T/4) máximo 24
  if (is.null(m)) m <- min(floor(T_error / 4), 24)
  
  resul_corr <- correlograma(serie_err_corr, m = m)
  
  # Métricas de error
  valores_orig <- serie$valores[idx_validos]
  
  MSE  <- sum(errores_limpios^2)              / T_error
  MAD  <- sum(abs(errores_limpios))           / T_error
  MAPE <- sum(abs(errores_limpios / valores_orig)) * 100 / T_error
  
  # Prueba t de media cero
  media_e <- mean(errores_limpios)
  sd_e    <- sd(errores_limpios)                       # desv. estándar muestral
  t_stat  <- media_e / (sd_e / sqrt(T_error))          # estadístico t
  gl_t    <- T_error - 1
  p_t     <- 2 * pt(abs(t_stat), df = gl_t, lower.tail = FALSE)  # bilateral
  cv_t    <- qt(0.975, df = gl_t)
  
  prueba_t <- tibble::tibble(
    media_errores      = round(media_e, 6),
    estadistico_t      = round(t_stat,  4),
    gl                 = gl_t,
    valor_critico_5pct = round(cv_t,    4),
    p_value            = round(p_t,     4),
    rechazo_H0         = p_t < 0.05       # H0: media = 0
  )
  
  # Ljung-Box (m − p grados de libertad)
  # ljung_box() recibe: r = ACF muestral, T_obs, m, p
  acf_err   <- resul_corr$acf                          #
  prueba_lb <- ljung_box(r = acf_err, T_obs = T_error, m = m, p = p)
  
  # Jarque-Bera 
  prueba_jb <- jarque_bera(e = errores_limpios)
  
  #  Durbin-Watson
  prueba_dw <- durbin_watson(e = errores_limpios)
  
  # resultados en consola
  
  cat(" VALIDACIÓN DE ERRORES —", modelo$metodo, "\n")
  
  
  cat(sprintf("  MSE  = %.4f\n", MSE))
  cat(sprintf("  MAD  = %.4f\n", MAD))
  cat(sprintf("  MAPE = %.4f %%\n", MAPE))
  
  cat("\n Prueba t (H0: media errores = 0)\n")
  print(prueba_t)
  
  cat("\n Ljung-Box (H0: no autocorrelación)\n")
  print(prueba_lb)
  
  cat("\n Jarque-Bera (H0: normalidad)\n")
  print(prueba_jb)
  
  cat("\n── Durbin-Watson (H0: sin autocorrelación)\n")
  print(prueba_dw)
  
  # Devolver resultados invisiblemente
  invisible(list(
    metricas   = list(MSE = MSE, MAD = MAD, MAPE = MAPE),
    prueba_t   = prueba_t,
    ljung_box  = prueba_lb,
    jarque_bera = prueba_jb,
    durbin_watson = prueba_dw,
    m          = m,
    p          = p,
    T_errores  = T_error
  ))
}


grafico_final <- function(serie, modelo, h_val, titulo_graf,
                          nombre_archivo) {
  
  T_total <- nrow(serie)
  n_est   <- T_total - h_val          # observaciones de estimación
  
  # Ajustados sobre tramo de estimación (sin NA de calentamiento)
  ajustados_df <- tibble(
    fecha   = serie$fecha[1:n_est],
    valores = modelo$ajustados %||% modelo$Ajustados,
    tipo    = "Ajustado"
  ) |> filter(!is.na(valores))
  
  # Pronósticos extramuestrales
  pron_vals  <- modelo$pronosticar(h_val)
  pron_df    <- tibble(
    fecha   = serie$fecha[(n_est + 1):T_total],
    valores = pron_vals,
    tipo    = "Pronóstico"
  )
  
  # Serie completa
  serie_df <- tibble(
    fecha   = serie$fecha,
    valores = serie$valores,
    tipo    = "Observado"
  )
  
  p <- ggplot() +
    geom_line(data = serie_df,
              aes(x = fecha, y = valores, colour = tipo),
              linewidth = 0.8) +
    geom_line(data = ajustados_df,
              aes(x = fecha, y = valores, colour = tipo),
              linewidth = 0.8, linetype = "dashed") +
    geom_line(data = pron_df,
              aes(x = fecha, y = valores, colour = tipo),
              linewidth = 1.0) +
    geom_point(data = pron_df,
               aes(x = fecha, y = valores, colour = tipo),
               size = 2) +
    scale_colour_manual(
      values = c("Observado"   = "#00B2EE",
                 "Ajustado"    = "#FF8C00",
                 "Pronóstico"  = "#CC0000")
    ) +
    labs(title  = titulo_graf,
         x      = "Fecha",
         y      = attr(serie, "unidad"),
         colour = NULL) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "top")
  
  ggsave(file.path("figs", nombre_archivo), p,
         width = 10, height = 5)
  print(p)
}

calc_mase <- function(errores_modelo, errores_ingenuo) {
  mean(abs(errores_modelo), na.rm = TRUE) /
    mean(abs(errores_ingenuo), na.rm = TRUE)
}


# ejemplo 1 media simple
nile_serie <- leer_serie(
  Nile,
  fuente = "Durbin, J. and Koopman, S. J. (2001). Time Series Analysis
             by State Space Methods. Oxford University Press.",
  unidad = "flujo de agua (10^8 m³)"
)
p_nile <- Graficar_serie(nile_serie, "Caudal anual del río Nilo 1871–1970")

# correlograma
corr_nile <- correlograma(nile_serie, m = 20)
# lb
lb_nile <- ljung_box(corr_nile$acf, T_obs = 100, m = 20, p = 0)
print(lb_nile)

# particion

T_nile      <- nrow(nile_serie)
h_nile      <- min(12L, floor(0.2 * T_nile))   # = 20


tramo_nile_est  <- slice_head(nile_serie, n = T_nile - h_nile)
tramo_nile_veri <- tail(nile_serie, h_nile)


# ajuste 
modelo_nile <- ajusta_media(tramo_nile_est$valores)

# ingenuo

ingenuo_nile_val <- rep(tramo_nile_est$valores[nrow(tramo_nile_est)],
                        h_nile)

# metricas estimacion para el tramo

res_est_nile <- modelo_nile$residuales
MSE_nile_est  <- mean(res_est_nile^2,  na.rm = TRUE)
MAD_nile_est  <- mean(abs(res_est_nile), na.rm = TRUE)
MAPE_nile_est <- mean(abs(res_est_nile /
                            tramo_nile_est$valores) * 100, na.rm = TRUE)
cat(sprintf("\nMétricas estimación — MSE=%.2f  MAD=%.2f  MAPE=%.2f%%\n",
            MSE_nile_est, MAD_nile_est, MAPE_nile_est))


# pronosticos y validacion

pron_nile   <- modelo_nile$pronosticar(h_nile)
err_val_nile <- tramo_nile_veri$valores - pron_nile
err_ing_nile <- tramo_nile_veri$valores - ingenuo_nile_val

MSE_val_nile  <- mean(err_val_nile^2)
MAD_val_nile  <- mean(abs(err_val_nile))
MAPE_val_nile <- mean(abs(err_val_nile / tramo_nile_veri$valores) * 100)

# mase ingenuo de estimacion

ingenuo_est_nile <- rep(tramo_nile_est$valores[1],
                        length(tramo_nile_est$valores))
mad_ing_est_nile <- mean(abs(diff(tramo_nile_est$valores)))
MASE_nile  <- MAD_val_nile  / mad_ing_est_nile
MASE_ing   <- mean(abs(err_ing_nile)) / mad_ing_est_nile

cat(sprintf("Métricas validación  — MSE=%.2f  MAD=%.2f  MAPE=%.2f%%\n",
            MSE_val_nile, MAD_val_nile, MAPE_val_nile))
cat(sprintf("MASE modelo=%.4f   MASE ingenuo=%.4f\n", MASE_nile, MASE_ing))

val_nile <- validar_errores(modelo_nile, nile_serie, m = 10, p = 0)
# ajuste air


# graficos finales 
modelo_nile_gf        <- modelo_nile
modelo_nile_gf$ajustados <- modelo_nile$Ajustados
grafico_final(nile_serie, modelo_nile_gf, h_nile,
              "Nile: Media simple — ajuste y pronóstico",
              "ej1_nile_final.png")




# ejemplo 2 media movil airmiles

airmiles_serie <- leer_serie(
  airmiles,
  fuente = "F.A.A. Statistical Handbook of Aviation.",
  unidad = "millas aereas (millones)"
)


# Grafico 

p_air <- Graficar_serie(airmiles_serie,
                        "Millas de pasajeros en vuelos comerciales 1937–1960")
ggsave("figs/ej2_airmiles_serie.png", p_air, width = 9, height = 4)
corr_air <- correlograma(airmiles_serie, m = min(floor(24/4), 24))


# lb

lb_air <- ljung_box(corr_air$acf, T_obs = 24, m = 6, p = 0)
print(lb_air)




# nile_serie <- leer_serie(Nile)
# Graficar_serie(nile_serie, "Nile")  # aleatorio con cambio estructural
# remove(nile_serie)

# USAccDeaths_serie <- leer_serie(USAccDeaths)  # tipo 2
# Graficar_serie(USAccDeaths_serie, "USAccDeaths") # mm o holt
# remove(USAccDeaths_serie)


# LakeHuron_serie <- leer_serie(LakeHuron) # aleatorio con cambio estructural
# Graficar_serie(LakeHuron_serie, "LakeHuron")
# remove(LakeHuron_serie)

# WWWusage_serie <- leer_serie(WWWusage) # rara aleatoria?
# Graficar_serie(WWWusage_serie, "WWWusage") # esta en dias
# remove(WWWusage_serie)


# austres_serie <- leer_serie(austres) # tipo 3 tendencia tan sapa
# Graficar_serie(austres_serie, "austres")# lineal de lejos
# remove(austres_serie)

# UKDriverDeaths_serie <- leer_serie(UKDriverDeaths) # tipo 2 estacional con cambio estructural
# Graficar_serie(UKDriverDeaths_serie, "UKDriverDeaths") 
# remove(UKDriverDeaths_serie)


# nottem_serie <- leer_serie(nottem) # tipo 1 estacional en media
# Graficar_serie(nottem_serie, "nottem") # le sirve los 4 metodos 
# remove(nottem_serie)

# discoveries_serie <- leer_serie(discoveries) # aleatoria
# Graficar_serie(discoveries_serie, "discoveries")
# remove(discoveries_serie)
