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


# particion

T_air   <- nrow(airmiles_serie)
h_air   <- min(12L, floor(0.2 * T_air))   # = 4.8 ~ 4)
tramo_air_est  <- slice_head(airmiles_serie, n = T_air - h_air)
tramo_air_veri <- tail(airmiles_serie, h_air)

# optimizacion de resultado

rejilla_mm_air <- data.frame(k = 3:12)
opt_air        <- optimizar(tramo_air_est$valores, "mm", rejilla_mm_air)
k_opt_air      <- opt_air$optimo$k # = 3
cat(sprintf("\nk óptimo MM = %d  (MSE = %.2f)\n",
            k_opt_air, opt_air$optimo$mse))


p_mse_air <- graficar_mse_1d(opt_air, "k", "MSE vs k — Media Móvil (airmiles)")

# ajuste
modelo_air_mm <- ajustar_mm(tramo_air_est$valores, k = k_opt_air)

# Ingenuo
ingenuo_air_val <- rep(tail(tramo_air_est$valores, 1), h_air)

# metricas

res_air <- modelo_air_mm$residuales
MSE_air_est  <- mean(res_air^2,  na.rm = TRUE)
MAD_air_est  <- mean(abs(res_air), na.rm = TRUE)
MAPE_air_est <- mean(abs(res_air /
                           tramo_air_est$valores) * 100, na.rm = TRUE)

pron_air     <- modelo_air_mm$pronosticar(h_air)
err_val_air  <- tramo_air_veri$valores - pron_air
err_ing_air  <- tramo_air_veri$valores - ingenuo_air_val

mad_ing_est_air <- mean(abs(diff(tramo_air_est$valores)))
MASE_air_mm  <- mean(abs(err_val_air))  / mad_ing_est_air
MASE_air_ing <- mean(abs(err_ing_air))  / mad_ing_est_air

cat(sprintf("Estimación — MSE=%.2f  MAD=%.2f  MAPE=%.2f%%\n",
            MSE_air_est, MAD_air_est, MAPE_air_est))
cat(sprintf("Validación — MASE modelo=%.4f   MASE ingenuo=%.4f\n",
            MASE_air_mm, MASE_air_ing))

# validacion de error

val_air <- validar_errores(modelo_air_mm, airmiles_serie, m = 6, p = 0)

# grafico final 
modelo_air_mm$ajustados <- modelo_air_mm$Ajustados
grafico_final(airmiles_serie, modelo_air_mm, h_air,
              sprintf("airmiles: Media Móvil k=%d — ajuste y pronóstico",
                      k_opt_air),
              "ej2_airmiles_final.png")


# ejemplo 3 aleatorio 

# lectura
www_serie <- leer_serie(
  WWWusage,
  fuente = "Durbin, J. and Koopman, S. J. (2001). Time Series Analysis
             by State Space Methods. Oxford University Press.",
  unidad = "numero de usuarios"
)



# graficas
p_www <- Graficar_serie(www_serie, "Usuarios WWW por minuto (100 minutos)")
ggsave("figs/ej3_www_serie.png", p_www, width = 9, height = 4)
corr_www <- correlograma(www_serie, m = 24)

# lb
lb_www <- ljung_box(corr_www$acf, T_obs = 100, m = 24, p = 0)
print(lb_www)

# particion

T_www   <- nrow(www_serie)
h_www   <- min(12L, floor(0.2 * T_www))   # = 12, floor = 20
tramo_www_est  <- slice_head(www_serie, n = T_www - h_www)
tramo_www_veri <- tail(www_serie, h_www)

# optimo para ses
rejilla_ses_www <- data.frame(alpha = seq(0.02, 0.98, by = 0.02))
opt_www         <- optimizar(tramo_www_est$valores, "ses", rejilla_ses_www)
alpha_opt_www   <- opt_www$optimo$alpha

p_mse_www <- graficar_mse_1d(opt_www, "alpha",
                             "MSE vs alpha — SES (WWWusage)")


# ajuste 

modelo_www_ses <- ajustar_ses(tramo_www_est$valores, alpha = alpha_opt_www)

# Ingenuo
ingenuo_www_val <- rep(tail(tramo_www_est$valores, 1), h_www)

# 7.
res_www <- modelo_www_ses$residuales
MSE_www_est  <- mean(res_www^2,  na.rm = TRUE)
MAD_www_est  <- mean(abs(res_www), na.rm = TRUE)
MAPE_www_est <- mean(abs(res_www /
                           tramo_www_est$valores) * 100, na.rm = TRUE)

pron_www     <- modelo_www_ses$pronosticar(h_www)
err_val_www  <- tramo_www_veri$valores - pron_www
err_ing_www  <- tramo_www_veri$valores - ingenuo_www_val

mad_ing_est_www <- mean(abs(diff(tramo_www_est$valores)))
MASE_www_ses <- mean(abs(err_val_www))  / mad_ing_est_www
MASE_www_ing <- mean(abs(err_ing_www))  / mad_ing_est_www

# validacio 
val_www <- validar_errores(modelo_www_ses, www_serie, m = 12, p = 1)

# grafico final
modelo_www_ses$ajustados <- modelo_www_ses$Ajustados
grafico_final(www_serie, modelo_www_ses, h_www,
              sprintf("WWWusage: SES alpha=%.2f — ajuste y pronóstico",
                      alpha_opt_www),
              "ej3_www_final.png")

#ejemplo 4 doble media movil aunque parece aleatoria

lake_serie <- leer_serie(
  LakeHuron,
  fuente = "Brockwell, P. J. and Davis, R. A. (1991). Time Series and
             Forecasting Methods. 2nd ed. Springer.",
  unidad = "nivel de agua en pies"
)

# Grafico y correlograma
p_lake <- Graficar_serie(lake_serie,
                         "Nivel del Lago Hurón 1875–1972")
ggsave("figs/ej4_lake_serie.png", p_lake, width = 9, height = 4)
corr_lake <- correlograma(lake_serie, m = 20)


# lb
lb_lake <- ljung_box(corr_lake$acf, T_obs = 98, m = 20, p = 0)
print(lb_lake)


#  Particion
T_lake   <- nrow(lake_serie)
h_lake   <- min(12L, floor(0.2 * T_lake))   # = 12
tramo_lake_est  <- slice_head(lake_serie, n = T_lake - h_lake)
tramo_lake_veri <- tail(lake_serie, h_lake)

# optimizacion DMM

rejilla_dmm_lake <- data.frame(k = 3:12)
opt_lake         <- optimizar(tramo_lake_est$valores, "dmm", rejilla_dmm_lake)
k_opt_lake       <- opt_lake$optimo$k

p_mse_lake <- graficar_mse_1d(opt_lake, "k",
                              "MSE vs k — DMM (LakeHuron)")
ggsave("figs/ej4_lake_mse_k.png", p_mse_lake, width = 7, height = 4)


# Ajuste
modelo_lake_dmm <- ajustar_dmm(tramo_lake_est$valores, k = k_opt_lake)

# Ingenuo
ingenuo_lake_val <- rep(tail(tramo_lake_est$valores, 1), h_lake)


# Métricas
res_lake <- modelo_lake_dmm$residuales
MSE_lake_est  <- mean(res_lake^2,  na.rm = TRUE)
MAD_lake_est  <- mean(abs(res_lake), na.rm = TRUE)
MAPE_lake_est <- mean(abs(res_lake /
                            tramo_lake_est$valores) * 100, na.rm = TRUE)

pron_lake     <- modelo_lake_dmm$pronosticar(h_lake)
err_val_lake  <- tramo_lake_veri$valores - pron_lake
err_ing_lake  <- tramo_lake_veri$valores - ingenuo_lake_val

mad_ing_est_lake <- mean(abs(diff(tramo_lake_est$valores)))
MASE_lake_dmm  <- mean(abs(err_val_lake))  / mad_ing_est_lake
MASE_lake_ing  <- mean(abs(err_ing_lake))  / mad_ing_est_lake


# Validación de errores
val_lake <- validar_errores(modelo_lake_dmm, lake_serie, m = 10, p = 0)

# Gráfico final
modelo_lake_dmm$ajustados <- modelo_lake_dmm$Ajustados
grafico_final(lake_serie, modelo_lake_dmm, h_lake,
              sprintf("LakeHuron: DMM k=%d — ajuste y pronóstico",
                      k_opt_lake),
              "ej4_lake_final.png")




# ejemplo 5 lineal 



co2_serie <- leer_serie(
  co2,
  fuente = "Keeling, C. D. and Whorf, T. P., Scripps Institution of
             Oceanography (SIO), University of California, La Jolla.",
  unidad = "CO2 (ppm)"
)

#  Gráfico y correlograma
p_co2 <- Graficar_serie(co2_serie,
                        "Concentración atmosférica ")
ggsave("figs/ej5_co2_serie.png", p_co2, width = 9, height = 4)
corr_co2 <- correlograma(co2_serie, m = 24)

# lb

lb_co2 <- ljung_box(corr_co2$acf, T_obs = 468, m = 24, p = 0)
print(lb_co2)


#  4. Partición
T_co2   <- nrow(co2_serie)
h_co2   <- min(12, floor(0.2 * T_co2))   # = 12
tramo_co2_est  <- slice_head(co2_serie, n = T_co2 - h_co2)
tramo_co2_veri <- tail(co2_serie, h_co2)



#ajuste

modelo_co2_lin <- ajustar_tendencia(tramo_co2_est$valores, "lineal")

# Ingenuo
ingenuo_co2_val <- rep(tail(tramo_co2_est$valores, 1), h_co2)


# Métricas
res_co2 <- modelo_co2_lin$residuales
MSE_co2_est  <- mean(res_co2^2)
MAD_co2_est  <- mean(abs(res_co2))
MAPE_co2_est <- mean(abs(res_co2 / tramo_co2_est$valores) * 100)

pron_co2      <- modelo_co2_lin$pronosticar(h_co2)
err_val_co2   <- tramo_co2_veri$valores - pron_co2
err_ing_co2   <- tramo_co2_veri$valores - ingenuo_co2_val

mad_ing_est_co2 <- mean(abs(diff(tramo_co2_est$valores)))
MASE_co2_lin  <- mean(abs(err_val_co2))  / mad_ing_est_co2
MASE_co2_ing  <- mean(abs(err_ing_co2))  / mad_ing_est_co2

# Validación de errores
val_co2 <- validar_errores(modelo_co2_lin, co2_serie, m = 24, p = 2)


grafico_final(co2_serie, modelo_co2_lin, h_co2,
              "CO₂ Mauna Loa: tendencia lineal — ajuste y pronóstico",
              "ej5_co2_final.png")


# ejemplo 6 cuadratica

air_serie <- leer_serie(
  AirPassengers,
  fuente = "Box, G. E. P., Jenkins, G. M. and Reinsel, G. C. (1976).
             Time Series Analysis, Forecasting and Control. 3rd ed.
             Holden-Day. Series G.",
  unidad = "pasajeros (miles)"
)

# Grafico y correlograma
p_airp <- Graficar_serie(air_serie,
                         "Pasajeros mensuales internacionales 1949–1960")
ggsave("figs/ej6_airpass_serie.png", p_airp, width = 9, height = 4)
corr_airp <- correlograma(air_serie, m = 24)

# lb 
lb_airp <- ljung_box(corr_airp$acf, T_obs = 144, m = 24, p = 0)
print(lb_airp)


# particion


T_airp  <- nrow(air_serie)
h_airp  <- min(12L, floor(0.2 * T_airp))   # = 12
tramo_airp_est  <- slice_head(air_serie, n = T_airp - h_airp)
tramo_airp_veri <- tail(air_serie, h_airp)

# ajuste

modelo_airp_cua <- ajustar_tendencia(tramo_airp_est$valores, "cuadratica")

# Ingenuo
ingenuo_airp_val <- rep(tail(tramo_airp_est$valores, 1), h_airp)


# Métricas
res_airp <- modelo_airp_cua$residuales
MSE_airp_est  <- mean(res_airp^2)
MAD_airp_est  <- mean(abs(res_airp))
MAPE_airp_est <- mean(abs(res_airp / tramo_airp_est$valores) * 100)

pron_airp     <- modelo_airp_cua$pronosticar(h_airp)
err_val_airp  <- tramo_airp_veri$valores - pron_airp
err_ing_airp  <- tramo_airp_veri$valores - ingenuo_airp_val

mad_ing_est_airp <- mean(abs(diff(tramo_airp_est$valores)))
MASE_airp_cua  <- mean(abs(err_val_airp))  / mad_ing_est_airp
MASE_airp_ing  <- mean(abs(err_ing_airp))  / mad_ing_est_airp


# Validacion
val_airp <- validar_errores(modelo_airp_cua, air_serie, m = 12, p = 3)

# grafico final
grafico_final(air_serie, modelo_airp_cua, h_airp,
              "AirPassengers: tendencia cuadrática — ajuste y pronóstico",
              "ej6_airpass_final.png")

