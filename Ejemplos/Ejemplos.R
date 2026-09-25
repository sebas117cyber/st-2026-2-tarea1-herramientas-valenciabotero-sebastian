# paquetes permitidps

library(tidyverse);library(purrr);library(patchwork);library(tsibble)


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


# ejemplo 1 # cuadratica

# protocolo
air_serie <- leer_serie(AirPassengers, "Box, G. E. P., Jenkins, G. M. and Reinsel, G. C. (1976) Time Series Analysis, Forecasting and Control. Third Edition. Holden-Day. Series G.", 
                        "pasageros(miles)")    # metodo lineal puede
air_grafico <- Graficar_serie(air_serie, 
                              "Total de pasageros mensuales, 1949 to 1960."
                              ) # tipo 4 tendencia estacional 

air_correlograma <- correlograma(air_serie)
# ljun-box serie
air_ljun <- ljung_box(air_serie$valores, 144, 12)


#tramos
air_h_tramo <- min(12, (0.2*144)) # por tener tendencia # ciclo es 12
print(air_h_tramo)

tramo_estimado_air <- air_serie |> 
  slice_head(n = nrow(air_serie) - air_h_tramo)
tramo_air_veri <- tail(air_serie, 12)



# ajuste air

ingenuo_air <- ajusta_media(air_serie$valores)

modelo_tramoES <- ajustar_tendencia(tramo_estimado_air$valores, "cuadratica")

validacio_tramoES <- validar_errores(modelo_tramoES, air_serie, m = 12, p = 3)

# la diferencia entre los pronosticados y validacion

Air_pronos <- modelo_tramoES$pronosticar(12)
print(Air_pronos)
print(tramo_air_veri$valores)



# ejemplo  2
airmiles_serie <- leer_serie(airmiles) # cuadratica posibe
Graficar_serie(airmiles_serie, "airmiles") # tendencia clara
airmiles_correlograma <- correlograma(airmiles_serie)



# ejemplo 3
co2_serie <- leer_serie(co2) # tipo 4 tendencia con estacional 
Graficar_serie(co2_serie, "co2") # metodo lineal



# ejemplo 4

lynx_serie <- leer_serie(lynx) # tipo 2 con ciclos
Graficar_serie(lynx_serie, "lynx") # dmm?



# ejemplo 5
UKgas_serie <- leer_serie(UKgas) # tipo 4
Graficar_serie(UKgas_serie, "UKgas") # parece exponencial

# ejemplo 6
sunspotyear_serie <- leer_serie(sunspot.year) # tipo 2 con ciclo
Graficar_serie(sunspotyear_serie, "sunspot.year") # holt winters




# ejemplo 7
JohnsonJohnson_serie <- leer_serie(JohnsonJohnson) # tipo 4
Graficar_serie(JohnsonJohnson_serie, "JohnsonJohnson") # parece exponencial
remove(JohnsonJohnson_serie)




# ejemplo 8
Seatbelts_data <- as_tsibble(Seatbelts[, "DriversKilled"])
datos_ts <- as.ts(Seatbelts_data)
Seatbelts_serie <- leer_serie(datos_ts)
Graficar_serie(Seatbelts_serie, "Seatbelts muertos") # estacional tipo 2

remove(Seatbelts_data)
remove(datos_ts)


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
