# esto es del d punto 1


ljung_box <- function(r, T_obs, m, p = 0) {
  
  # r = (residuos, valores, errores), T_obs= total de observaciones
  # m = periodo estacional p = 
  
  # Estadistico Qm
  h     <- 1:m
  Qm    <- T_obs * (T_obs + 2) * sum(r[h]^2 / (T_obs - h))
  
  # Grados de libertad
  df    <- m - p
  
  # Valor critico al 5%
  cv    <- qchisq(0.95, df) # pchisq= da un cuantil especifico de la chicuadrado con df grados de libertas
  
  # Valor p
  pval  <- 1 - pchisq(Qm, df)
  
  #resultado como tabla
  tibble::tibble(
    estadistico  = round(Qm,   4),
    df           = df,
    valor_critico_5pct = round(cv, 4),
    p_value      = pval,
    rechazo_H0   = pval < 0.05
  )
}

# prueba jarque 

jarque_bera <- function(e) {
  
  # e = errores
  
  N  <- length(e)
  me <- mean(e)
  
  # Momentos centrados
  m2 <- mean((e - me)^2)
  m3 <- mean((e - me)^3)
  m4 <- mean((e - me)^4)
  
  # Asimetria
  A  <- m3 / m2^(3/2)
  #  kurtosis
  K  <- m4 / m2^2
  
  # Estadistico JB
  JB   <- (N / 6) * (A^2 + (K - 3)^2 / 4)
  
  # Grados de libertad
  df   <- 2
  
  # Valor critico al 5%
  cv   <- qchisq(0.95, df)
  
  # Valor p
  pval <- 1 - pchisq(JB, df)
  
  tibble::tibble(
    asimetria    = round(A,    4),
    curtosis     = round(K,    4),
    estadistico  = round(JB,   4),
    df           = df,
    valor_critico_5pct = round(cv, 4),
    p_value      = round(pval, 4),
    rechazo_H0   = pval < 0.05
  )
}


durbin_watson <- function(e) {
  
  # errores
  
  n   <- length(e)
  
  # Diferencias consecutivas
                              # sumas al cudrados    
  num <- sum(diff(e)^2)       # suma de las diferencias entre las e_t
  den <- sum(e^2)             # suma de los e_t 
  
  d   <- num / den
  
  # Valor critico aproximado (region de no rechazo: 1.5 < d < 2.5)
  # Limites de Savin-White para referencia general (n grande)
  dl  <- 1.5
  du  <- 2.5
  
  tibble::tibble(
    estadistico        = round(d, 4),
    df                 = NA_real_,
    valor_critico_5pct = NA_real_,   # depende de n y k; requiere tablas
    p_value            = NA_real_,   # distribucion no estandar
    interpretacion     = dplyr::case_when(
      d < dl  ~ "Autocorrelacion positiva",
      d > (4 - dl) ~ "Autocorrelacion negativa",
      d >= du & d <= (4 - du) ~ "Sin autocorrelacion",
      TRUE    ~ "Region de indecision"
    )
  )
}
