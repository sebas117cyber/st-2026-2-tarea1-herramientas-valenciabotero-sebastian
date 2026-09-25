# metodos de pronostico 

library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)




# y = datos(serie ts)


# 1 media simple 


# 1 media simple
ajusta_media <- function(y) {
  
  # Validaciones
  stopifnot(!anyNA(y), is.numeric(y), length(y) >= 2)
  
  T_obs <- length(y)
  
  # Vector de valores ajustados (pronósticos un paso adelante)
  fitted <- numeric(T_obs)  
  
  # Calentamiento: el primer pronóstico es Y1
  fitted[1] <- y[1]
  
  # Actualización recursiva:
  for (t in 2:T_obs) {
    fitted[t] <- fitted[t - 1] + (1 / t) * (y[t] - fitted[t - 1])
  }
  
  # Residuos: y_t - fitted_{t-1} (error de prediccion un paso antes)
  residuo <- numeric(T_obs)
  residuo[1] <- 0  # sin error en t=1
  for (t in 2:T_obs) {
    residuo[t] <- y[t] - fitted[t - 1]
  }
  
  # yhat: pronostico de Y_t hecho con info hasta t-1; calentamiento = NA
  yhat    <- c(NA, fitted[-T_obs])   # yhat[t] = fitted[t-1]
  
  # pronosticar: devuelve vector de h pronosticos extramuestrales
  media_final <- fitted[T_obs]       # Pronóstico extramuestral: media acumulada final
  pronosticar <- function(h) {
    stopifnot(is.numeric(h), length(h) == 1, h >= 1)
    rep(media_final, h)              # igual para cualquier h
  }
  
  # Resultado como lista
  resul_ajust_med <- list(
    metodo         = "Media simple",
    yhat           = yhat,                          
    pronosticar    = pronosticar,                   
    parametros     = list(media_final = media_final), 
    Ajustados      = fitted,
    residuales     = residuo,
    Extra_muestral = media_final,
    T_observado    = T_obs
  )
  
  return(resul_ajust_med)
}



# 2 media movil MM
ajustar_mm <- function(y, k) { 
  
  # Validaciones
  stopifnot(!anyNA(y), is.numeric(y), length(y) >= 2)
  stopifnot(is.numeric(k), length(k) == 1, k >= 2, k < length(y))
  
  T_obs <- length(y)
  
  if (k <= 2) stop("k debe ser mayor a 2")
  
  # valores ajustados y residuos
  fitted   <- rep(NA, T_obs) # rep(NA, )hace que se rellenen los faltantes con NA
  residuos <- rep(NA, T_obs)
  
  # Calentamiento: k periodos, no hay pronostico hasta t = k+1
  # A partir de t = k+1, el pronostico es el promedio de las k obs anteriores
  for (j in (k + 1):T_obs) {
    fitted[j]   <- mean(y[(j - k):(j - 1)])   # promedio de las k obs previas
    residuos[j] <- y[j] - fitted[j]            # error de prediccion
  }
  
  # Pronostico extramuestral: promedio de las ultimas k observaciones
  Extra_muesMM <- mean(y[(T_obs - k + 1):T_obs])
  
  # yhat: pronostico de Y_t hecho con info hasta t-1; calentamiento = NA
  yhat <- fitted                                 # ya tiene NA en calentamiento
  
  # pronosticar: devuelve vector de h pronosticos extramuestrales
  pronosticar <- function(h) {
    stopifnot(is.numeric(h), length(h) == 1, h >= 1)
    rep(Extra_muesMM, h)                         # igual para cualquier h
  }
  
  # Resultado
  Resul_MM <- list(
    metodo         = paste("Media movil k =", k),
    yhat           = yhat,                        # nuevo
    pronosticar    = pronosticar,                 # nuevo
    parametros     = list(k = k,                  # nuevo
                          MM_final = Extra_muesMM),
    k              = k,
    Ajustados      = fitted,
    residuales     = residuos,
    Extra_muestral = Extra_muesMM,
    T_observados   = T_obs
  )
  
  return(Resul_MM)
}



# 3 SES suavizamiento exponencial
ajustar_ses <- function(y, alpha) {
  
  # Validaciones
  stopifnot(!anyNA(y), is.numeric(y), length(y) >= 2)
  stopifnot(is.numeric(alpha), length(alpha) == 1, alpha > 0, alpha < 1)
  
  T_obs <- length(y)
  
  # Validacion de alpha
  if (alpha <= 0 | alpha >= 1) stop("alpha debe estar entre 0 y 1")
  
  # Vectores de fitted y residuos
  fitted    <- numeric(T_obs)
  residuals <- numeric(T_obs)
  
  # Calentamiento: Yhat_2 = Y_1
  fitted[1]    <- y[1]
  residuals[1] <- 0
  
  # Correccion de error: Yhat_{t+1} = Yhat_t + alpha * e_t
  for (t in 2:T_obs) {
    e_t           <- y[t] - fitted[t - 1]        # error en t
    fitted[t]     <- fitted[t - 1] + alpha * e_t # correccion de error
    residuals[t]  <- e_t
  }
  
  # Verificacion: Promedio ponderado (debe ser identico)
  fitted_pp    <- numeric(T_obs)
  fitted_pp[1] <- y[1]
  
  for (t in 2:T_obs) {
    fitted_pp[t] <- alpha * y[t - 1] + (1 - alpha) * fitted_pp[t - 1]
  }
  
  # Comprobacion numerica de que ambas formas coinciden
  coinciden <- all(abs(fitted - fitted_pp) < 1e-10)
  
  # Pronostico extramuestral: alpha*Y_T + (1-alpha)*Yhat_T
  forecast_value <- alpha * y[T_obs] + (1 - alpha) * fitted[T_obs]
  
  diferencia <- mean(fitted) - mean(fitted_pp)
  
  # yhat: pronostico de Y_t hecho con info hasta t-1; calentamiento = NA
  yhat <- c(NA, fitted[-T_obs])                  # yhat[t] = fitted[t-1]
  
  # pronosticar: devuelve vector de h pronosticos extramuestrales
  pronosticar <- function(h) {
    stopifnot(is.numeric(h), length(h) == 1, h >= 1)
    rep(forecast_value, h)                        # igual para cualquier h
  }
  
  # Resultado
  Resultado_SES <- list(
    metodo              = paste("SES alpha =", alpha),
    yhat                = yhat,                   # nuevo
    pronosticar         = pronosticar,            # nuevo
    parametros          = list(alpha = alpha,     # nuevo
                               nivel_final = fitted[T_obs]),
    alpha               = alpha,
    Ajustados           = fitted,
    Ajustados_ponderado = fitted_pp,
    diferencia_         = diferencia,
    residuales          = residuals,
    Extra_muestral      = forecast_value,
    formas_coinciden    = coinciden,
    T_observados        = T_obs
  )
  
  return(Resultado_SES)
}



# 4 doble media movil
ajustar_dmm <- function(y, k) {
  
  # Validaciones
  stopifnot(!anyNA(y), is.numeric(y), length(y) >= 2)
  stopifnot(is.numeric(k), length(k) == 1, k >= 2, length(y) >= 2*k - 1)
  
  T_obs <- length(y)
  
  if (k <= 2) stop("k debe ser mayor a 2")
  
  MM    <- rep(NA, T_obs)   # media movil simple
  DMM   <- rep(NA, T_obs)   # doble media movil
  E_hat <- rep(NA, T_obs)   # nivel estimado
  b_hat <- rep(NA, T_obs)   # pendiente estimada
  fitted   <- rep(NA, T_obs) # hace que se rellenen los faltantes con NA
  residuos <- rep(NA, T_obs)
  
  # ciclo para calcular las MM
  for (t in (2*k - 1):T_obs) {
    
    MM[t]  <- mean(y[(t - k + 1):t])
    # media movil de la media movil
    DMM[t] <- mean(MM[(t - k + 1):t])
    # trayectorias y pendiente
    E_hat[t] <- 2 * MM[t] - DMM[t]
    b_hat[t] <- (2 / (k - 1)) * (MM[t] - DMM[t])
    
    # pronostico h=1
    if (t > 2*k - 1) {
      fitted[t]   <- E_hat[t - 1] + b_hat[t - 1] * 1
      residuos[t] <- y[t] - fitted[t]
    }
  }
  
  # yhat: pronostico de Y_t hecho con info hasta t-1; calentamiento = NA
  yhat <- fitted                                  # ya tiene NA en calentamiento
  
  # Pronostico extramuestral
  forecast_fn <- function(h) E_hat[T_obs] + b_hat[T_obs] * h
  
  # pronosticar: devuelve vector de h pronosticos extramuestrales
  pronosticar <- function(h) {
    stopifnot(is.numeric(h), length(h) == 1, h >= 1)
    forecast_fn(1:h)                              # vector Y_{T+1},...,Y_{T+h}
  }
  
  list(
    metodo      = paste("Doble media movil k =", k),
    yhat        = yhat,                           
    pronosticar = pronosticar,                    
    parametros  = list(k      = k,                
                       E_final = E_hat[T_obs],
                       b_final = b_hat[T_obs]),
    k           = k,
    MM          = MM,       # trayectoria MM_t(k)
    DMM         = DMM,      # trayectoria DMM_t(k)
    E_hat       = E_hat,    # trayectoria nivel
    b_hat       = b_hat,    # trayectoria pendiente
    Ajustados   = fitted,
    residuales  = residuos,
    Extra_muestral = forecast_fn,
    T_observados = T_obs
  )
}



# 5 regresion lineal, 6 tendencia cuadratica, 7 tendencia exponencial
ajustar_tendencia <- function(y, tipo = c("lineal", "cuadratica", "exponencial"),
                              corregir_sesgo = FALSE) {
  
  # Validaciones
  stopifnot(!anyNA(y), is.numeric(y), length(y) >= 3)
  
  tipo   <- match.arg(tipo)
  T_obs  <- length(y)
  t_vect <- seq_len(T_obs)
  
  Tabla_modelo <- function(modelo, beta = NULL, residuos = NULL, nombres_coef = NULL) {
    if (is.list(modelo)) {
      X            <- modelo$X
      beta         <- modelo$beta
      residuos     <- modelo$residuos
      nombres_coef <- modelo$nombres_coef
    } else {
      X <- modelo
    }
    
    beta     <- as.vector(beta)
    residuos <- as.vector(residuos)
    
    n <- nrow(X)
    k <- ncol(X)
    
    # barlett
    lag    <- floor(4 * (n / 100)^(2 / 9))
    sigma2 <- sum(residuos^2) / (n - k)
    
    Xe <- X * residuos
    S  <- crossprod(Xe)
    
    for (l in seq_len(lag)) {
      if (l >= n) break
      w       <- 1 - l / (lag + 1)
      Gamma_l <- crossprod(Xe[(l + 1):n, , drop = FALSE],
                           Xe[1:(n - l), , drop = FALSE])
      S       <- S + w * (Gamma_l + t(Gamma_l))
    }
    
    XtX_inv <- solve(crossprod(X))
    V_rob   <- XtX_inv %*% S %*% XtX_inv
    
    se_ord <- sqrt(pmax(diag(sigma2 * XtX_inv), 0))
    se_rob <- sqrt(pmax(diag(V_rob), 0))
    
    t_stat <- beta / se_rob
    p_val  <- 2 * pt(abs(t_stat), df = n - k, lower.tail = FALSE)
    
    y_usado <- as.vector(X %*% beta + residuos)
    ss_tot  <- sum((y_usado - mean(y_usado))^2)
    ss_res  <- sum(residuos^2)
    r2      <- 1 - ss_res / ss_tot
    
    dw <- sum(diff(residuos)^2) / ss_res
    
    tabla <- dplyr::tibble(
      coeficiente   = nombres_coef,
      estimacion    = beta,
      se_ordinario  = se_ord,
      se_robusto    = se_rob,
      t_estadistico = t_stat,
      p_valor       = p_val
    )
    
    list(tabla = tabla, r2 = r2, sigma2 = sigma2, dw = dw)
  }
  
  # Tendencia Lineal
  if (tipo == "lineal") {
    # matriz lineal
    X_lineal    <- cbind(1, t_vect)
    beta_lineal <- solve(crossprod(X_lineal), crossprod(X_lineal, y))
    beta_lineal <- as.vector(beta_lineal)
    
    # valores ajustados y residuos
    ajustados_Li <- as.vector(X_lineal %*% beta_lineal)
    residuo_lin  <- y - ajustados_Li
    
    resul <- list(
      X            = X_lineal,
      beta         = beta_lineal,
      residuos     = residuo_lin,
      nombres_coef = c("beta0", "beta1")
    )
    
    fitted   <- ajustados_Li
    residuos <- residuo_lin
    
    # pronosticar: devuelve vector de h pronosticos extramuestrales
    pronosticar <- function(h) {
      stopifnot(is.numeric(h), length(h) == 1, h >= 1)
      t_fut <- (T_obs + 1):(T_obs + h)           # tiempos futuros
      beta_lineal[1] + beta_lineal[2] * t_fut
    }
    
    parametros <- list(beta0 = beta_lineal[1], beta1 = beta_lineal[2])
  }
  
  # Tendencia Cuadrática
  if (tipo == "cuadratica") {
    # matriz cuadratica
    X_cuadra  <- cbind(1, t_vect, t_vect^2) # toma secuencia y combina
    beta_cuad <- solve(crossprod(X_cuadra), crossprod(X_cuadra, y))
    beta_cuad <- as.vector(beta_cuad)
    
    ajustados_cua <- as.vector(X_cuadra %*% beta_cuad)
    residuo_cua   <- y - ajustados_cua
    
    resul <- list(
      X            = X_cuadra,
      beta         = beta_cuad,
      residuos     = residuo_cua,
      nombres_coef = c("beta0", "beta1", "beta2")
    )
    
    fitted   <- ajustados_cua
    residuos <- residuo_cua
    
    # pronosticar: devuelve vector de h pronosticos extramuestrales
    pronosticar <- function(h) {
      stopifnot(is.numeric(h), length(h) == 1, h >= 1)
      t_fut <- (T_obs + 1):(T_obs + h)           # tiempos futuros
      beta_cuad[1] + beta_cuad[2] * t_fut + beta_cuad[3] * t_fut^2
    }
    
    parametros <- list(beta0 = beta_cuad[1], beta1 = beta_cuad[2],
                       beta2 = beta_cuad[3])
  }
  
  # Tendencia Exponencial
  if (tipo == "exponencial") {
    if (any(y <= 0)) stop("Hay valores negativos o ceros en la serie")
    
    #  matriz para exponencial
    trans_y <- log(y)           # transformacion de los datos
    X_exp   <- cbind(1, t_vect)
    beta_ln <- solve(crossprod(X_exp), crossprod(X_exp, trans_y)) 
    beta_ln <- as.vector(beta_ln)
    
    residuos_trans <- trans_y - as.vector(X_exp %*% beta_ln)
    sigma2_trans   <- sum(residuos_trans^2) / (T_obs - 2)
    
    # parametros llevados a escala
    beta_0or <- exp(beta_ln[1])
    beta_1or <- exp(beta_ln[2])
    
    # respuesta de corregir de sesgo
    fact_sesgo <- if (corregir_sesgo) exp(sigma2_trans / 2) else 1
    
    ajustados_exp <- beta_0or * beta_1or^t_vect * fact_sesgo   
    residuos_exp  <- y - ajustados_exp
    
    resul <- list(
      X            = X_exp,
      beta         = beta_ln,
      residuos     = residuos_trans,
      nombres_coef = c("beta0", "beta1")
    )
    
    fitted   <- ajustados_exp
    residuos <- residuos_exp
    
    # pronosticar: devuelve vector de h pronosticos extramuestrales
    pronosticar <- function(h) {
      stopifnot(is.numeric(h), length(h) == 1, h >= 1)
      t_fut <- (T_obs + 1):(T_obs + h)           # tiempos futuros
      beta_0or * beta_1or^t_fut * fact_sesgo
    }
    
    parametros <- list(beta0_ln = beta_ln[1], beta1_ln = beta_ln[2],
                       beta0_or = beta_0or,   beta1_or = beta_1or)
  }
  
  # Mostrar en consola
  res <- Tabla_modelo(resul)
  
  print(res$tabla, digits = 4)
  cat(sprintf("\nR²            = %.6f", res$r2))
  cat(sprintf("\nσ²            = %.6f", res$sigma2))
  cat(sprintf("\nDurbin-Watson = %.4f\n", res$dw))
  
  # yhat: pronostico de Y_t hecho con info hasta t-1; para tendencia = ajustados
  yhat <- fitted                                  # regresion usa toda la muestra
  
  invisible(list(
    tipo        = tipo,
    yhat        = yhat,                           
    pronosticar = pronosticar,                    
    parametros  = parametros,                     
    ajustados   = fitted,     
    residuales  = residuos,   
    tabla       = res$tabla,
    r2          = res$r2,
    sigma2      = res$sigma2,
    dw          = res$dw
  ))
}



# 8 Holt-winters aditivo (lineal)
ajustar_holt <- function(y, alpha, beta) {
  
  # Validaciones
  stopifnot(!anyNA(y), is.numeric(y), length(y) >= 2)
  stopifnot(is.numeric(alpha), length(alpha) == 1, alpha > 0, alpha < 1)
  stopifnot(is.numeric(beta),  length(beta)  == 1, beta  > 0, beta  < 1)
  
  T_obs <- length(y)
  
  # Validaciones
  if (alpha <= 0 | alpha >= 1) stop("alpha debe estar entre 0 y 1")
  if (beta  <= 0 | beta  >= 1) stop("beta debe estar entre 0 y 1")
  
  # Vectores para trayectorias
  L      <- numeric(T_obs)   # nivel
  Tend   <- numeric(T_obs)   # pendiente (T_hat)
  fitted <- numeric(T_obs)   # Yhat_t
  
  # Ecuaciones estandar
  # Calentamiento
  L[1]      <- y[1]
  Tend[1]   <- 0
  fitted[1] <- y[1]   # Yhat_1 = Y_1 (sin pronostico real)
  
  for (t in 2:T_obs) {
    Yhat_t    <- L[t-1] + Tend[t-1]                          # pronostico para t
    fitted[t] <- Yhat_t
    L[t]      <- alpha * y[t] + (1 - alpha) * Yhat_t         # nivel
    Tend[t]   <- beta * (L[t] - L[t-1]) + (1 - beta) * Tend[t-1]  # pendiente
  }
  
  residuals <- y - fitted
  
  # Correccion de error 
  L2    <- numeric(T_obs)
  Tend2 <- numeric(T_obs)
  
  L2[1]    <- y[1]
  Tend2[1] <- 0
  
  for (t in 2:T_obs) {
    e_t      <- y[t] - (L2[t-1] + Tend2[t-1])               # error
    L2[t]    <- L2[t-1] + Tend2[t-1] + alpha * e_t           # nivel corregido
    Tend2[t] <- Tend2[t-1] + alpha * beta * e_t              # pendiente corregida
  }
  
  # Verificacion numerica: ambas formas deben coincidir
  coinciden_L    <- all(abs(L - L2)       < 1e-10)
  coinciden_Tend <- all(abs(Tend - Tend2) < 1e-10)
  
  # Pronostico extramuestral: Y_{T+h} = L_T + T_T * h
  forecast_fn <- function(h) L[T_obs] + Tend[T_obs] * h
  
  # yhat: pronostico de Y_t hecho con info hasta t-1; calentamiento = NA
  yhat <- c(NA, fitted[-T_obs])                  # yhat[t] = fitted[t-1]
  
  # pronosticar: devuelve vector de h pronosticos extramuestrales
  pronosticar <- function(h) {
    stopifnot(is.numeric(h), length(h) == 1, h >= 1)
    forecast_fn(1:h)                              # vector Y_{T+1},...,Y_{T+h}
  }
  
  list(
    metodo          = paste("parametros alpha =", alpha, "beta =", beta),
    yhat            = yhat,                        
    pronosticar     = pronosticar,                
    parametros      = list(alpha   = alpha,       
                           beta    = beta,
                           L_final = L[T_obs],
                           T_final = Tend[T_obs]),
    alpha           = alpha,
    beta            = beta,
    L               = L,          # trayectoria nivel
    Tend            = Tend,        # trayectoria pendiente
    ajustados       = fitted,
    residuales      = residuals,
    Extra_muestral  = forecast_fn,
    coinciden_L     = coinciden_L,
    coinciden_Tend  = coinciden_Tend,
    T_observaciones = T_obs
  )
}


# # Rejillas definidas 
# rejilla_mm   <- data.frame(k     = 2:12)
# rejilla_dmm  <- data.frame(k     = 2:12)
# rejilla_ses  <- data.frame(alpha = seq(0.02, 0.98, by = 0.02))
# 
# # Rejilla Holt: producto de alpha x beta
# alpha_seq    <- seq(0.05, 0.95, by = 0.05)
# 
# # crea un data.frame de todas las combinaciones de entre alpha y beta
# rejilla_holt <- expand.grid(alpha = alpha_seq, beta = alpha_seq) 




optimizar <- function(y, metodo, rejilla) {
  
  # funsion MSE
  calc_mse <- function(params) {
    
    # Extraer parametros segun metodo
    res <- switch(metodo,
                  
                  "mm" = {
                    k <- params[1]
                    ajustar_mm(y, k)
                  },
                  
                  "dmm" = {
                    k <- params[1]
                    ajustar_dmm(y, k)
                  },
                  
                  "ses" = {
                    alpha <- params[1]
                    ajustar_ses(y, alpha)
                  },
                  
                  "holt" = {
                    alpha <- params[1]
                    beta  <- params[2]
                    ajustar_holt(y, alpha, beta)
                  },
                  
                  stop("metodo no reconocido: use 'mm', 'dmm', 'ses' o 'holt'")
    )
    
    # MSE: promedio de residuos^2 
    residuos <- res$residuales
    mean(residuos^2, na.rm = TRUE) #ignorando NA
  }
  
  # Calcular MSE para cada fila de la rejilla
  rejilla$mse <- apply(rejilla, 1, calc_mse)
  
  # Fila con menor MSE
  optimo <- rejilla[which.min(rejilla$mse), ] # which.min = saca el indice donde esta el menor mse en rejilla
  
  list(
    rejilla = rejilla,
    optimo  = optimo
  )
}


# 
# Grafico de mse vs alpha
graficar_mse_1d <- function(resultado, nombre_param, titulo) {
  
  ggplot2::ggplot(resultado$rejilla,
                  ggplot2::aes(x = .data[[nombre_param]], y = mse)) +
    ggplot2::geom_line(color = "steelblue", linewidth = 1) +
    ggplot2::geom_point(data = resultado$optimo,
                        ggplot2::aes(x = .data[[nombre_param]], y = mse),
                        color = "red", size = 3) +
    ggplot2::labs(title = titulo,
                  x     = nombre_param,
                  y     = "MSE") +
    ggplot2::theme_minimal()
}


# grafico para holt

graficar_mse_holt <- function(resultado, titulo) {
  
  ggplot2::ggplot(resultado$rejilla,
                  ggplot2::aes(x = alpha, y = beta, fill = mse)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(option = "C") +
    ggplot2::geom_point(data    = resultado$optimo,
                        ggplot2::aes(x = alpha, y = beta),
                        color   = "red",
                        size    = 4,
                        inherit.aes = FALSE) +
    ggplot2::labs(title = titulo,
                  x     = "alpha",
                  y     = "beta",
                  fill  = "MSE") +
    ggplot2::theme_minimal()
}

# # Uso de graficos
# graficar_mse_1d(res_mm,  "k",     "MSE - Media Movil")
# graficar_mse_1d(res_dmm, "k",     "MSE - Doble Media Movil")
# graficar_mse_1d(res_ses, "alpha", "MSE - SES")
# graficar_mse_holt(res_holt,         "MSE - Holt (alpha, beta)")

