# metodos de pronostico 

library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)




# y = datos(serie ts)


# 1 media simple 


ajusta_media <- function(y) {
  
  T_obs <- length(y)
  
  # Vector de fitted values (pronósticos un paso adelante)
  fitted <- numeric(T_obs)  
  
  # Calentamiento: el primer pronóstico es Y1
  fitted[1] <- y[1]
  
  # Actualización recursiva:
  for (t in 2:T_obs) {
    fitted[t] <- fitted[t - 1] + (1 / t) * (y[t] - fitted[t - 1])
  }
  
  # Pronóstico extramuestral: media acumulada final (igual para cualquier h)
  forecast_value <- fitted[T]
  
  # Residuos: y_t - fitted_{t-1} (error de prediccion un paso antes)
  residuo <- numeric(T_obs)
  residuo[1] <- 0  # sin error en t=1
  for (t in 2:T_obs) {
    residuo[t] <- y[t] - fitted[t - 1]
  }
  
  # Resultado como lista
  resul_ajust_med <- list(
    metodo         = "Media simple",
    fitted         = fitted,
    residuales      = residuo,
    forecast_value = forecast_value,  # Y_{T+h} para cualquier h >= 1
    T_observado    = T_obs
  )
  
  return(resul_ajust_med)
}



# 2 media movil MM


ajustar_mm <- function(y, k) { 
  
  T_obs <- length(y)
  
  if (k<=2) stop("k debe ser mayor a 2")
  
  # Fitted values y residuos
  fitted    <- rep(NA, T_obs) # rep(NA, )hace que se rellenen los faltantes con NA
  residuos <- rep(NA, T_obs)
  
  # Calentamiento: k periodos, no hay pronostico hasta t = k+1
  # A partir de t = k+1, el pronostico es el promedio de las k obs anteriores
  for (j in (k + 1):T_obs) {
    fitted[j]   <- mean(y[(j - k):(j - 1)])   # promedio de las k obs previas
    residuos[j] <- y[j] - fitted[j]            # error de prediccion
  }
  
  # Pronostico extramuestral: promedio de las ultimas k observaciones
  Extra_muesMM <- mean(y[(T_obs - k + 1):T_obs])
  
  # Resultado
  Resul_MM <- list(
    metodo         = paste("Media movil k =", k),
    k              = k,
    fitted         = fitted,
    residuales     = residuos,
    Extra_muestral = Extra_muesMM,   # Y_{T+h} para cualquier h >= 1
    T_observados   = T_obs
  )
  
  return(Resul_MM)
  
}



# 3 SES suavizamiento exponencial


ajustar_ses <- function(y, alpha) {
  
  T_obs <- length(y)
  
  # Validacion de alpha
  if (alpha <= 0 | alpha >= 1) stop("alpha debe estar entre 0 y 1")
  
  # Vectores de fitted y residuos
  fitted    <- numeric(T_obs)
  residuals <- numeric(T_obs)
  
  # Calentamiento: Yhat_2 = Y_1
  fitted[1]    <- y[1]
  residuals[1] <- 0
  
  # Forma 1 - Correccion de error: Yhat_{t+1} = Yhat_t + alpha * e_t
  for (t in 2:T_obs) {
    e_t        <- y[t] - fitted[t - 1]          # error en t
    fitted[t]  <- fitted[t - 1] + alpha * e_t   # correccion de error
    residuals[t] <- e_t
  }
  
  # Verificacion: Forma 2 - Promedio ponderado (debe ser identico)
  fitted_pp    <- numeric(T_obs)
  fitted_pp[1] <- y[1]
  
  for (t in 2:T_obs) {
    fitted_pp[t] <- alpha * y[t - 1] + (1 - alpha) * fitted_pp[t - 1]
  }
  
  # Comprobacion numerica de que ambas formas coinciden
  coinciden <- all(abs(fitted - fitted_pp) < 1e-3)
  
  # Pronostico extramuestral: alpha*Y_T + (1-alpha)*Yhat_T
  forecast_value <- alpha * y[T_obs] + (1 - alpha) * fitted[T_obs]
  
  diferencia <- mean(fitted) - mean(fitted_pp)
  
  # Resultado
  Resultado_SES <- list(
    metodo          = paste("SES alpha =", alpha),
    alpha           = alpha,
    fitted          = fitted,
    fitted_ponderado = fitted_pp,
    diferena         = diferencia,
    residuals       = residuals,
    forecast_value  = forecast_value,    # Y_{T+h} para cualquier h >= 1
    formas_coinciden = coinciden,
    T_observados     = T_obs
  )
  
  
}



# 4 doble media movil


ajustar_dmm <- function(y, k) {
  
  T_obs <- length(y)
  
  if (k<=2) stop("k debe ser mayor a 2")
  
  MM    <- rep(NA, T_obs)   # media movil simple
  DMM   <- rep(NA, T_obs)   # doble media movil
  E_hat <- rep(NA, T_obs)   # nivel estimado
  b_hat <- rep(NA, T_obs)   # pendiente estimada
  fitted    <- rep(NA, T_obs) # hace que se rellenen los faltantes con NA
  residuos <- rep(NA, T_obs)
  
  # ciclo para calcular las MM
  for (t in  (2*k-1):T_obs) {
    
    MM[t] <- mean(y[(t - k + 1):t])
    # media movil de la media movil
    DMM[t] <- mean(MM[(t - k + 1):t])
    # trayectorias y pendiente
    E_hat[t] <- 2 * MM[t] - DMM[t]
    b_hat[t] <- (2 / (k - 1)) * (MM[t] - DMM[t])
    
    # pronostico h=1
    if (t > 2*k - 1) {
      fitted[t]    <- E_hat[t - 1] + b_hat[t - 1] * 1
      residuos[t] <- y[t] - fitted[t]
    }
    
  }
  
  # Pronostico extramuestral
  forecast_fn <- function(h) E_hat[T_obs] + b_hat[T_obs] * h
  
  list(
    metodo         = paste("Doble media movil k =", k),
    k              = k,
    MM             = MM,       # trayectoria MM_t(k)
    DMM            = DMM,      # trayectoria DMM_t(k)
    E_hat          = E_hat,    # trayectoria nivel
    b_hat          = b_hat,    # trayectoria pendiente
    fitted         = fitted,
    residuals      = residuos,
    forecast_fn    = forecast_fn,   # funcion: ingresa h, devuelve pronostico
    T_observados   = T_obs
  )
  
}


# 5 regresion lineal, 6 tendencia cuadratica, 7 tendencia exponencial

ajustar_tendencia <- function(y, tipo = c("lineal", "cuadratica", "exponencial"),
                              corregir_sesgo = FALSE) {
  
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
    X_lineal    <- cbind(1, t_vect)
    beta_lineal <- solve(crossprod(X_lineal), crossprod(X_lineal, y))
    beta_lineal <- as.vector(beta_lineal)
    
    ajustados_Li <- as.vector(X_lineal %*% beta_lineal)
    residuo_lin  <- y - ajustados_Li
    
    resul   <- list(
      X            = X_lineal,
      beta         = beta_lineal,
      residuos     = residuo_lin,
      nombres_coef = c("beta0", "beta1")
    )
    # FIX 3: definir fitted y residuos para el invisible()
    fitted   <- ajustados_Li
    residuos <- residuo_lin
  }
  
  # Tendencia Cuadrática
  if (tipo == "cuadratica") {
    X_cuadra  <- cbind(1, t_vect, t_vect^2)
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
    # FIX 3: definir fitted y residuos para el invisible()
    fitted   <- ajustados_cua
    residuos <- residuo_cua
  }
  
  # Tendencia Exponencial
  if (tipo == "exponencial") {
    if (any(y <= 0)) stop("Hay valores negativos o ceros en la serie")
    
    #  matriz para exponencial
    trans_y <- log(y)
    X_exp   <- cbind(1, t_vect)
    beta_ln <- solve(crossprod(X_exp), crossprod(X_exp, trans_y)) 
    beta_ln <- as.vector(beta_ln)
    
    residuos_trans <- trans_y - as.vector(X_exp %*% beta_ln)
    sigma2_trans   <- sum(residuos_trans^2) / (T_obs - 2)
    
    # 
    beta_0or <- exp(beta_ln[1])
    beta_1or <- exp(beta_ln[2])
    
    fact_sesgo <- if (corregir_sesgo) exp(sigma2_trans / 2) else 1
    
  
    ajustados_exp <- beta_0or * beta_1or^t_vect * fact_sesgo   
    residuos_exp  <- y - ajustados_exp
    
    resul <- list(
      X            = X_exp,
      beta         = beta_ln,
      residuos     = residuos_trans,
      nombres_coef = c("beta0", "beta1")
    )
    # FIX 3: definir fitted y residuos para el invisible()
    fitted   <- ajustados_exp
    residuos <- residuos_exp
  }
  
  # Mostrar en consola
  res <- Tabla_modelo(resul)
  
  print(res$tabla, digits = 4)
  cat(sprintf("\nR²            = %.6f", res$r2))
  cat(sprintf("\nσ²            = %.6f", res$sigma2))
  cat(sprintf("\nDurbin-Watson = %.4f\n", res$dw))
  
  invisible(list(
    tipo       = tipo,
    fitted     = fitted,     # FIX 3: ahora definido en cada bloque
    residuos   = residuos,   # FIX 3: ahora definido en cada bloque
    parametros = res$tabla,
    r2         = res$r2,
    sigma2     = res$sigma2,
    dw         = res$dw
  ))
}





# 8 Holt-winters aditivo (lineal)
  
ajustar_holt <- function(y, alpha, beta) {
  
  T_obs <- length(y)
  
  # Validaciones
  if (alpha <= 0 | alpha >= 1) stop("alpha debe estar entre 0 y 1")
  if (beta  <= 0 | beta  >= 1) stop("beta debe estar entre 0 y 1")
  
  # Vectores para trayectorias
  L      <- numeric(T_obs)   # nivel
  Tend   <- numeric(T_obs)   # pendiente (T_hat)
  fitted <- numeric(T_obs)   # Yhat_t
  
  # Forma 1: Ecuaciones estandar
  # Calentamiento
  L[1]      <- y[1]
  Tend[1]   <- 0
  fitted[1] <- y[1]   # Yhat_1 = Y_1 (sin pronostico real)
  
  for (t in 2:T_obs) {
    Yhat_t  <- L[t-1] + Tend[t-1]        # pronostico para t
    fitted[t] <- Yhat_t
    L[t]    <- alpha * y[t] + (1 - alpha) * Yhat_t          # nivel
    Tend[t] <- beta * (L[t] - L[t-1]) + (1 - beta) * Tend[t-1]  # pendiente
  }
  
  residuals <- y - fitted
  
  # Forma 2: Correccion de error (verificacion)
  L2    <- numeric(T_obs)
  Tend2 <- numeric(T_obs)
  
  L2[1]    <- y[1]
  Tend2[1] <- 0
  
  for (t in 2:T_obs) {
    e_t      <- y[t] - (L2[t-1] + Tend2[t-1])     # error
    L2[t]    <- L2[t-1] + Tend2[t-1] + alpha * e_t          # nivel corregido
    Tend2[t] <- Tend2[t-1] + alpha * beta * e_t             # pendiente corregida
  }
  
  # Verificacion numerica: ambas formas deben coincidir
  coinciden_L    <- all(abs(L - L2)       < 1e-10)
  coinciden_Tend <- all(abs(Tend - Tend2) < 1e-10)
  
  # Pronostico extramuestral: Y_{T+h} = L_T + T_T * h
  forecast_fn <- function(h) L[T_obs] + Tend[T_obs] * h
  
  list(
    metodo          = paste("parametros alpha =", alpha, "beta =", beta),
    alpha           = alpha,
    beta            = beta,
    L               = L,          # trayectoria nivel
    Tend            = Tend,        # trayectoria pendiente
    fitted          = fitted,
    residuals       = residuals,
    forecast_fn     = forecast_fn, # ingresa h, devuelve pronostico
    coinciden_L     = coinciden_L,
    coinciden_Tend  = coinciden_Tend,
    T_observaciones = T_obs
  )
  
  
}
  
  
