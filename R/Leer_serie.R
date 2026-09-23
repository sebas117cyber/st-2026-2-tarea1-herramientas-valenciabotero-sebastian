# funsion leer serie



library(tibble)
library(dplyr)

leer_serie <- function(x, fuente = "desconocida", unidad = "sin unidad") {
  
  # fuente y unidad = desconocida si no se provee
  
  # LEER LOS DATOS SEGUN EL TIPO DE ENTRADA
  
  #  CASO objeto ts
  if (inherits(x, "ts")) {
    
    freq   <- frequency(x)        # periodos estcional por año (s = 1, 4, 12 …) 
    valores <- as.numeric(x)
    n       <- length(valores)    # numero de datos 
    
    # time(x) devuelve el año en formato decimal, ej: 1949.000, 1949.083 …
    tiempo   <- as.numeric(time(x))
    anios    <- floor(tiempo)
    fraccion <- tiempo - anios
    
    # Periodo dentro del año (1, 2, 3 … freq)
    periodo <- round(fraccion * freq) + 1L # 1L = un perido s
    periodo <- pmin(pmax(periodo, 1L), freq)   # limitar entre 1 y freq
    
    # Convertir a fecha Date segun la frecuencia
    if (freq == 12) {
      # Mensual: primer dia del mes
      fechas <- as.Date(paste(anios, periodo, "01", sep = "-"))
      
    } else if (freq == 4) {
      # Trimestral: primer dia del trimestre (mes 1, 4, 7 o 10)
      mes_inicio <- (periodo - 1L) * 3L + 1L
      fechas <- as.Date(paste(anios, mes_inicio, "01", sep = "-"))
      
    } else if (freq == 1) {
      # Anual: para el 1 de enero
      fechas <- as.Date(paste(anios, "01", "01", sep = "-"))
      
    } else if (freq == 52) {
      # Semanal: 1 de enero + (semana - 1) * 7 dias
      fechas <- as.Date(paste(anios, "01", "01", sep = "-")) + (periodo - 1L) * 7L
      
    } else {
      # Cualquier otra frecuencia: dia del año
      fechas <- as.Date(paste(anios, "01", "01", sep = "-")) + (periodo - 1L)
    }
    
    # ── CASO B: el usuario entrega la ruta a un CSV ──────────────
  } else if (is.character(x)) {
    
    if (!file.exists(x)) stop("No se encontro el archivo: ", x)
    
    raw <- read.csv(x, stringsAsFactors = FALSE)
    
    # Verificar que existan las columnas requeridas
    if (!all(c("fecha", "valor") %in% names(raw))) {
      stop("El CSV debe tener columnas llamadas 'fecha' y 'valor'.")
    }
    
    # Intentar convertir fechas (primero ISO, luego DD/MM/YYYY)
    fechas <- suppressWarnings(as.Date(raw$fecha))
    if (anyNA(fechas)) {
      fechas <- suppressWarnings(as.Date(raw$fecha, format = "%d/%m/%Y"))
    }
    if (anyNA(fechas)) stop("No se pudieron convertir las fechas. Use YYYY-MM-DD o DD/MM/YYYY.")
    
    valores <- as.numeric(raw$valor)
    n       <- length(valores)
    
    # Detectar frecuencia desde la diferencia mediana entre fechas (en dias)
    diffs_dias    <- as.numeric(diff(fechas))
    mediana_dias  <- median(diffs_dias)
    
    freq <- if (mediana_dias < 3) {
      365L        # diaria
    } else if (mediana_dias < 14) {
      52L         # semanal
    } else if (mediana_dias < 60) {
      12L         # mensual
    } else if (mediana_dias < 200) {
      4L          # trimestral
    } else {
      1L          # anual
    }
    
  } else {
    stop("`x` debe ser un objeto ts o la ruta a un CSV.")
  }
  
  # ============================================================
  # BLOQUE 2: VERIFICAR QUE LAS FECHAS SEAN CORRECTAS
  # ============================================================
  
  diffs <- as.numeric(diff(fechas))
  
  # 2a. Fechas estrictamente crecientes
  if (any(diffs <= 0)) {
    stop("Las fechas no son estrictamente crecientes. Revise los datos.")
  }
  
  # 2b. Fechas equiespaciadas (tolerancia de ±3 dias para meses/trimestres)
  espaciado_esperado <- 365.25 / freq
  tolerancia         <- if (freq >= 52) 0 else 3
  
  if (any(abs(diffs - espaciado_esperado) > tolerancia)) {
    stop(
      "Las fechas no estan equiespaciadas segun la frecuencia detectada (",
      freq, " periodos/año). ",
      "Espaciado esperado: ~", round(espaciado_esperado, 1), " dias."
    )
  }
  
  # ============================================================
  # BLOQUE 3: CONSTRUIR EL TIBBLE DE SALIDA
  # ============================================================
  
  resultado <- tibble(
    t     = seq_len(n),   # indice entero desde 1
    fecha = fechas,       # clase Date
    y     = valores       # valores numericos
  )
  
  # Adjuntar atributos informativos
  attr(resultado, "frecuencia") <- freq
  attr(resultado, "fuente")     <- fuente
  attr(resultado, "unidad")     <- unidad
  
  # Mensaje de resumen
  etiqueta <- c("1" = "Anual", "4" = "Trimestral",
                "12" = "Mensual", "52" = "Semanal", "365" = "Diaria")
  
  message(
    "Serie cargada.\n",
    "  Observaciones : ", n,                                        "\n",
    "  Desde         : ", format(min(fechas)),                      "\n",
    "  Hasta         : ", format(max(fechas)),                      "\n",
    "  Frecuencia    : ", etiqueta[as.character(freq)],
    " (", freq, " periodos/año)\n",
    "  Fuente        : ", fuente,                                   "\n",
    "  Unidad        : ", unidad
  )
  
  return(resultado)
}
