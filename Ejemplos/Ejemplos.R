# paquetes permitidps

library(tidyverse);library(purrr);library(patchwork);library(tsibble)

# ejemplo 1

# protocolo
air_serie <- leer_serie(AirPassengers, "Box, G. E. P., Jenkins, G. M. and Reinsel, G. C. (1976) Time Series Analysis, Forecasting and Control. Third Edition. Holden-Day. Series G.", 
                        "pasageros(miles)")    # metodo lineal puede
air_grafico <- Graficar_serie(air_serie, 
                              "Monthly totals of international airline passengers, 1949 to 1960."
                              ) # tipo 4 tendencia estacional 

air_correlograma <- correlograma(air_serie)

air_ljun <- ljung_box(air_serie$valores, 144, 12)

air_h_tramo <- min(12, (0.2*144)) # por tener tendencia
print(air_val_tramo)
tramo_air <- tail(air_serie, 12)


# ajuste
modelo_air <- ajustar_tendencia(air_serie$valores, "cuadratica")
MSE_modelo_air  <- sum(modelo_air$residuales^2)
MAD_modelo_air  <- sum(abs(modelo_air$residuales)) / 144
MAPE_modelo_air <- sum(abs(modelo_air$residuales / air_serie$valores))*100 / 144



Air_pronos <- modelo_air$pronosticar(3)


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
