# paquetes permitidps



# datos ejemplo 1

datos_1 <- AirPassengers
show(datos_1)

library(tidyverse); library(tsibble); library(feasts); library(fable)
serie <- as_tsibble(AirPassengers) # 1. indice temporal explicito
serie |> autoplot(value) # 2. SIEMPRE graficar primero
serie |> gg_season(value) # 3. mirar la estacionalidad
serie |> model(STL(log(value))) |> components() |> autoplot() # 4. anatomia
serie |> ACF(value, lag_max = 36) |> autoplot() # 5. la firma de la dependencia
entrena <- serie |> filter_index(~ "1958 dic") # 6. particion sin mirar el futuro
ajuste <- entrena |> model(naive = NAIVE(value), snaive = SNAIVE(value))
ajuste |> forecast(h = 24) |> accuracy(serie) # 7. medir el error, no opinar
