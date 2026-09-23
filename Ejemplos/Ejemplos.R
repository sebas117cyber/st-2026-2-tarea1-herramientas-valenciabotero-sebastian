# paquetes permitidps

library(tidyverse);library(purrr);library(patchwork)

# datos ejemplo 

Prueba_leer <- leer_serie(nottem)
attr(Prueba_leer, "fuente")


Graficar_serie(Prueba_leer, titulo = "Temperatura en Nottingham")

correlograma(Prueba_leer)
