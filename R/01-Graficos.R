# Grafico de serie



Graficar_serie <- function(serie, titulo) {
  
  colnames(serie) <- c("t", "fecha", "valores")
  # nivel general
  nivel_general <- mean(serie$valores, na.rm = TRUE)
  
  # Construir el gráfico
  p <- ggplot(serie, aes(x = fecha, y = valores)) +
    
    geom_line(colour = "#00B2EE", linewidth = 0.8) +
    
    # muestra del nivel
    geom_hline(
      aes(yintercept = nivel_general, linetype = "Nivel general"),
      colour    = "#8B2323",
      linewidth = 0.8
    ) +
    
    # Leyenda manual para la línea de referencia
    scale_linetype_manual(
      name   = "Referencia",
      values = c("Nivel general" = "dashed")
    ) +
    
    # Etiquetas del gráfico
    labs(
      title    = titulo,
      subtitle = NULL,
      x        = "Fecha",
      y        = "Valor"
    ) +
    
    # Tema limpio
    theme_minimal(base_size = 13)
  
  # Mostrar el gráfico
  print(p)
}




