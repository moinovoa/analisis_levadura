# Librerías ----
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(plotly)
library(nnet)
library(patchwork)

# Creación base inicial ----
ruta <-  "/Users/aloolrey/Desktop/Simulación Est/Proyecto_1/yeast.xlsx"
datos <- read_excel(ruta)
datos <- datos %>%
  mutate(across(where(is.character), ~ na_if(.x, "NA"))) %>%
  mutate(across(c(mcg, gvh, alm, mit, erl, pox, vac, nuc), ~as.numeric(.x))) %>%
  mutate(name = as.factor(name))

nas <- datos %>%
  filter(if_any(everything(), ~is.na(.x)))
colSums(is.na(nas))

# Gráficos iniciales
p1 <- ggplot(datos, aes(x = alm)) + 
  geom_histogram(fill = "steelblue", color = "white") + 
  theme_minimal()

p2 <- ggplot(datos, aes(x = mit)) + 
  geom_histogram(fill = "steelblue", color = "white") + 
  theme_minimal()

p3 <- ggplot(datos, aes(x = erl)) + 
  geom_histogram(fill = "steelblue", color = "white") + 
  theme_minimal()

p4 <- ggplot(datos, aes(x = nuc)) + 
  geom_histogram(fill = "steelblue", color = "white") + 
  theme_minimal()

(p1 + p2) / (p3 + p4)

# Imputación
datos_ajustados <- datos %>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), median(., na.rm = TRUE), .)))

# Estadísticas ----
summary(datos_ajustados)

# Regresión logística multinomial ----
modelo <- multinom(name ~ mcg + gvh + alm + mit + erl + pox + vac + nuc, data = datos_ajustados)
summary(modelo)

## N muestras  ----
set.seed(12345)

coefs_list <- list()
N = 100

for(i in 1:N) {
  muestra <- datos_ajustados %>% sample_n(nrow(.), replace = TRUE)
  
  m <- multinom(name ~ mcg + gvh + alm + mit + erl + pox + vac + nuc, 
                data = muestra, trace = FALSE)
  
  # extraer coeficientes
  coefs <- summary(m)$coefficients
  
  # guardar SOLO variables de interés
  coefs_list[[i]] <- coefs[, c("mcg", "gvh", "erl", "pox")]
}

## Densidades ----
coefs_ds <- bind_rows(lapply(coefs_list, function(x) {
  df <- as.data.frame(x)
  df$categoria <- rownames(x)
  return(df)
}), .id = "iteracion")

coefs_long <- coefs_ds %>%
  pivot_longer(cols = c(mcg, gvh, erl, pox),
               names_to = "variable",
               values_to = "valor")

ggplot(coefs_long, aes(x = valor)) +
  geom_density(fill = "lightblue", alpha = 0.4) +
  facet_grid(categoria ~ variable, scales = "free") +
  theme_minimal()

### Intervalos de confianza ----
intervalos <- coefs_long %>%
  group_by(categoria, variable) %>%
  summarise(
    media = mean(valor),
    inf = quantile(valor, 0.025),
    sup = quantile(valor, 0.975),
    .groups = 'drop'
  ) %>%
  ungroup()

ggplot(coefs_long, aes(x = valor)) +
  geom_density(fill = "lightblue", alpha = 0.4) + 
  facet_grid(categoria ~ variable, scales = "free") +
  
  geom_vline(data = intervalos, aes(xintercept = inf), 
             color = "red", linetype = "dashed") +
  geom_vline(data = intervalos, aes(xintercept = sup), 
             color = "red", linetype = "dashed") +
  
  geom_vline(data = intervalos, aes(xintercept = media), 
             color = "darkblue") +
  
  theme_minimal() +
  labs(title = "Densidades de los Parámetros con IC 95%",
       subtitle = "Líneas punteadas indican el intervalo de confianza del 95%",
       x = "Valor de la Beta",
       y = "Densidad")

# Monte Carlo (Validación Train/Test) ----
M = 1000
resultados <- numeric(M)
n = nrow(datos_ajustados)

for(i in 1:M) {
  idx_train <- sample(1:n, size = 0.8 * n)
  
  train <- datos_ajustados[idx_train, ]
  test  <- datos_ajustados[-idx_train, ]
  
  modelo_mc <- multinom(name ~ mcg + gvh + alm + mit + erl + pox + vac + nuc,
                        data = train, trace = FALSE)
  
  pred <- predict(modelo_mc, newdata = test)
  
  acc <- mean(pred == test$name)
  
  resultados[i] <- acc
}

# Resultado final Monte Carlo
media_accuracy <- mean(resultados)
sd_accuracy <- sd(resultados)

cat("\n--- Resultados Monte Carlo (Train/Test) ---\n")
cat("Accuracy promedio:", round(media_accuracy * 100, 2), "%\n")
cat("Desviación estándar:", round(sd_accuracy * 100, 2), "%\n")


# Datos perturbados (Análisis de Sensibilidad) ----
sigmas <- seq(0.005, 0.035, by = 0.01)
aciertos_robustez <- numeric(length(sigmas))

cat("\n--- Resultados de Análisis de Sensibilidad ---\n")
for(i in 1:length(sigmas)) {
  
  sigma_ruido <- sigmas[i]
  
  datos_perturbados <- datos_ajustados %>% 
    mutate(across(c(mcg, gvh, alm, mit, erl, pox, vac, nuc), 
                  ~ .x + rnorm(n(), mean = 0, sd = sigma_ruido)))
  
  modelo_perturbado <- multinom(name ~ mcg + gvh + alm + mit + erl + pox + vac + nuc, 
                                data = datos_perturbados, trace = FALSE)
  
  predicciones <- predict(modelo_perturbado, newdata = datos_ajustados)
  
  aciertos_robustez[i] <- mean(predicciones == datos_ajustados$name)
  
  cat("Sigma =", sigma_ruido, 
      "- Acierto =", round(aciertos_robustez[i] * 100, 2), "%\n")
}

resultados_robustez <- data.frame(
  sigma = sigmas,
  acierto = aciertos_robustez,
  porcentaje = aciertos_robustez * 100
)

# Gráfica de Análisis de Sensibilidad
ggplot(resultados_robustez, aes(x = sigma, y = porcentaje)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "darkred", size = 3) +
  labs(
    title = "Robustez del modelo ante distintos niveles de ruido",
    subtitle = "Análisis de sensibilidad variando la desviación estándar del error (sigma)",
    x = "Nivel de ruido (sigma)",
    y = "Porcentaje de acierto (%)"
  ) +
  theme_minimal()