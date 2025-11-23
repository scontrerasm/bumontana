  library(tidyverse)
  library(ggplot2)
  library(vegan)    
  library(tibble)

dat <- read.csv(
  "Reflectance data.csv",
  check.names = TRUE,
  stringsAsFactors = FALSE
)
dat <- dat[, !is.na(names(dat)) & names(dat) != "", drop = FALSE]

dat <- dat %>%
  mutate(
    Poblacion = factor(trimws(Poblacion)),
    Sexo      = factor(trimws(Sexo))
  )

# Etiquetas de parches (ajusta si tus prefijos difieren)
patch_labels <- c(
  Nuca   = "Nuca",
  DorsoA = "Dorso_alto",
  DorsoM = "Dorso_medio",
  DorsoB = "Dorso_bajo",
  CobA   = "Cobertoras_alares"
)
patches <- names(patch_labels)

norm_pattern <- paste0("^(", paste(patches, collapse="|"), ")_(R|G|B|UV)_Norm$")

long_norm <- dat %>%
  pivot_longer(
    cols = matches(norm_pattern),
    names_to  = c("Patch","Channel"),
    values_to = "Value",
    names_pattern = "^([A-Za-z]+)_(R|G|B|UV)_Norm$"
  )

wide_by_patch <- long_norm %>%
  mutate(
    Patch   = factor(Patch, levels = patches, labels = unname(patch_labels)),
    Channel = factor(Channel, levels = c("R","G","B","UV"))
  ) %>%
  pivot_wider(
    id_cols    = c(Especimen, Poblacion, Sexo, Patch),
    names_from = Channel,
    values_from= Value,
    values_fill= NA_real_
  ) %>%
  drop_na(R,G,B)  # para análisis RGB

# ============================================================
# A) Varianza entre individuos por población (promediando parches)
#    -> mide variación visible (RGB) entre individuos dentro de cada población
# ============================================================
by_ind <- wide_by_patch %>%
  group_by(Especimen, Poblacion, Sexo) %>%
  summarise(across(c(R,G,B), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop")

# Tabla de varianzas por población y canal
var_RGB_by_pop <- by_ind %>%
  pivot_longer(c(R,G,B), names_to = "Canal", values_to = "Valor") %>%
  group_by(Poblacion, Canal) %>%
  summarise(Varianza = var(Valor, na.rm = TRUE),
            Media    = mean(Valor, na.rm = TRUE),
            n = dplyr::n(),
            .groups = "drop")

cat("\n[A] Varianza entre individuos (promedio por individuo):\n")
print(var_RGB_by_pop)

# Gráfico de barras de varianza por canal
ggplot(var_RGB_by_pop, aes(Poblacion, Varianza, fill = Poblacion)) +
  geom_col(width = 0.7, alpha = 0.85) +
  facet_wrap(~ Canal, nrow = 1) +
  labs(title = "Varianza RGB entre individuos por población",
       y = "Varianza (promedio por individuo)", x = "") +
  theme_minimal() + theme(legend.position = "none")

# ============================================================
# B) Dispersión multivariada en RGB POR PARCHE
#    -> ¿alguna población es más "dispersa" en el espacio RGB dentro de un parche?
# ============================================================
for (pt in levels(wide_by_patch$Patch)) {
  d_pt <- wide_by_patch %>%
    filter(Patch == pt) %>%
    drop_na(R,G,B)
  if (nrow(d_pt) < 5) next
  
  rgb_mat <- as.matrix(d_pt[,c("R","G","B")])
  disp    <- betadisper(dist(rgb_mat), d_pt$Poblacion)
  
  cat("\n[B] Dispersión RGB por parche -", pt, "\n")
  print(permutest(disp))  # prueba de homogeneidad de dispersiones
  
  # Plot rápido (base R) — opcional
  # plot(disp, main = paste0("Dispersión RGB - ", pt))
}

# ============================================================
# C) Variación dentro de individuo entre parches
#    -> si hay ruido de captura/iluminación, suele verse como mayor SD entre parches
# ============================================================
# Desv. estándar de R,G,B a través de parches para CADA individuo
sd_within_ind <- wide_by_patch %>%
  group_by(Especimen, Poblacion) %>%
  summarise(
    sd_R = sd(R, na.rm = TRUE),
    sd_G = sd(G, na.rm = TRUE),
    sd_B = sd(B, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(starts_with("sd_"),
               names_to = "Canal", values_to = "SD_parches") %>%
  mutate(Canal = recode(Canal, sd_R="R", sd_G="G", sd_B="B"))

cat("\n[C] SD dentro de individuo entre parches (por población y canal):\n")
print(sd_within_ind %>%
        group_by(Poblacion, Canal) %>%
        summarise(Media_SD = mean(SD_parches, na.rm = TRUE),
                  Mediana_SD = median(SD_parches, na.rm = TRUE),
                  n_ind = dplyr::n(),
                  .groups = "drop"))

# Violin/box: distribución de SD_intra por población
ggplot(sd_within_ind, aes(Poblacion, SD_parches, fill = Poblacion)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.2, outlier.shape = 21) +
  facet_wrap(~ Canal, nrow = 1) +
  labs(title = "Variación dentro de individuo entre parches (SD)",
       y = "SD a través de parches", x = "") +
  theme_minimal() + theme(legend.position = "none")

# ============================================================
# Notas de interpretación:
# - Si en [A] una población (p.ej. cucullata) muestra varianzas más altas en R/G/B,
#   hay mayor heterogeneidad entre sus individuos en RGB.
# - Si en [B] el permutest de betadisper da p pequeño en un parche,
#   esa población tiene mayor dispersión multivariada (posible ruido o variación real).
# - Si en [C] su SD dentro de individuo es alta, podría apuntar a inconsistencias
#   entre fotos/parches (iluminación) o a variación real entre regiones del plumaje.
# ============================================================
