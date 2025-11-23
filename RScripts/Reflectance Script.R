
library(tidyverse)
library(ggplot2)
library(vegan)      
library(emmeans)    
library(tibble)
library(grid)     

pollos <- read.csv("Reflectance data.csv", check.names = TRUE, stringsAsFactors = FALSE)
pollos <- pollos[, !is.na(names(pollos)) & names(pollos) != "", drop = FALSE]

pollos <- pollos %>%
  rename(Specimen   = Especimen, Population = Poblacion, Sex = Sexo) %>%
  mutate(Population = factor(trimws(Population)), Sex = factor(recode(trimws(Sex), "H" = "F", "M" = "M", .default = trimws(Sex))))


patch_codes  <- c("Nuca","DorsoA","DorsoM","DorsoB","CobA")
patch_labels <- c(Nuca = "Nape", DorsoA = "Upper back", DorsoM = "Mid back", DorsoB = "Lower back", CobA   = "Wing coverts")

# Patrón de las columnas; <Patch>_(R|G|B|UV)_Norm
norm_pattern <- paste0("^(", paste(patch_codes, collapse="|"), ")_(R|G|B|UV)_Norm$")
long_norm <- pollos %>%
  pivot_longer(cols = matches(norm_pattern), names_to  = c("Patch","Channel"),  values_to = "Value", names_pattern = "^([A-Za-z]+)_(R|G|B|UV)_Norm$")

wide_by_patch <- long_norm %>%
  mutate(Patch   = factor(Patch, levels = patch_codes, labels = unname(patch_labels)), Channel = factor(Channel, levels = c("R","G","B","UV"))) %>%
  pivot_wider(id_cols = c(Specimen, Population, Sex, Patch), names_from  = Channel, values_from = Value, values_fill = NA_real_) %>%
  mutate(VIS_sum   = R + G + B, UV_ratio  = UV / pmax(VIS_sum, 1e-6), UV_minusV = UV - (VIS_sum/3))

# Promedio por individuos
by_ind <- wide_by_patch %>%
  group_by(Specimen, Population, Sex) %>%
  summarise(across(c(R, G, B, UV, UV_ratio, UV_minusV), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  drop_na(R, G, B, UV)

stopifnot(nrow(by_ind) >= 5)

#PCA vis+UV
X_rgbu <- scale(select(by_ind, R, G, B, UV))
pca_rgbu <- prcomp(X_rgbu, center = TRUE, scale. = TRUE)
cat("\n[PCA RGB+UV] Explained variance:\n"); print(summary(pca_rgbu))

scores_rgbu  <- as.data.frame(pca_rgbu$x) %>%
  bind_cols(select(by_ind, Specimen, Population, Sex))
loads_rgbu <- as.data.frame(pca_rgbu$rotation[, 1:2]) %>%
  rownames_to_column("Channel")

ggplot(scores_rgbu, aes(PC1, PC2, color = Population)) +
  geom_point(aes(shape = Sex), size = 3, alpha = 0.9) +
  stat_ellipse(aes(group = Population), level = 0.68, linetype = 2) +
  geom_segment(data = loads_rgbu,
               aes(x = 0, y = 0, xend = PC1*3, yend = PC2*3),
               arrow = arrow(length = unit(0.25,"cm")), color = "black") +
  geom_text(data = loads_rgbu,
            aes(x = PC1*3.2, y = PC2*3.2, label = Channel),
            color = "black") +
  labs(title = "Global PCA by individual",
       subtitle = "Channels: R, G, B, UV (normalized)") +
  theme_minimal()

#PCA Sólo visible 

X_rgb <- scale(select(by_ind, R, G, B))
pca_rgb <- prcomp(X_rgb, center = TRUE, scale. = TRUE)
cat("\n[PCA RGB] Explained variance:\n"); print(summary(pca_rgb))

scores_rgb <- as.data.frame(pca_rgb$x) %>%
  bind_cols(select(by_ind, Specimen, Population, Sex))
loads_rgb <- as.data.frame(pca_rgb$rotation[, 1:2]) %>%
  rownames_to_column("Channel")

ggplot(scores_rgb, aes(PC1, PC2, color = Population)) +
  geom_point(aes(shape = Sex), size = 3, alpha = 0.9) +
  stat_ellipse(aes(group = Population), level = 0.68, linetype = 2) +
  geom_segment(data = loads_rgb,
               aes(x = 0, y = 0, xend = PC1*3, yend = PC2*3),
               arrow = arrow(length = unit(0.25,"cm")), color = "black") +
  geom_text(data = loads_rgb,
            aes(x = PC1*3.2, y = PC2*3.2, label = Channel),
            color = "black") +
  labs(title = "PCA, RGB only (averaged by individual)",
       subtitle = "Channels: R, G, B (normalized)") +
  theme_minimal()

#PermANOVA y Betadisper 

rgb_mat <- as.matrix(select(by_ind, R, G, B))
dist_rgb_ind <- dist(rgb_mat, method = "euclidean")
perm_rgb <- adonis2(dist_rgb_ind ~ Population, data = by_ind, permutations = 999)
cat("\n[PERMANOVA RGB]\n"); print(perm_rgb)
disp_rgb <- betadisper(dist_rgb_ind, by_ind$Population)
cat("\n[betadisper RGB]\n"); print(permutest(disp_rgb))

uv_vec <- as.matrix(select(by_ind, UV))
dist_uv_ind <- dist(uv_vec, method = "euclidean")
perm_uv <- adonis2(dist_uv_ind ~ Population, data = by_ind, permutations = 999)
cat("\n[PERMANOVA UV]\n"); print(perm_uv)
disp_uv <- betadisper(dist_uv_ind, by_ind$Population)
cat("\n[betadisper UV]\n"); print(permutest(disp_uv))

#NMDSs 
#UV
set.seed(.Random.seed)
nmds_uv <- metaMDS(uv_vec, distance = "euclidean", k = 2, trymax = 100)
cat("\n[NMDS UV] Stress:\n"); print(nmds_uv$stress)
scores_nmds_uv <- as.data.frame(scores(nmds_uv)) %>%
  bind_cols(select(by_ind, Population, Sex))
ggplot(scores_nmds_uv, aes(NMDS1, NMDS2, color = Population)) +
  geom_point(aes(shape = Sex), size = 3, alpha = 0.9) +
  stat_ellipse(aes(group = Population), level = 0.68, linetype = 2) +
  labs(title = "NMDS – UV only", subtitle = "Euclidean distance on UV") +
  theme_minimal()
#Vis
set.seed(.Random.seed)
nmds_rgb <- metaMDS(rgb_mat, distance = "euclidean", k = 2, trymax = 100)
cat("\n[NMDS RGB] Stress:\n"); print(nmds_rgb$stress)

sc_rgb_sites <- scores(nmds_rgb, display = "sites") 
sc_rgb_df    <- as.data.frame(sc_rgb_sites)

by_ind_labels <- by_ind %>%
  dplyr::ungroup() %>%
  dplyr::select(Specimen, Population, Sex)

stopifnot(nrow(sc_rgb_df) == nrow(by_ind_labels))
scores_nmds_rgb <- cbind(sc_rgb_df, by_ind_labels)

# Graficar
ggplot(scores_nmds_rgb, aes(NMDS1, NMDS2, color = Population)) +
  geom_point(aes(shape = Sex), size = 3, alpha = 0.9) +
  stat_ellipse(aes(group = Population), level = 0.68, linetype = 2) +
  labs(title = "NMDS – RGB", subtitle = "Euclidean distance on R,G,B") +
  theme_minimal()

#ANOVA solo UV y Tukey de reflectancias UV 
a_uv <- aov(UV ~ Population, data = by_ind)
cat("\n[ANOVA UV]\n"); print(summary(a_uv))
cat("\n[EMMs UV] (Tukey)\n"); print(emmeans(a_uv, pairwise ~ Population))

ggplot(by_ind, aes(Population, UV, fill = Population)) +
  geom_boxplot(alpha = 0.75, outlier.shape = 21) +
  labs(title = "UV reflectance by population (mean per individual)", y = "UV (normalized)", x = "Population") +
  theme_minimal() + theme(legend.position = "none")

#Distancias entre poblaciones en espacio visible y UV 
centroids <- by_ind %>%
  group_by(Population) %>%
  summarise(R = mean(R), G = mean(G), B = mean(B), UV = mean(UV), .groups="drop")

D_rgb_pop <- dist(centroids[, c("R","G","B")], method = "euclidean")
D_uv_pop  <- dist(centroids[, "UV", drop = FALSE], method = "euclidean")
cat("\n[Population distances RGB]:\n"); print(as.matrix(D_rgb_pop))
cat("\n[Population distances UV]:\n");  print(as.matrix(D_uv_pop))

hc_rgb <- hclust(D_rgb_pop, method = "average")
hc_uv  <- hclust(D_uv_pop,  method = "average")
op <- par(mfrow = c(1,2))
plot(hc_rgb, main = "Population distances (RGB)", xlab = "", sub = "", labels = centroids$Population)
plot(hc_uv,  main = "Population distances (UV)",  xlab = "", sub = "", labels = centroids$Population)
par(op)

#Diagnósticos de las varianzas

#Entre individuos
by_ind_rgb <- by_ind %>% select(Specimen, Population, Sex, R, G, B)
var_RGB_by_pop <- by_ind_rgb %>%
  pivot_longer(c(R,G,B), names_to = "Channel", values_to = "Value") %>%
  group_by(Population, Channel) %>%
  summarise(Variance = var(Value, na.rm = TRUE), Mean = mean(Value, na.rm = TRUE), n = dplyr::n(), .groups = "drop")
cat("\n[Variance among individuals RGB]:\n"); print(var_RGB_by_pop)
var_RGB_by_pop <- var_RGB_by_pop %>%
  mutate(Channel = factor(Channel, levels = c("R","G","B")))


ggplot(var_RGB_by_pop, aes(Population, Variance, fill = Population)) +
  geom_col(width = 0.7, alpha = 0.85) +
  facet_wrap(~ Channel, nrow = 1) +
  labs(title = "RGB variance among individuals by population",
       y = "Variance (per individual mean)", x = "") +
  theme_minimal() + theme(legend.position = "none")

#Dispersión por Parche
for (pt in levels(wide_by_patch$Patch)) 
{d_pt <- wide_by_patch %>%
  filter(Patch == pt) %>%
  drop_na(R,G,B)
if (nrow(d_pt) < 5) next
m_pt <- as.matrix(d_pt[, c("R","G","B")])
disp_pt <- betadisper(dist(m_pt), d_pt$Population)
cat("\n[betadisper RGB] Patch:", pt, "\n")
print(permutest(disp_pt))}

## Desviación intra-individual entre parches
sd_within_ind <- wide_by_patch %>%
  group_by(Specimen, Population) %>%
  summarise(sd_R = sd(R, na.rm = TRUE), sd_G = sd(G, na.rm = TRUE), sd_B = sd(B, na.rm = TRUE), .groups = "drop") %>%
  pivot_longer(starts_with("sd_"), names_to = "Channel", values_to = "SD_patches") %>% 
  mutate(Channel = recode(Channel, sd_R="R", sd_G="G", sd_B="B"))

cat("\n[SD within individual across patches]:\n")
print(sd_within_ind %>%
        group_by(Population, Channel) %>%
        summarise(Mean_SD = mean(SD_patches, na.rm = TRUE), Median_SD = median(SD_patches, na.rm = TRUE), n_ind = dplyr::n(), .groups = "drop"))
sd_within_ind <- sd_within_ind %>%
  mutate(Channel = factor(Channel, levels = c("R","G","B")))
ggplot(sd_within_ind, aes(Population, SD_patches, fill = Population)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.2, outlier.shape = 21) +
  facet_wrap(~ Channel, nrow = 1) +
  labs(title = "Within-individual variation across patches (SD)",
       y = "SD across patches", x = "") +
  theme_minimal() + theme(legend.position = "none")
#Como heatmap 
D_rgb_pop <- dist(centroids[, c("R","G","B")], method = "euclidean")
D_uv_pop  <- dist(centroids[, "UV", drop = FALSE], method = "euclidean")

# Convertir a matrices
mat_rgb <- as.matrix(D_rgb_pop)
mat_uv  <- as.matrix(D_uv_pop)

# Derretir para ggplot
df_rgb <- melt(mat_rgb, varnames = c("Pop1","Pop2"), value.name = "Distance")
df_uv  <- melt(mat_uv,  varnames = c("Pop1","Pop2"), value.name = "Distance")

# Revisar
head(df_rgb)

# Añadir etiquetas para clarificar (opcional)
names(df_rgb) <- c("Pop1","Pop2","Distance")
names(df_uv)  <- c("Pop1","Pop2","Distance")

# Heatmap RGB
ggplot(df_rgb, aes(Pop1, Pop2, fill = Distance)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "green", high = "red") +
  labs(title = "Heatmap – Distancias RGB",
       x = "Population", y = "Population") +
  theme_minimal()

# Heatmap UV
ggplot(df_uv, aes(Pop1, Pop2, fill = Distance)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Heatmap – Distancias UV",
       x = "Population", y = "Population") +
  theme_minimal()
