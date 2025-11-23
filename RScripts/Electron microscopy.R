timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
logfile <- paste0("SEM_resultados_", timestamp, ".txt")

zz <- file(logfile, open = "wt")
sink(zz, type = "output", split = TRUE)
sink(zz, type = "message", append = TRUE)
cat("Salida y mensajes registrados en:", logfile, "\n\n")
pkgs <- c("tidyverse","lme4","lmerTest","broom","ggplot2",
          "performance","car","rstatix","FactoMineR","factoextra",
          "vegan","rlang","emmeans","ggrepel")
invisible(lapply(pkgs, function(p){
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}))

# No entrar a modo debug si hay errores
options(error = NULL)

# ------------------------------------------------
# 1) Cargar datos (CSV ancho) y normalizar campos
# ------------------------------------------------
# Requiere columnas: Specimen, Population, Sex, Patch, y variables morfométricas
datos <- read.csv("Electron microscopy data.csv", check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(all(c("Specimen","Population","Sex","Patch") %in% names(datos)))

datos <- datos %>%
  mutate(
    Specimen   = factor(Specimen),
    Population = factor(Population),
    Sex        = factor(Sex),
    Patch      = factor(Patch)
  )

# Variables morfométricas (ajusta si usas otros nombres)
vars <- c("Barb_width","Barbule_length","Barbule_width","Barbule_angle","Barbule_spacing")
vars <- intersect(vars, names(datos))
if (length(vars) == 0) stop("No se encontraron columnas de variables morfométricas.")

# Convertir a numéricas con cuidado
datos <- datos %>% mutate(across(all_of(vars), ~ suppressWarnings(as.numeric(.))))

# Derivado opcional: densidad lineal de barbillas
if ("Barbule_spacing" %in% vars) {
  datos <- datos %>% mutate(Barbule_density = 1/Barbule_spacing)
  vars  <- unique(c(vars, "Barbule_density"))
}

# --------------------------------------------------------------------
# 2) Función robusta para ajustar modelos por variable dentro de un Patch
# --------------------------------------------------------------------
fit_var_models <- function(subdat, v){
  sub <- droplevels(
    subdat %>%
      dplyr::select(Specimen, Population, Sex, all_of(v)) %>%
      dplyr::rename(y = !!rlang::sym(v))
  ) %>% dplyr::filter(!is.na(y))
  
  if (NROW(sub) < 5) {
    cat("\n[Salto]", v, ": <5 filas no-NA en este parche.\n")
    return(invisible(NULL))
  }
  
  fixed_terms <- c()
  if (nlevels(sub$Population) >= 2) fixed_terms <- c(fixed_terms, "Population")
  if (nlevels(sub$Sex)        >= 2) fixed_terms <- c(fixed_terms, "Sex")
  fixed_part <- if (length(fixed_terms) > 0) paste(fixed_terms, collapse = " + ") else "1"
  
  rand_part <- if (dplyr::n_distinct(sub$Specimen) >= 2) "(1|Specimen)" else NULL
  rhs <- paste(c(fixed_part, rand_part), collapse = " + ")
  fml <- stats::as.formula(paste("y ~", rhs))
  
  cat("\nVariable:", v, "\nFórmula usada:", deparse(fml), "\n")
  
  if (!is.null(rand_part)) {
    out <- try(
      lme4::lmer(fml, data = sub, REML = TRUE, na.action = na.exclude,
                 control = lme4::lmerControl(optimizer = "bobyqa",
                                             calc.derivs = FALSE,
                                             optCtrl = list(maxfun = 1e5))),
      silent = TRUE
    )
    if (inherits(out, "try-error")) {
      cat("[Aviso] lmer falló; usando LM simple.\n")
      out <- stats::lm(stats::as.formula(paste("y ~", fixed_part)), data = sub)
    }
  } else {
    out <- stats::lm(stats::as.formula(paste("y ~", fixed_part)), data = sub)
  }
  
  if (inherits(out, "lmerMod")) {
    print(anova(out))
    print(summary(out))
  } else {
    print(anova(out))
    print(summary(out))
  }
  invisible(out)
}

# ----------------------------------------------------------
# 3) Bucle principal: analizar CADA PATCH por separado
# ----------------------------------------------------------
for (p in levels(datos$Patch)) {
  
  cat("\n========================\nPARCHE:", p, "\n========================\n")
  subdat <- datos %>% dplyr::filter(Patch == p) %>% droplevels()
  if (nrow(subdat) == 0) { cat("Sin filas para", p, "\n"); next }
  
  # --- QC: violines + jitter facetados ---
  sub_long <- subdat %>%
    dplyr::select(Specimen, Population, Sex, dplyr::all_of(vars)) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(vars), names_to = "Variable", values_to = "Valor")
  
  g_qc <- ggplot(sub_long, aes(x = Population, y = Valor, fill = Population)) +
    geom_violin(trim = FALSE, alpha = 0.3) +
    geom_jitter(width = 0.1, alpha = 0.6, size = 1) +
    facet_wrap(~Variable, scales = "free_y") +
    labs(title = paste("Per - population distribution -", p), y = "Valor", x = "Población") +
    theme_minimal()
  print(g_qc)
  
  # --- Resumen por individuo (para inspección) ---
  by_ind <- subdat %>%
    dplyr::group_by(Specimen, Population, Sex) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(vars),
                    list(mean = ~mean(.x, na.rm = TRUE),
                         sd   = ~sd(.x,   na.rm = TRUE),
                         n    = ~sum(!is.na(.x))),
                    .names = "{.col}_{.fn}"),
      .groups = "drop"
    )
  cat("\nResumen por individuo (primeras filas):\n")
  print(utils::head(by_ind, 10))
  
  # --- Modelos por variable (robustos) ---
  cat("\n--- MODELOS (robustos) ---\n")
  invisible(lapply(vars, function(v) fit_var_models(subdat, v)))
  
  # --- ANOVA / KRUSKAL por variable (robusto) ---
  cat("\n--- ANOVA / KRUSKAL por variable ---\n")
  for (v in vars){
    d <- subdat %>%
      dplyr::select(Population, Sex, Specimen, !!rlang::sym(v)) %>%
      dplyr::rename(y = !!rlang::sym(v)) %>%
      tidyr::drop_na(y) %>%
      droplevels()
    
    if (nrow(d) < 5) next
    
    fixed_terms <- c()
    if (nlevels(d$Population) >= 2) fixed_terms <- c(fixed_terms, "Population")
    if (nlevels(d$Sex)        >= 2) fixed_terms <- c(fixed_terms, "Sex")
    if (length(fixed_terms) == 0){
      cat("\nVariable:", v, "— sin factores con ≥2 niveles en este parche; salto.\n")
      next
    }
    
    fixed_part <- paste(fixed_terms, collapse = " + ")
    fml <- stats::as.formula(paste("y ~", fixed_part))
    cat("\nVariable:", v, " | Fórmula:", deparse(fml), "\n")
    
    if ("Population" %in% fixed_terms) {
      cat("Shapiro por Population:\n")
      print(d %>% dplyr::group_by(Population) %>%
              dplyr::summarise(p_shapiro = tryCatch(shapiro.test(y)$p.value, error = \(e) NA_real_), .groups="drop"))
      cat("Levene por Population:\n")
      print(car::leveneTest(y ~ Population, data = d))
    }
    if ("Sex" %in% fixed_terms) {
      cat("Shapiro por Sex:\n")
      print(d %>% dplyr::group_by(Sex) %>%
              dplyr::summarise(p_shapiro = tryCatch(shapiro.test(y)$p.value, error = \(e) NA_real_), .groups="drop"))
      cat("Levene por Sex:\n")
      print(car::leveneTest(y ~ Sex, data = d))
    }
    
    aov_fit <- aov(fml, data = d)
    print(summary(aov_fit))
    
    if ("Population" %in% fixed_terms) {
      cat("Kruskal-Wallis (Population):\n")
      print(rstatix::kruskal_test(y ~ Population, data = d))
    }
  }
  
  # --- CONTRASTES PAREADOS (Population) ---
  cat("\n--- CONTRASTES PAREADOS (Population) ---\n")
  for (v in vars) {
    d <- subdat %>%
      dplyr::select(Population, Sex, Specimen, !!rlang::sym(v)) %>%
      dplyr::rename(y = !!rlang::sym(v)) %>%
      tidyr::drop_na(y) %>%
      droplevels()
    
    if (nrow(d) < 5 || nlevels(d$Population) < 2) next
    
    grp_n <- table(d$Population)
    n_groups_ge2 <- sum(grp_n >= 2)
    
    lev <- tryCatch(car::leveneTest(y ~ Population, data = d), error = function(e) NULL)
    p_lev <- if (!is.null(lev)) lev[1, "Pr(>F)"] else NA_real_
    
    cat("\nVariable:", v, "| tamaños por grupo:", paste(names(grp_n), grp_n, sep="=", collapse=", "), "\n")
    
    # Tukey si homocedasticidad y ≥2 grupos con n≥2
    if (!is.na(p_lev) && p_lev > 0.05 && n_groups_ge2 >= 2) {
      fixed_terms <- c("Population")
      if (nlevels(d$Sex) >= 2) fixed_terms <- c(fixed_terms, "Sex")
      fml_tuk <- stats::as.formula(paste("y ~", paste(fixed_terms, collapse = " + ")))
      
      aov_fit <- aov(fml_tuk, data = d)
      cat("Contrastes Tukey (emmeans) — Levene p =", round(p_lev, 3),
          "| fórmula:", deparse(fml_tuk), "\n")
      tuk <- emmeans::emmeans(aov_fit, specs = ~ Population)
      print(emmeans::contrast(tuk, method = "pairwise", adjust = "tukey"))
      next
    }
    
    # Games–Howell si todos los grupos tienen n≥2
    if (all(grp_n >= 2)) {
      cat("Contrastes Games–Howell — (varianzas desiguales o desequilibrio)\n")
      gh <- rstatix::games_howell_test(d, y ~ Population) %>% dplyr::arrange(group1, group2)
      print(gh)
    } else {
      cat("No se puede aplicar Tukey ni Games–Howell (hay grupos con n=1).\n")
    }
    
    # Dunn (BH) tras Kruskal como plan C
    kw <- rstatix::kruskal_test(d, y ~ Population)
    if (!is.na(kw$p) && kw$p < 0.1 && nlevels(d$Population) >= 2) {
      cat("Post-hoc Dunn (BH) tras Kruskal (p =", round(kw$p, 4), ")\n")
      dunn <- rstatix::dunn_test(d, y ~ Population, p.adjust.method = "BH") %>%
        dplyr::arrange(group1, group2)
      print(dunn)
    } else {
      cat("Kruskal sin señal clara o datos insuficientes para Dunn; se omite.\n")
    }
  }
  
  # --- PCA por individuo (medias) → BIPLOT: puntos + hulls + vectores ---
  cat("\n--- PCA por individuo (medias) ---\n")
  for_pca <- subdat %>%
    dplyr::group_by(Specimen, Population, Sex) %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(vars), ~mean(.x, na.rm = TRUE)), .groups = "drop")
  
  num_mat <- for_pca %>%
    dplyr::select(dplyr::all_of(vars)) %>%
    dplyr::select(where(~sum(!is.na(.)) > 1)) %>%
    dplyr::select(where(~sd(.x, na.rm = TRUE) > 1e-9))
  
  if (ncol(num_mat) >= 2 && nrow(num_mat) >= 3) {
    pca_res <- FactoMineR::PCA(num_mat, graph = FALSE, scale.unit = TRUE)
    
    # Coordenadas de individuos
    coords <- as.data.frame(pca_res$ind$coord) %>%
      dplyr::mutate(Specimen = for_pca$Specimen,
                    Population = droplevels(for_pca$Population),
                    Sex = for_pca$Sex)
    
    pc1_var <- round(100 * pca_res$eig[1, "percentage of variance"], 1)
    pc2_var <- round(100 * pca_res$eig[2, "percentage of variance"], 1)
    
    # Polígonos por población
    hulls <- coords %>%
      dplyr::group_by(Population) %>%
      dplyr::filter(dplyr::n() >= 3) %>%
      dplyr::slice(chull(Dim.1, Dim.2)) %>%
      dplyr::ungroup()
    
    # Cargas de variables
    var_coord <- as.data.frame(pca_res$var$coord)
    if (nrow(var_coord) >= 1) {
      var_coord$Variable <- rownames(var_coord)
      xr <- range(coords$Dim.1, na.rm = TRUE); yr <- range(coords$Dim.2, na.rm = TRUE)
      vx <- max(abs(var_coord$Dim.1), na.rm = TRUE); vy <- max(abs(var_coord$Dim.2), na.rm = TRUE)
      s <- 0.8 * min( (diff(xr) / (2*vx)), (diff(yr) / (2*vy)) )
      if (!is.finite(s) || s <= 0) s <- 1
      var_plot <- var_coord %>% dplyr::mutate(PC1 = Dim.1 * s, PC2 = Dim.2 * s)
    } else {
      var_plot <- NULL
    }
    
    # Biplot
    g_biplot <- ggplot() +
      geom_polygon(data = hulls,
                   aes(x = Dim.1, y = Dim.2, group = Population, fill = Population),
                   alpha = 0.15, linewidth = 0.4) +
      geom_point(data = coords,
                 aes(x = Dim.1, y = Dim.2, color = Population),
                 size = 2, alpha = 0.9) +
      { if (!is.null(var_plot))
        geom_segment(data = var_plot,
                     aes(x = 0, y = 0, xend = PC1, yend = PC2),
                     arrow = arrow(length = unit(0.02, "npc")), linewidth = 0.5)
        else NULL } +
      { if (!is.null(var_plot))
        ggrepel::geom_text_repel(data = var_plot,
                                 aes(x = PC1, y = PC2, label = Variable),
                                 size = 3, max.overlaps = Inf)
        else NULL } +
      labs(title = paste("PCA biplot -", p),
           x = paste0("PC1 (", pc1_var, "%)"),
           y = paste0("PC2 (", pc2_var, "%)")) +
      theme_minimal()
    print(g_biplot)
    
  } else {
    cat("No hay suficientes variables con varianza para PCA en", p, "\n")
  }
  
  # --- PERMANOVA (medias por espécimen, z-score, euclidiana) + betadisper ---
  cat("\n--- PERMANOVA + betadisper (multivariante) ---\n")
  X  <- for_pca %>% dplyr::select(dplyr::all_of(vars))
  Xz <- X %>% dplyr::mutate(dplyr::across(everything(), ~ as.numeric(scale(.))))
  keep_cols <- which(colSums(!is.na(Xz)) > 0 & apply(Xz, 2, function(z) sd(z, na.rm=TRUE) > 1e-9))
  Xz <- Xz[, keep_cols, drop = FALSE]
  keep_rows <- which(rowSums(!is.na(Xz)) > 0)
  Xz <- Xz[keep_rows, , drop = FALSE]
  meta <- for_pca[keep_rows, c("Specimen","Population","Sex"), drop = FALSE] %>% droplevels()
  
  ok_perma <- nrow(Xz) >= 6 && ncol(Xz) >= 2 && nlevels(meta$Population) >= 2
  if (!ok_perma) {
    cat("Datos insuficientes para PERMANOVA en", p, "\n")
  } else {
    set.seed(123)
    per <- vegan::adonis2(Xz ~ Population + Sex, data = meta,
                          permutations = 999, method = "euclidean", by = "margin")
    cat("\nResultado PERMANOVA (adonis2) —", p, "\n")
    print(per)
    
    grp_n <- table(meta$Population)
    if (sum(grp_n >= 2) >= 2) {
      D  <- stats::dist(Xz, method = "euclidean")
      bd <- vegan::betadisper(D, group = meta$Population, type = "centroid")
      cat("\nbetadisper (dispersiones) por Population —", p, "\n")
      suppressWarnings({
        print(anova(bd))
        print(vegan::permutest(bd, permutations = 999))
      })
    } else {
      cat("betadisper omitido: se requieren ≥2 grupos con n≥2 (tamaños:",
          paste(names(grp_n), grp_n, sep="=", collapse=", "), ")\n")
    }
  }
  
  # --- NMDS (opcional) con POLÍGONOS ---
  cat("\n--- NMDS (opcional) ---\n")
  num_mat_nmds <- for_pca %>% dplyr::select(dplyr::all_of(vars))
  if (ncol(num_mat_nmds) >= 1 && nrow(num_mat_nmds) >= 3) {
    num_mat_nmds <- num_mat_nmds %>% dplyr::mutate(across(everything(), ~scale(.)[,1]))
    keep_cols <- which(colSums(!is.na(num_mat_nmds)) > 0 &
                         apply(num_mat_nmds, 2, function(z) sd(z, na.rm=TRUE) > 1e-9))
    num_mat_nmds <- num_mat_nmds[, keep_cols, drop = FALSE]
    num_mat_nmds <- num_mat_nmds[rowSums(!is.na(num_mat_nmds)) > 0, , drop = FALSE]
    
    if (ncol(num_mat_nmds) >= 2 && nrow(num_mat_nmds) >= 4 &&
        dplyr::n_distinct(for_pca$Population) >= 2) {
      
      set.seed(123)
      nmds <- try(vegan::metaMDS(num_mat_nmds, distance = "euclidean", k = 2,
                                 trymax = 200, autotransform = FALSE, trace = FALSE),
                  silent = TRUE)
      
      if (!inherits(nmds, "try-error")) {
        nmds_scores <- as.data.frame(vegan::scores(nmds, display = "sites")) %>%
          dplyr::mutate(Specimen = for_pca$Specimen,
                        Population = droplevels(for_pca$Population),
                        Sex = for_pca$Sex)
        
        hulls_nmds <- nmds_scores %>%
          dplyr::group_by(Population) %>%
          dplyr::filter(dplyr::n() >= 3) %>%
          dplyr::slice(chull(NMDS1, NMDS2)) %>%
          dplyr::ungroup()
        
        g_nmds <- ggplot(nmds_scores, aes(NMDS1, NMDS2, color = Population, fill = Population)) +
          geom_polygon(data = hulls_nmds, aes(group = Population), alpha = 0.15, linewidth = 0.4) +
          geom_point(size = 2, alpha = 0.9) +
          labs(title = paste("NMDS -", p),
               x = "NMDS1", y = "NMDS2",
               subtitle = paste("Stress =", round(nmds$stress, 3))) +
          theme_minimal()
        print(g_nmds)
      } else {
        cat("NMDS no se pudo ajustar de forma estable en", p, "\n")
      }
    } else {
      cat("No hay suficientes datos/varianza para NMDS en", p, "\n")
    }
  } else {
    cat("No hay suficientes datos para NMDS en", p, "\n")
  }
  
  # --- Boxplots consolidados (1 figura por parche con todos los rasgos) ---
  sub_long_bp <- subdat %>%
    dplyr::select(Population, dplyr::all_of(vars)) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(vars),
                        names_to = "Variable", values_to = "Valor") %>%
    dplyr::filter(!is.na(Valor)) %>%
    droplevels()
  
  g_bp <- ggplot(sub_long_bp, aes(x = Population, y = Valor, fill = Population)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    geom_jitter(width = 0.12, alpha = 0.6, size = 0.9) +
    facet_wrap(~ Variable, scales = "free_y", ncol = 3) +
    labs(title = paste(p,"Boxplots"),
         x = "Population", y = NULL) +
    theme_minimal() +
    theme(legend.position = "bottom",
          strip.text = element_text(face = "bold"))
  print(g_bp)
}

sink(type = "message")
sink(type = "output")
close(zz)
cat("Archivo de salida cerrado correctamente.\n")