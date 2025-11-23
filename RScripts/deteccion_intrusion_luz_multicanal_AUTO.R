# deteccion_intrusion_luz_multicanal_AUTO.R
# Autor: ChatGPT
# Descripción: Versión multicanal (RGB + UV) con mapeo manual para tus columnas:
#   - Negro_R/G/B/UV
#   - Blanco_R/G/B/UV
# ID de foto: "Especimen"
#
# Ejecuta con: source("deteccion_intrusion_luz_multicanal_AUTO.R")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(MASS)        # cov.rob
  library(stringr)
  library(purrr)
  library(tidyselect)
})

# ===================== PARÁMETROS =====================
input_path <- "Reflectance data.csv"
channels <- c("R","G","B","UV")

# Mapeo manual basado en tu CSV
manual_map <- list(
  R  = list(white="Blanco_R", black="Negro_R"),
  G  = list(white="Blanco_G", black="Negro_G"),
  B  = list(white="Blanco_B", black="Negro_B"),
  UV = list(white="Blanco_UV", black="Negro_UV")
)

# Columna ID detectada
id_col <- "Especimen"

# Umbrales (ajustables)
alpha_mahal <- 0.01
k_mad <- 3
min_contrast_rel <- 0.10
auto_normalize_255 <- TRUE
# =====================================================

mad_flag <- function(x, k = 3) {
  m <- stats::median(x, na.rm = TRUE)
  md <- stats::mad(x, constant = 1.4826, na.rm = TRUE)
  if (is.na(md) || md == 0) {
    rep(FALSE, length(x))
  } else {
    abs(x - m) > k * md
  }
}

# -------- Cargar datos --------
if (!file.exists(input_path)) {
  stop(paste0("No se encontró el archivo: ", input_path))
}
df <- suppressMessages(readr::read_csv(input_path, show_col_types = FALSE))
if (!is.data.frame(df) || nrow(df) == 0) stop("El archivo está vacío o no es una tabla válida.")
nms <- names(df)

# Validar columnas del mapeo
for (ch in channels) {
  if (is.null(manual_map[[ch]])) stop("Falta manual_map para canal ", ch)
  wb <- manual_map[[ch]]
  if (!all(c(wb$white, wb$black) %in% nms)) {
    stop("En canal ", ch, " no existen columnas: ", wb$white, " y/o ", wb$black, " en el CSV.")
  }
}
if (!(id_col %in% nms)) stop("Columna ID '", id_col, "' no existe en el CSV.")

# Función para evaluar un canal
eval_channel <- function(dat, id_col, wcol, kcol, channel,
                         alpha_mahal, k_mad, min_contrast_rel, auto_normalize_255) {
  sub <- dat %>% 
    dplyr::select(tidyselect::all_of(c(id_col, wcol, kcol))) %>%
    dplyr::rename(ID = !!id_col, BLANCO = !!wcol, NEGRO = !!kcol)

  sub <- sub %>% dplyr::filter(!is.na(BLANCO), !is.na(NEGRO))
  if (nrow(sub) == 0) return(NULL)

  BL <- sub$BLANCO; NK <- sub$NEGRO

  if (auto_normalize_255) {
    max_val <- max(BL, NK, na.rm = TRUE)
    min_val <- min(BL, NK, na.rm = TRUE)
    if (is.finite(max_val) && max_val > 1.5 && max_val <= 300 && min_val >= 0) {
      BL <- BL / 255
      NK <- NK / 255
      sub$normalized_255 <- TRUE
    } else {
      sub$normalized_255 <- FALSE
    }
  } else {
    sub$normalized_255 <- FALSE
  }

  BLc <- pmin(pmax(BL, 0), 1)
  NKc <- pmin(pmax(NK, 0), 1)

  contraste      <- BLc - NKc
  rango_conjunto <- max(c(BLc, NKc), na.rm = TRUE) - min(c(BLc, NKc), na.rm = TRUE)
  contraste_rel  <- ifelse(rango_conjunto > 0, contraste / rango_conjunto, NA_real_)

  flag_blanco_mad <- mad_flag(BLc, k = k_mad)
  flag_negro_mad  <- mad_flag(NKc, k = k_mad)

  xy <- cbind(BLc, NKc)
  rob <- tryCatch(MASS::cov.rob(xy), error = function(e) NULL)
  if (is.null(rob) || any(!is.finite(rob$center)) || any(!is.finite(rob$cov))) {
    center <- colMeans(xy, na.rm = TRUE)
    covmat <- stats::cov(xy, use = "complete.obs")
  } else {
    center <- rob$center; covmat <- rob$cov
  }
  md2 <- tryCatch(stats::mahalanobis(xy, center = center, cov = covmat, tol = 1e-12),
                  error = function(e) rep(NA_real_, nrow(xy)))
  p_mahal <- stats::pchisq(md2, df = 2, lower.tail = FALSE)
  flag_mahal <- ifelse(is.na(p_mahal), FALSE, p_mahal < alpha_mahal)

  flag_contraste <- (contraste_rel < min_contrast_rel) | (BLc <= NKc)

  razones <- mapply(function(fb, fn, fm, fc) {
    rs <- c()
    if (fb) rs <- c(rs, "BLANCO fuera (MAD)")
    if (fn) rs <- c(rs, "NEGRO fuera (MAD)")
    if (fm) rs <- c(rs, "Bivariado fuera (Mahalanobis)")
    if (fc) rs <- c(rs, "Contraste bajo/invertido")
    if (length(rs) == 0) "OK" else paste(rs, collapse = " | ")
  }, fb = flag_blanco_mad, fn = flag_negro_mad, fm = flag_mahal, fc = flag_contraste, SIMPLIFY = TRUE)

  sub %>%
    dplyr::mutate(
      canal = channel,
      BLANCO_norm = BLc,
      NEGRO_norm  = NKc,
      contraste    = contraste,
      contraste_rel = contraste_rel,
      md2 = md2,
      p_mahal = p_mahal,
      flag_blanco_mad = flag_blanco_mad,
      flag_negro_mad  = flag_negro_mad,
      flag_mahal      = flag_mahal,
      flag_contraste  = flag_contraste,
      repetir_foto_canal = flag_blanco_mad | flag_negro_mad | flag_mahal | flag_contraste,
      razones = razones
    ) %>%
    dplyr::select(ID, canal, BLANCO_norm, NEGRO_norm, contraste, contraste_rel,
                  md2, p_mahal, flag_blanco_mad, flag_negro_mad, flag_mahal, flag_contraste,
                  repetir_foto_canal, razones, normalized_255)
}

# Ejecutar por canal con mapeo manual
results_list <- purrr::imap(manual_map, ~ eval_channel(df, id_col, .x$white, .x$black, .y,
                                          alpha_mahal, k_mad, min_contrast_rel, auto_normalize_255))
results_list <- results_list[!vapply(results_list, is.null, logical(1))]
if (length(results_list) == 0) stop("No hay resultados (quizá todos los canales estaban vacíos).")

res_long <- dplyr::bind_rows(results_list)

# Guardar long
out_long <- "repetir_fotos_flags_multicanal_long.csv"
readr::write_csv(res_long, out_long, na = "")

# Resumen por foto (agregado sobre canales)
resumen <- res_long %>%
  dplyr::group_by(ID) %>%
  dplyr::summarize(
    canales_evaluados = paste(sort(unique(canal)), collapse = ","),
    repetir_foto_global = any(repetir_foto_canal, na.rm = TRUE),
    canales_con_flag = paste(sort(unique(canal[repetir_foto_canal])), collapse = ","),
    razones_top = paste(head(unique(razones[repetir_foto_canal & razones != "OK"]), 3), collapse = " | ")
  ) %>%
  dplyr::ungroup()

out_res <- "repetir_fotos_flags_multicanal_resumen.csv"
readr::write_csv(resumen, out_res, na = "")

message("\n===== RESUMEN GLOBAL =====")
message("Total fotos: ", dplyr::n_distinct(res_long$ID))
message("Con al menos un canal con problema: ", sum(resumen$repetir_foto_global, na.rm = TRUE))

# ---------- Gráficos ----------
for (ch in unique(res_long$canal)) {
  sub <- res_long %>% dplyr::filter(canal == ch)
  if (nrow(sub) == 0) next
  p <- ggplot2::ggplot(sub, ggplot2::aes(x = NEGRO_norm, y = BLANCO_norm, shape = repetir_foto_canal)) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
    ggplot2::labs(x = paste0("NEGRO (norm) - ", ch), y = paste0("BLANCO (norm) - ", ch),
         title = paste0("Diagnóstico BLANCO vs NEGRO - Canal ", ch)) +
    ggplot2::theme_minimal()
  ggplot2::ggsave(paste0("diag_scatter_blanco_vs_negro_", ch, ".png"), p, width = 7, height = 5, dpi = 150)
}

prop_ch <- res_long %>%
  dplyr::group_by(canal) %>%
  dplyr::summarize(prop_repetir = mean(repetir_foto_canal, na.rm = TRUE), n = dplyr::n())

p2 <- ggplot2::ggplot(prop_ch, ggplot2::aes(x = canal, y = prop_repetir)) +
  ggplot2::geom_col() +
  ggplot2::labs(x = "Canal", y = "Proporción marcada", title = "Proporción de fotos marcadas por canal") +
  ggplot2::theme_minimal()
ggplot2::ggsave("diag_proporcion_flags_por_canal.png", p2, width = 7, height = 5, dpi = 150)

p3 <- ggplot2::ggplot(res_long, ggplot2::aes(x = canal, y = contraste_rel)) +
  ggplot2::geom_violin(trim = FALSE) +
  ggplot2::geom_jitter(width = 0.1, height = 0, alpha = 0.4) +
  ggplot2::geom_hline(yintercept = min_contrast_rel, linetype = 2) +
  ggplot2::labs(x = "Canal", y = "Contraste relativo", title = "Distribución de contraste relativo por canal") +
  ggplot2::theme_minimal()
ggplot2::ggsave("diag_contraste_relativo_por_canal.png", p3, width = 7, height = 5, dpi = 150)

message("\nArchivos generados:")
message(" - ", normalizePath(out_long))
message(" - ", normalizePath(out_res))
message(" - PNGs por canal + overview")
message("\nListo. Revisa los CSV y gráficos.")