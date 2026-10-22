# =============================================================================
# 05 — APPENDIX S3, FIGURES S1–S73: selected model of every species
# =============================================================================
# Run from the repository root (paths are relative):
#   data/habitat_static/  dados_absolutos_<k>anos*.xlsx  (k = 0, 2, 3, 4, 5)
#   data/habitat_change/  dados_delta<l>anos*.xlsx       (l = 1, ..., 5)
#   results/              selected_models.xlsx and model_parameters.xlsx
#                         (script 01) and outputs
#
# One figure per species with the same four panels as Figure 4:
#   (a) partial effect of static habitat amount; (b) partial effect of habitat
#   change; (c) QQ plot of the scaled residuals (DHARMa) with the P-value of
#   the dispersion test; (d) dAIC profile across buffer radii.
# The dispersion P-value shown is the one stored in selected_models.xlsx
# (column P_dispersion), i.e. the value used to exclude 8 species from
# species-level inference.
#
# Outputs (results/):
#   Appendix_S3_figures.pdf          73 A4 pages, one per species
#   Appendix_S3_png/FigureS<n>.png   each figure as PNG (300 dpi, 18 x 22.5 cm)
#   Appendix_S3_figure_list.xlsx     figure number x species
# =============================================================================

library(openxlsx)
library(dplyr)
library(stringr)
library(ggplot2)
library(visreg)
library(DHARMa)
library(cowplot)

pasta_base  <- "results"
pasta_abs   <- file.path("data", "habitat_static")
pasta_delta <- file.path("data", "habitat_change")
arquivo_melhores <- "selected_models.xlsx"    # selected model per species (script 01)
arquivo_todas    <- "model_parameters.xlsx"   # all 175 models per species (sheet "all_models")

escalas_espaciais <- c(200, 600, 1000, 2000, 5000, 10000, 20000)
anos_estaticos <- c(0, 2, 3, 4, 5)
anos_delta     <- 1:5

# ---- reading and joining by site + year (as in script 01) ----------------------
anos_do_arquivo <- function(x) as.numeric(str_extract(x, "[0-9]+(?=anos)"))
ler_planilha <- function(pasta, prefixo, ano) {
  arq <- list.files(pasta, pattern = "\\.xlsx$")
  arq <- arq[str_detect(arq, prefixo) & anos_do_arquivo(arq) %in% ano]
  if (length(arq) != 1) stop("Expected 1 file for ", prefixo, " year ", ano)
  read.xlsx(file.path(pasta, arq))
}
criar_chave <- function(tab) paste(round(tab$Longitude_x, 5), round(tab$Latitude_y, 5), tab$Year_finish2, sep = "_")
habitat_por_local <- function(tab, tipo, ano) {
  tab$chave <- criar_chave(tab)
  cols <- if (tipo == "estatico") paste0("buffer_size", escalas_espaciais, "_absolute_Year_finish2_minus_", ano, "years") else
    paste0("buffer_size", escalas_espaciais, "_delta_Year_finish2-year")
  h <- tab[, c("chave", cols)]; h <- h[!duplicated(h$chave), ]
  colnames(h) <- c("chave", as.character(escalas_espaciais)); h
}
base <- ler_planilha(pasta_abs, "dados_absolutos_", 0) %>%
  select(Species, Longitude_x, Latitude_y, Year_finish2, PE_PAI, PE_total_effort)
base$chave <- criar_chave(base)
hab.estatico <- list(); hab.delta <- list()
for (a in anos_estaticos) { h <- habitat_por_local(ler_planilha(pasta_abs, "dados_absolutos_", a), "estatico", a)
  hab.estatico[[as.character(a)]] <- h[match(base$chave, h$chave), as.character(escalas_espaciais)] }
for (a in anos_delta) { h <- habitat_por_local(ler_planilha(pasta_delta, "dados_delta", a), "delta", a)
  hab.delta[[as.character(a)]] <- h[match(base$chave, h$chave), as.character(escalas_espaciais)] }

montar_dados <- function(especie, ano.est, ano.del, escala) {
  idx <- which(base$Species == especie)
  df <- data.frame(PE_PAI = base$PE_PAI[idx],
                   Static = hab.estatico[[as.character(ano.est)]][idx, as.character(escala)],
                   Change = hab.delta[[as.character(ano.del)]][idx, as.character(escala)],
                   Effort = base$PE_total_effort[idx])
  df <- na.omit(df); df$Effort[df$Effort == 0] <- 0.001; df
}

melhores <- read.xlsx(file.path(pasta_base, arquivo_melhores), sheet = "all_species")
todas    <- read.xlsx(file.path(pasta_base, arquivo_todas), sheet = "all_models")

# ---- theme and colours ------------------------------------------------------------
tema <- theme_classic(base_size = 8) +
  theme(plot.title = element_text(size = 8, face = "italic", hjust = 0),
        plot.subtitle = element_text(size = 7),
        axis.title = element_text(size = 7), axis.text = element_text(size = 7, color = "black"),
        plot.margin = margin(3, 9, 3, 3))   # right margin keeps the last x label (e.g. 0.75) from being cut
cor_estatico <- "#2166AC"; cor_delta <- "#1B7837"

# ---- the four panels of one species -------------------------------------------------
paineis_especie <- function(sp) {
  m <- melhores[melhores$species == sp, ]
  df <- montar_dados(sp, m$static_lag_k, m$change_interval_l, m$radius_m)
  # model on the original scale of the predictors (same fit and AIC as the standardized model)
  mod <- glm(PE_PAI ~ Static + Change + Effort, family = Gamma(link = "log"), data = df)
  nome <- gsub("_", " ", sp)

  p1 <- visreg(mod, "Static", gg = TRUE, plot = TRUE,
               line = list(col = cor_estatico, size = 0.8),
               fill = list(fill = "grey85"),
               points = list(size = 0.9, col = "grey35", alpha = 0.6)) +
    labs(x = "Forest cover (%)",
         y = "Partial effect (log abundance)", title = nome,
         subtitle = paste0("radius = ", m$radius_m / 1000, " km; k = ", m$static_lag_k, " yr")) + tema

  p2 <- visreg(mod, "Change", gg = TRUE, plot = TRUE,
               line = list(col = cor_delta, size = 0.8),
               fill = list(fill = "grey85"),
               points = list(size = 0.9, col = "grey35", alpha = 0.6)) +
    labs(x = "Habitat change (% points)",
         y = "Partial effect (log abundance)", title = " ",
         subtitle = paste0("l = ", m$change_interval_l, " yr")) + tema

  set.seed(1)
  sim <- simulateResiduals(mod, n = 1000)
  p_disp <- m$P_dispersion   # same P-value used to exclude species (selected_models.xlsx)
  res <- sort(sim$scaledResiduals)
  qq <- data.frame(esperado = (seq_along(res) - 0.5) / length(res), observado = res)
  p3 <- ggplot(qq, aes(esperado, observado)) +
    geom_abline(slope = 1, intercept = 0, color = "#B2182B", linewidth = 0.5) +
    geom_point(size = 0.9, color = "grey30") +
    annotate("text", x = 0.02, y = 0.98, hjust = 0, vjust = 1, size = 2.5,
             label = paste0("dispersion P = ", formatC(p_disp, format = "f", digits = 2)),
             color = ifelse(p_disp >= 0.05, cor_delta, "#B2182B")) +
    scale_x_continuous(breaks = c(0, 0.5, 1)) +   # fewer ticks: 0.00-1.00 in 0.25 steps crowded at 7 pt
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "Expected", y = "Observed (scaled residuals)", title = " ", subtitle = " ") + tema

  # dAIC profile: one grey line per combination of k and l, selected combination in black
  a <- todas[todas$species == sp, ]
  a$dAIC <- a$AIC - min(a$AIC)
  a$comb <- paste(a$static_lag_k, a$change_interval_l)
  a$sel <- a$static_lag_k == m$static_lag_k & a$change_interval_l == m$change_interval_l
  p4 <- ggplot(a, aes(radius_m, dAIC, group = comb)) +
    geom_line(data = a[!a$sel, ], color = "grey75", linewidth = 0.3) +
    geom_line(data = a[a$sel, ], color = "black", linewidth = 0.8) +
    geom_point(data = a[a$sel, ], color = "black", size = 0.9) +
    geom_vline(xintercept = m$radius_m, linetype = "dashed", linewidth = 0.4) +
    geom_hline(yintercept = 2, linetype = "dotted", linewidth = 0.3, color = "grey40") +
    scale_x_log10(breaks = c(200, 1000, 5000, 20000), labels = c("0.2", "1", "5", "20")) +
    labs(x = "Buffer radius (km)", y = expression(Delta * "AIC"), title = " ",
         subtitle = paste0("k = ", m$static_lag_k, " yr, l = ", m$change_interval_l, " yr")) + tema
  list(p1, p2, p3, p4)
}


# ---- one figure per species (Figure S1 ... S73) ----------------------------------
# alphabetical order, as in the table of Appendix S3
todas_sp <- sort(melhores$species)
pasta_png <- file.path(pasta_base, "Appendix_S3_png")
dir.create(pasta_png, showWarnings = FALSE)

montar_pagina <- function(n, sp) {
  m  <- melhores[melhores$species == sp, ]
  pp <- paineis_especie(sp)
  pp <- lapply(pp, function(p) p + labs(title = NULL))   # species name goes in the caption
  painel <- plot_grid(plotlist = pp, ncol = 2, align = "hv", axis = "tblr",
                      labels = c("a", "b", "c", "d"), label_size = 9, label_fontface = "bold")

  # line 1: "Figure S<n>." in bold + species name in italics
  titulo <- bquote(bold(.(paste0("Figure S", n, ".")))~italic(.(gsub("_", " ", sp))))
  # following lines: caption text, wrapped to the page width
  excl <- if (isTRUE(m$P_dispersion < 0.05)) " This species was excluded from species-level inference (dispersion test P < 0.05)." else ""
  texto <- paste0("Selected model: radius = ", format(m$radius_m, big.mark = ","), " m, static lag k = ",
                  m$static_lag_k, " yr, change interval l = ", m$change_interval_l, " yr (n = ", m$n, " sites). ",
                  "(a) Partial effect of static habitat amount and (b) of habitat change (lines: fitted ",
                  "effect; bands: 95% confidence interval; points: partial residuals). (c) QQ plot of ",
                  "simulated scaled residuals (DHARMa) with the P-value of the dispersion test. ",
                  "(d) \u0394AIC across buffer radii for each combination of k and l (black: selected ",
                  "combination; dashed line: selected radius; dotted line: \u0394AIC = 2).", excl)
  texto <- paste(strwrap(texto, width = 115), collapse = "\n")

  legenda <- ggdraw() +
    draw_label(titulo, x = 0, y = 1, hjust = 0, vjust = 1, size = 9) +
    draw_label(texto, x = 0, y = 0.84, hjust = 0, vjust = 1, size = 8, lineheight = 1.1)

  plot_grid(legenda, painel, ncol = 1, rel_heights = c(0.22, 1))
}

# 1) build the pages
paginas <- list()
lista_fig <- data.frame(figure = character(), species = character(), stringsAsFactors = FALSE)
n <- 0
for (sp in todas_sp) {
  pg <- tryCatch(montar_pagina(n + 1, sp),
                 error = function(e) { message("Error in ", sp, ": ", conditionMessage(e)); NULL })
  if (is.null(pg)) next
  n <- n + 1
  paginas[[n]] <- pg
  lista_fig[n, ] <- c(paste0("Figure S", n), gsub("_", " ", sp))
  cat("Figure S", n, ": ", sp, "\n", sep = "")
}

# 2) PDF: one A4 page per species (~1.5 cm side margins)
pdf(file.path(pasta_base, "Appendix_S3_figures.pdf"), width = 8.27, height = 11.69)
for (pg in paginas) print(ggdraw() + draw_plot(pg, x = 0.07, y = 0.20, width = 0.86, height = 0.75))
dev.off()

# 3) PNG: figure only (18 x 22.5 cm), to insert into Word (Insert > Pictures)
for (i in seq_along(paginas)) {
  ggsave(file.path(pasta_png, paste0("FigureS", i, ".png")), paginas[[i]],
         width = 18, height = 22.5, units = "cm", dpi = 300, bg = "white")
}

write.xlsx(lista_fig, file.path(pasta_base, "Appendix_S3_figure_list.xlsx"))
cat("\nTotal number of figures:", n, "(should be 73)\n")
