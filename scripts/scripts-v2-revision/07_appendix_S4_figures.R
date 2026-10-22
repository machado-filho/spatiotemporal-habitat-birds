# =============================================================================
# 07 — APPENDIX S4, FIGURES S1–S3
# =============================================================================
# Run from the repository root after scripts 01 and 02 (paths are relative).
# Inputs (results/):
#   permutation_results.xlsx   sheets null_distribution, species_73, species_65
#   collinearity.xlsx          sheets selected_models, r_by_lags
# Outputs (results/):
#   FigureS1.png  null distribution of the number of species with P < 0.05
#                 in the selected model (999 permutations) and observed value
#   FigureS2.png  coefficient of each predictor in the full model vs. the model
#                 without the other habitat predictor (65 species)
#   FigureS3.png  median |r| between static and change predictors for each
#                 combination of static lag k and change interval l
# All text is >= 7 pt at final size.
# =============================================================================

library(openxlsx)
library(ggplot2)
library(cowplot)

pasta_base <- "results"

nulo  <- read.xlsx(file.path(pasta_base, "permutation_results.xlsx"), sheet = "null_distribution")
sp73  <- read.xlsx(file.path(pasta_base, "permutation_results.xlsx"), sheet = "species_73")
sp65  <- read.xlsx(file.path(pasta_base, "permutation_results.xlsx"), sheet = "species_65")
colin <- read.xlsx(file.path(pasta_base, "collinearity.xlsx"), sheet = "selected_models")
r_lags <- read.xlsx(file.path(pasta_base, "collinearity.xlsx"), sheet = "r_by_lags")

# ---- common style -------------------------------------------------------------------
tema <- theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 7, colour = "black"),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 7.5),
        legend.text = element_text(size = 7),
        plot.margin = margin(4, 6, 4, 4))
cor_obs <- "#B2182B"

# adds the "Figure S<n>" label in a strip above the plot
com_rotulo <- function(p, rotulo, altura = 0.93) {
  ggdraw() +
    draw_plot(p, x = 0, y = 0, width = 1, height = altura) +
    draw_label(rotulo, x = 0.01, y = 0.99, hjust = 0, vjust = 1, fontface = "bold", size = 10)
}

# ---- Figure S1: null distribution ------------------------------------------------------
# observed value = number of species (out of 73) whose selected model had P < 0.05
n_obs <- sum(sp73$any_P_below_0.05_obs, na.rm = TRUE)
v <- nulo$n_species_P_below_0.05
cat("Observed:", n_obs, "| null mean:", round(mean(v), 1), "| null range:", min(v), "-", max(v),
    "| P =", round((sum(v >= n_obs) + 1) / (length(v) + 1), 4), "\n")

s1 <- ggplot(nulo, aes(n_species_P_below_0.05)) +
  geom_histogram(breaks = seq(min(v) - 0.5, 74.5, 1), fill = "#9AA5B1", colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = n_obs, colour = cor_obs, linewidth = 0.8) +
  annotate("text", x = n_obs - 0.8, y = Inf, vjust = 1.3, hjust = 1, size = 2.6, colour = cor_obs,
           label = paste0("observed\n(", n_obs, ")"), lineheight = 0.9) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  coord_cartesian(xlim = c(25, 74)) +
  labs(x = "Number of species with a significant habitat\npredictor in the selected model (P < 0.05)",
       y = paste0("Frequency (", nrow(nulo), " permutations)")) +
  tema
ggsave(file.path(pasta_base, "FigureS1.png"), com_rotulo(s1, "Figure S1"),
       width = 11.5, height = 8.5, units = "cm", dpi = 300, bg = "white")

# ---- Figure S2: coefficients with and without the other predictor -------------------------
# 65 species with adequate residuals; black outline = significant after permutation + FDR
m <- merge(colin[colin$P_dispersion >= 0.05, ],
           sp65[, c("species", "P_FDR_static", "P_FDR_change")], by = "species")

monta_painel <- function(xc, yc, fdr, titulo) {
  d <- data.frame(x = m[[xc]], y = m[[yc]], r = abs(m$r_pearson), sig = m[[fdr]] < 0.05, painel = titulo)
  lim <- range(c(d$x, d$y), na.rm = TRUE) + c(-0.1, 0.1)
  list(pontos = d, limites = data.frame(x = lim, y = lim, painel = titulo))
}
pa <- monta_painel("coef_static_alone", "coef_static_full", "P_FDR_static", "(a) Static habitat amount")
pb <- monta_painel("coef_change_alone", "coef_change_full", "P_FDR_change", "(b) Habitat change")
pontos  <- rbind(pa$pontos, pb$pontos)
limites <- rbind(pa$limites, pb$limites)   # same range on x and y within each panel

s2 <- ggplot(pontos, aes(x, y)) +
  geom_blank(data = limites) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey45", linewidth = 0.3) +
  geom_point(data = pontos[!pontos$sig, ], aes(fill = r), shape = 21, size = 1.6,
             colour = "transparent", alpha = 0.85) +
  geom_point(data = pontos[pontos$sig, ], aes(fill = r), shape = 21, size = 2.3,
             colour = "black", stroke = 0.5) +
  scale_fill_viridis_c(name = "|r| between\npredictors", limits = c(0, 1)) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  facet_wrap(~ painel, scales = "free") +
  labs(x = "Coefficient without the other predictor", y = "Coefficient in the full model") +
  tema +
  theme(aspect.ratio = 1,
        strip.background = element_blank(),
        strip.text = element_text(size = 8, hjust = 0, face = "bold"),
        panel.spacing = unit(0.8, "cm"))
ggsave(file.path(pasta_base, "FigureS2.png"), com_rotulo(s2, "Figure S2"),
       width = 18, height = 9, units = "cm", dpi = 300, bg = "white")
cat("Figure S2: significant species — static:", sum(pa$pontos$sig, na.rm = TRUE),
    "| change:", sum(pb$pontos$sig, na.rm = TRUE), "\n")

# ---- Figure S3: median |r| by combination of k and l ----------------------------------------
s3 <- ggplot(r_lags, aes(factor(change_interval_l), factor(static_lag_k), fill = median_abs_r)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", median_abs_r),
                colour = ifelse(median_abs_r > 0.6, "dark", "light")), size = 2.6, show.legend = FALSE) +
  scale_colour_manual(values = c(dark = "black", light = "white")) +
  scale_fill_viridis_c(name = "Median |r|", limits = c(0, 1)) +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  labs(x = "Habitat change interval, l (years)", y = "Static habitat lag, k (years)") +
  tema + theme(axis.line = element_blank(), aspect.ratio = 1)
ggsave(file.path(pasta_base, "FigureS3.png"), com_rotulo(s3, "Figure S3"),
       width = 11, height = 9, units = "cm", dpi = 300, bg = "white")
