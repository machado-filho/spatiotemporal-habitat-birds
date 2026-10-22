# =============================================================================
# 02 — COLLINEARITY BETWEEN STATIC HABITAT AMOUNT AND HABITAT CHANGE
# =============================================================================
# Run from the repository root (paths are relative):
#   data/habitat_static/  dados_absolutos_<k>anos*.xlsx  (k = 0, 2, 3, 4, 5)
#   data/habitat_change/  dados_delta<l>anos*.xlsx       (l = 1, ..., 5)
#   results/              outputs written by the scripts
#
# Stand-alone script: reads the habitat files, joins habitat by the
# site + year key (as in 01_models_and_permutation.R) and, for the selected
# model of each species, computes:
#   - Pearson correlation between static and change predictors, and VIF;
#   - coefficient of static habitat in the full model vs. a model with static
#     habitat only;
#   - coefficient of habitat change in the full model vs. a model with
#     habitat change only.
# It also computes the static x change correlation for every combination of
# lags, intervals and radii.
# Requires results/selected_models.xlsx (script 01).
#
# Outputs (results/):
#   collinearity.xlsx               selected_models, r_by_lags, r_by_radius,
#                                   r_all_combinations
#   collinearity_coefficients.png   static coefficient, full vs. alone
#   collinearity_by_lags.png        median |r| for each combination of k and l
# =============================================================================

library(openxlsx)
library(dplyr)
library(stringr)
library(car)
library(ggplot2)

pasta_base  <- "results"
pasta_abs   <- file.path("data", "habitat_static")
pasta_delta <- file.path("data", "habitat_change")

# selected models (output of script 01)
arquivo_melhores <- "selected_models.xlsx"
aba_melhores     <- "all_species"      # "all_species" (73) or "residuals_ok" (65)

escalas_espaciais <- c(200, 600, 1000, 2000, 5000, 10000, 20000)
anos_estaticos <- c(0, 2, 3, 4, 5)
anos_delta     <- 1:5

# -----------------------------------------------------------------------------
# Functions (same as in script 01)
# -----------------------------------------------------------------------------
anos_do_arquivo <- function(x) as.numeric(str_extract(x, "[0-9]+(?=anos)"))

ler_planilha <- function(pasta, prefixo, ano) {
  arq <- list.files(pasta, pattern = "\\.xlsx$")
  arq <- arq[str_detect(arq, prefixo) & anos_do_arquivo(arq) %in% ano]
  if (length(arq) != 1) stop("Expected 1 file for ", prefixo, " year ", ano,
                             ", found: ", paste(arq, collapse = ", "))
  read.xlsx(file.path(pasta, arq))
}

criar_chave <- function(tab) {
  paste(round(tab$Longitude_x, 5), round(tab$Latitude_y, 5), tab$Year_finish2, sep = "_")
}

habitat_por_local <- function(tab, tipo, ano) {
  tab$chave <- criar_chave(tab)
  if (tipo == "estatico") {
    cols <- paste0("buffer_size", escalas_espaciais, "_absolute_Year_finish2_minus_", ano, "years")
  } else {
    cols <- paste0("buffer_size", escalas_espaciais, "_delta_Year_finish2-year")
  }
  h <- tab[, c("chave", cols)]
  h <- h[!duplicated(h$chave), ]
  colnames(h) <- c("chave", as.character(escalas_espaciais))
  h
}

padronizar <- function(x) as.numeric(scale(x))

# -----------------------------------------------------------------------------
# Reading the data and joining by site + year
# -----------------------------------------------------------------------------
base <- ler_planilha(pasta_abs, "dados_absolutos_", 0) %>%
  select(Species, Longitude_x, Latitude_y, Year_finish2, PE_PAI, PE_total_effort)
base$chave <- criar_chave(base)

hab.estatico <- list()
for (a in anos_estaticos) {
  h <- habitat_por_local(ler_planilha(pasta_abs, "dados_absolutos_", a), "estatico", a)
  hab.estatico[[as.character(a)]] <- h[match(base$chave, h$chave), as.character(escalas_espaciais)]
}
hab.delta <- list()
for (a in anos_delta) {
  h <- habitat_por_local(ler_planilha(pasta_delta, "dados_delta", a), "delta", a)
  hab.delta[[as.character(a)]] <- h[match(base$chave, h$chave), as.character(escalas_espaciais)]
}

montar_dados <- function(especie, ano.est, ano.del, escala) {
  idx <- which(base$Species == especie)
  df <- data.frame(PE_PAI = base$PE_PAI[idx],
                   estatico = hab.estatico[[as.character(ano.est)]][idx, as.character(escala)],
                   delta = hab.delta[[as.character(ano.del)]][idx, as.character(escala)],
                   PE_total_effort = base$PE_total_effort[idx])
  df <- na.omit(df)
  df$PE_total_effort[df$PE_total_effort == 0] <- 0.001
  df
}

ajustar <- function(formula, df) glm(formula, family = Gamma(link = "log"), data = df)
coef_p  <- function(mod, termo) {
  s <- coef(summary(mod))
  c(coef = s[termo, "Estimate"], p = s[termo, grep("^Pr", colnames(s))])
}

# quick check of the alignment (expected > 0.95)
cat("Correlation habitat year 0 vs year k (5000 m):\n")
print(sapply(as.character(anos_estaticos), function(a)
  round(cor(hab.estatico[["0"]][, "5000"], hab.estatico[[a]][, "5000"], use = "complete.obs"), 3)))

# -----------------------------------------------------------------------------
# 1) Selected models
# -----------------------------------------------------------------------------
melhores <- read.xlsx(file.path(pasta_base, arquivo_melhores), sheet = aba_melhores)

lista_res <- list()
for (i in 1:nrow(melhores)) {
  sp <- melhores$species[i]
  df <- montar_dados(sp, melhores$static_lag_k[i], melhores$change_interval_l[i], melhores$radius_m[i])
  if (nrow(df) < 5 || sd(df$estatico) == 0 || sd(df$delta) == 0) next
  df$estatico <- padronizar(df$estatico)
  df$delta    <- padronizar(df$delta)

  res <- tryCatch({
    m.conj  <- ajustar(PE_PAI ~ estatico + delta + PE_total_effort, df)   # full model
    m.est   <- ajustar(PE_PAI ~ estatico + PE_total_effort, df)           # static only
    m.delta <- ajustar(PE_PAI ~ delta + PE_total_effort, df)              # change only
    ce.conj <- coef_p(m.conj, "estatico"); cd.conj <- coef_p(m.conj, "delta")
    ce.so   <- coef_p(m.est, "estatico");  cd.so   <- coef_p(m.delta, "delta")
    v <- vif(m.conj)
    data.frame(
      species = sp, P_dispersion = melhores$P_dispersion[i],
      static_lag_k = melhores$static_lag_k[i], change_interval_l = melhores$change_interval_l[i],
      radius_m = melhores$radius_m[i],
      n = nrow(df),
      r_pearson = cor(df$estatico, df$delta),
      VIF_static = unname(v["estatico"]), VIF_change = unname(v["delta"]),
      coef_static_full  = unname(ce.conj["coef"]), P_static_full  = unname(ce.conj["p"]),
      coef_static_alone = unname(ce.so["coef"]),   P_static_alone = unname(ce.so["p"]),
      coef_change_full  = unname(cd.conj["coef"]), P_change_full  = unname(cd.conj["p"]),
      coef_change_alone = unname(cd.so["coef"]),   P_change_alone = unname(cd.so["p"]),
      # check: should match the coefficients in selected_models.xlsx
      coef_static_selected = melhores$coef_static[i], coef_change_selected = melhores$coef_change[i]
    )
  }, error = function(e) { message("Error in ", sp, ": ", conditionMessage(e)); NULL })
  lista_res[[sp]] <- res
}
colin <- do.call(rbind, lista_res); rownames(colin) <- NULL

classif <- function(coef, p) ifelse(p >= 0.05, "ns", ifelse(coef > 0, "positive", "negative"))
colin <- colin %>%
  mutate(collinear = abs(r_pearson) > 0.7 | pmax(VIF_static, VIF_change) > 5,
         static_full  = classif(coef_static_full, P_static_full),
         static_alone = classif(coef_static_alone, P_static_alone),
         change_full  = classif(coef_change_full, P_change_full),
         change_alone = classif(coef_change_alone, P_change_alone),
         matches_selected_model = abs(coef_static_full - coef_static_selected) < 1e-4 &
                                  abs(coef_change_full - coef_change_selected) < 1e-4)

# -----------------------------------------------------------------------------
# 2) Static x change correlation for every combination (no models)
# -----------------------------------------------------------------------------
lista_cor <- list(); cont <- 1
for (sp in unique(base$Species)) for (a in anos_estaticos) for (d in anos_delta) for (e in escalas_espaciais) {
  df <- montar_dados(sp, a, d, e)
  if (nrow(df) < 5 || sd(df$estatico) == 0 || sd(df$delta) == 0) next
  lista_cor[[cont]] <- data.frame(species = sp, static_lag_k = a, change_interval_l = d, radius_m = e,
                                  r = cor(df$estatico, df$delta))
  cont <- cont + 1
}
cor_todas <- do.call(rbind, lista_cor)
cor_por_anos <- cor_todas %>% group_by(static_lag_k, change_interval_l) %>%
  summarise(median_abs_r = median(abs(r)), prop_abs_r_above_0.7 = mean(abs(r) > 0.7), .groups = "drop")
cor_por_escala <- cor_todas %>% group_by(radius_m) %>%
  summarise(median_abs_r = median(abs(r)), prop_abs_r_above_0.7 = mean(abs(r) > 0.7), .groups = "drop")

# -----------------------------------------------------------------------------
# 3) Summary
# -----------------------------------------------------------------------------
ok <- colin$P_dispersion >= 0.05
cat("\n===== SELECTED MODELS — adequate residuals (", sum(ok), "species ) =====\n")
cat("Coefficients matching selected_models.xlsx:", sum(colin$matches_selected_model[ok]), "of", sum(ok), "\n")
cat("Median |r|:", round(median(abs(colin$r_pearson[ok])), 2),
    "| maximum |r|:", round(max(abs(colin$r_pearson[ok])), 2), "\n")
cat("Maximum VIF:", round(max(pmax(colin$VIF_static, colin$VIF_change)[ok]), 2), "\n")
cat("Species with |r| > 0.7 or VIF > 5:", sum(colin$collinear[ok]), "\n")
cat("\nStatic: alone (rows) x full model (columns)\n")
print(table(Alone = colin$static_alone[ok], Full = colin$static_full[ok]))
cat("\nChange: alone (rows) x full model (columns)\n")
print(table(Alone = colin$change_alone[ok], Full = colin$change_full[ok]))
cat("\nCorrelation by combination of k and l (all species and radii):\n")
print(as.data.frame(cor_por_anos))
print(as.data.frame(cor_por_escala))

cat("\n===== SOFTWARE VERSIONS (for the Methods) =====\n")
print(R.version.string); print(packageVersion("DHARMa")); print(packageVersion("car"))

# figures (the versions used in Appendix S4 are made by script 07)
g1 <- ggplot(colin[ok, ], aes(x = coef_static_alone, y = coef_static_full, color = abs(r_pearson))) +
  geom_hline(yintercept = 0, color = "grey60") + geom_vline(xintercept = 0, color = "grey60") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(size = 2.5) +
  scale_color_viridis_c(name = "|r|", limits = c(0, 1)) +
  labs(x = "Static habitat coefficient (model without habitat change)",
       y = "Static habitat coefficient (full model)") +
  theme_classic()
g2 <- ggplot(cor_por_anos, aes(x = factor(change_interval_l), y = factor(static_lag_k), fill = median_abs_r)) +
  geom_tile() + geom_text(aes(label = sprintf("%.2f", median_abs_r)), size = 3.5) +
  scale_fill_viridis_c(name = "Median |r|", limits = c(0, 1)) +
  labs(x = "Habitat change interval (years)", y = "Static habitat lag (years)") + theme_minimal()
print(g1); print(g2)
ggsave(file.path(pasta_base, "collinearity_coefficients.png"), g1, width = 16, height = 13, units = "cm", dpi = 300)
ggsave(file.path(pasta_base, "collinearity_by_lags.png"), g2, width = 14, height = 11, units = "cm", dpi = 300)

write.xlsx(list(selected_models = colin, r_by_lags = cor_por_anos,
                r_by_radius = cor_por_escala, r_all_combinations = cor_todas),
           file.path(pasta_base, "collinearity.xlsx"))
