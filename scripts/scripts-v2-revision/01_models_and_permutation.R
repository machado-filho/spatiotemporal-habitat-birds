# =============================================================================
# 01 — STATIC + CHANGE HABITAT MODELS ACROSS SPATIAL AND TEMPORAL SCALES,
#      MODEL SELECTION, FIGURE 3 AND PERMUTATION TEST
# =============================================================================
# Run from the repository root (paths are relative):
#   data/habitat_static/  dados_absolutos_<k>anos*.xlsx  (k = 0, 2, 3, 4, 5)
#   data/habitat_change/  dados_delta<l>anos*.xlsx       (l = 1, ..., 5)
#   results/              outputs written by the scripts
#
# For each species, a gamma GLM (log link)
#   abundance ~ z[H(t0 - k)] + z[dH(l)] + sampling effort
# is fitted for every combination of static lag k (0, 2, 3, 4, 5 years),
# change interval l (1-5 years) and buffer radius (7 radii, 200-20,000 m),
# i.e. 175 models per species. The model with the lowest AIC is selected.
# The whole selection procedure is then repeated on 999 permutations of the
# habitat predictors to obtain permutation P-values, corrected with FDR.
#
# Notes on data handling:
#  - Habitat values of each year (static and change) are joined to the
#    abundance table by a site + year key (Longitude_x, Latitude_y,
#    Year_finish2), never by row order.
#  - Abundance (PE_PAI) and sampling effort always come from the same base
#    table (static habitat file for k = 0).
#  - Static lags: 0, 2, 3, 4, 5 (there is no file for k = 1). Change
#    intervals: 1-5. A change interval of 0 is not used, because it is an
#    exact linear transformation of current habitat (r = 1), not a change.
#  - A single function (montar_dados) builds the data for any
#    species x lag x interval x radius and is used in every loop.
#  - Variable and function names are in Portuguese; comments, console
#    messages and output tables are in English.
#
# Outputs (results/):
#   model_parameters.xlsx       all_models: every model of every species;
#                               selected_models: selected model per species (73)
#   selected_models.xlsx        residuals_ok (65 species) and all_species (73)
#   permutation_p_values.rds    raw permutation P-values
#   permutation_results.xlsx    species_65, species_73, null_distribution
#   Figure3.tif                 standardized coefficients of the selected models
# =============================================================================

library(openxlsx)
library(stringr)
library(ggplot2)
library(dplyr)
library(DHARMa)
library(car)
library(visreg)
library(cowplot)

pasta_base  <- "results"
pasta_abs   <- file.path("data", "habitat_static")
pasta_delta <- file.path("data", "habitat_change")
dir.create(pasta_base, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# SCALES
# -----------------------------------------------------------------------------
# Spatial: roughly geometric progression, from the neighbourhood of the point
# (hundreds of metres) to the landscape scale (tens of km).
escalas_espaciais <- c(200, 600, 1000, 2000, 5000, 10000, 20000)

# Temporal: up to 5 years (~1-2 generations of small passerines;
# Bird et al. 2020, Conserv. Biol. 34:1252-1261).
anos_estaticos <- c(0, 2, 3, 4, 5)   # there is no file for k = 1
anos_delta     <- 1:5                # change interval 0 excluded (see header)

# To run only the four example species, set them here (NULL = all species):
especies_exemplo <- NULL
#especies_exemplo <- c("Xenops_rutilans", "Camptostoma_obsoletum", "Cnemotriccus_fuscatus", "Trichothraupis_melanops")

mostrar_graficos_AIC  <- FALSE   # TRUE = show the AIC profile of each species (press Enter to continue)
mostrar_diagnosticos  <- FALSE   # TRUE = show visreg + DHARMa plots of each species (press Enter to continue)

# -----------------------------------------------------------------------------
# HELPER FUNCTIONS
# -----------------------------------------------------------------------------
anos_do_arquivo <- function(x) as.numeric(str_extract(x, "[0-9]+(?=anos)"))

# reads the file of one year (the prefix distinguishes static from change files)
ler_planilha <- function(pasta, prefixo, ano) {
  arq <- list.files(pasta, pattern = "\\.xlsx$")
  arq <- arq[str_detect(arq, prefixo) & anos_do_arquivo(arq) %in% ano]
  if (length(arq) != 1) stop("Expected 1 file for ", prefixo, " year ", ano,
                             ", found: ", paste(arq, collapse = ", "))
  read.xlsx(file.path(pasta, arq))
}

# site + sampling-year key
criar_chave <- function(tab) {
  paste(round(tab$Longitude_x, 5), round(tab$Latitude_y, 5), tab$Year_finish2, sep = "_")
}

# habitat per site + year (one row per key), only for the radii used,
# with columns renamed to the radius value ("200", "600", ...)
habitat_por_local <- function(tab, tipo, ano) {
  tab$chave <- criar_chave(tab)
  if (tipo == "estatico") {
    cols <- paste0("buffer_size", escalas_espaciais, "_absolute_Year_finish2_minus_", ano, "years")
  } else {
    cols <- paste0("buffer_size", escalas_espaciais, "_delta_Year_finish2-year")
  }
  faltam <- setdiff(cols, colnames(tab))
  if (length(faltam) > 0) stop("Columns not found (", tipo, " ", ano, "): ", paste(faltam, collapse = ", "))
  h <- tab[, c("chave", cols)]
  # habitat of a given site + year must be identical in all rows
  n_chaves <- length(unique(h$chave))
  n_combinacoes <- nrow(distinct(h))
  if (n_combinacoes > n_chaves) {
    warning(tipo, " ", ano, ": ", n_combinacoes - n_chaves,
            " site-years with different habitat values among rows! Using the first row.")
  }
  h <- h[!duplicated(h$chave), ]
  colnames(h) <- c("chave", as.character(escalas_espaciais))
  h
}

# standardization (mean 0, SD 1) of the habitat predictors: does not change
# AIC or P-values, only makes static and change coefficients comparable
padronizar_preditores <- function(df) {
  df$preditor.current <- as.numeric(scale(df$preditor.current))
  df$preditor.delta   <- as.numeric(scale(df$preditor.delta))
  df
}

# -----------------------------------------------------------------------------
# READING AND JOINING THE DATA (done once)
# -----------------------------------------------------------------------------
# base table: abundance and effort (static habitat file, k = 0)
base <- ler_planilha(pasta_abs, "dados_absolutos_", 0) %>%
  select(Species, Longitude_x, Latitude_y, Year_finish2, PE_PAI, PE_total_effort)
base$chave <- criar_chave(base)

# habitat of each year, aligned to the rows of the base table by the key
hab.estatico <- list()
for (a in anos_estaticos) {
  h <- habitat_por_local(ler_planilha(pasta_abs, "dados_absolutos_", a), "estatico", a)
  hab.estatico[[as.character(a)]] <- h[match(base$chave, h$chave), as.character(escalas_espaciais)]
}
hab.delta <- list()
for (a in c(0, anos_delta)) {   # change interval 0 is read only for the checks below
  h <- habitat_por_local(ler_planilha(pasta_delta, "dados_delta", a), "delta", a)
  hab.delta[[as.character(a)]] <- h[match(base$chave, h$chave), as.character(escalas_espaciais)]
}

especies <- unique(base$Species)
if (!is.null(especies_exemplo)) especies <- especies_exemplo

# builds the data of one species x static lag x change interval x radius
# perm = shuffled order of the habitat predictors (permutation test only)
montar_dados <- function(especie, ano.est, ano.del, escala, perm = NULL) {
  idx <- which(base$Species == especie)
  cur <- hab.estatico[[as.character(ano.est)]][idx, as.character(escala)]
  del <- hab.delta[[as.character(ano.del)]][idx, as.character(escala)]
  if (!is.null(perm)) { cur <- cur[perm]; del <- del[perm] }
  df <- data.frame(PE_PAI = base$PE_PAI[idx],
                   preditor.current = cur,
                   preditor.delta = del,
                   PE_total_effort = base$PE_total_effort[idx])
  df <- na.omit(df)
  df$PE_total_effort[df$PE_total_effort == 0] <- 0.001
  df
}

ajustar_glm <- function(df) {
  glm(PE_PAI ~ preditor.current + preditor.delta + PE_total_effort,
      family = Gamma(link = "log"), data = df)
}

# =============================================================================
# DATA CHECKS — run and inspect before fitting the models
# =============================================================================
cat("\n===== DATA CHECKS =====\n")
# (a) rows of the base table without habitat (key not found) in each year
cat("Rows without habitat (key not found):\n")
print(sapply(hab.estatico, function(m) sum(is.na(m[, 1]))))
print(sapply(hab.delta,    function(m) sum(is.na(m[, 1]))))

# (b) habitat in year 0 vs years 2-5 (expected > 0.95: forest cover changes little in a few years)
cat("\nCorrelation habitat year 0 vs year k (5000 m radius):\n")
print(sapply(as.character(anos_estaticos), function(a)
  round(cor(hab.estatico[["0"]][, "5000"], hab.estatico[[a]][, "5000"], use = "complete.obs"), 3)))

# (c) is change l equal to H(year 0) - H(year l)? (expected ~1 for the years with a static file)
cat("\nCorrelation change l vs [H(0) - H(l)] (5000 m radius):\n")
for (a in intersect(anos_delta, anos_estaticos)) {
  dif <- hab.estatico[["0"]][, "5000"] - hab.estatico[[as.character(a)]][, "5000"]
  cat("  change", a, ":", round(cor(hab.delta[[as.character(a)]][, "5000"], dif, use = "complete.obs"), 3), "\n")
}

# (d) is change 0 identical to current habitat? (if r = 1 it is not a change -> excluded)
cat("\nCorrelation change 0 vs habitat year 0 (5000 m):",
    round(cor(hab.delta[["0"]][, "5000"], hab.estatico[["0"]][, "5000"], use = "complete.obs"), 3), "\n")

# (e) correlation between static and change predictors (expected low to moderate)
cat("\nCorrelation static vs change for each combination of k and l (5000 m):\n")
print(outer(anos_estaticos, anos_delta, Vectorize(function(a, d)
  round(cor(hab.estatico[[as.character(a)]][, "5000"], hab.delta[[as.character(d)]][, "5000"],
            use = "complete.obs"), 2))) %>%
  `dimnames<-`(list(paste0("k", anos_estaticos), paste0("l", anos_delta))))

# =============================================================================
# MAIN LOOP: all combinations
# =============================================================================
lista_AICs <- list()
contador <- 1

for (ano.est in anos_estaticos) {
  for (ano.del in anos_delta) {
    for (especie in especies) {
      for (escala in escalas_espaciais) {

        temp_df <- montar_dados(especie, ano.est, ano.del, escala)
        if (nrow(temp_df) < 5) next
        if (var(temp_df$preditor.current) == 0 || var(temp_df$preditor.delta) == 0) next
        temp_df <- padronizar_preditores(temp_df)

        mod.i <- tryCatch(ajustar_glm(temp_df), error = function(e) {
          message("Error: ", especie, " k ", ano.est, " l ", ano.del, " radius ", escala)
          NULL })
        if (is.null(mod.i)) next

        resul <- tryCatch(summary(mod.i), error = function(e) NULL)
        teste <- tryCatch(Anova(mod.i, type = "III"), error = function(e) NULL)
        if (is.null(resul) || is.null(teste)) next

        lista_AICs[[contador]] <- data.frame(
          static_lag_k = ano.est,
          change_interval_l = ano.del,
          radius_m = escala,
          n = nrow(temp_df),
          AIC = AIC(mod.i),
          species = especie,
          coef_change = resul$coefficients["preditor.delta", "Estimate"],
          coef_static = resul$coefficients["preditor.current", "Estimate"],
          P_change = teste["preditor.delta", "Pr(>Chisq)"],      # type III likelihood-ratio test
          P_static = teste["preditor.current", "Pr(>Chisq)"],
          pseudo_R2 = (resul$null.deviance - resul$deviance) / resul$null.deviance
        )
        contador <- contador + 1
      }
    }
  }
  print(paste("Static lag", ano.est, "done"))
}

df_AICsABS <- do.call(rbind, lista_AICs)
# the best model is chosen by AIC only (no prior filtering by P)
df_AICsABS$any_P_below_0.05 <- df_AICsABS$P_change < 0.05 | df_AICsABS$P_static < 0.05
df_AICsABS$lags_k_l <- paste0(df_AICsABS$static_lag_k, ".", df_AICsABS$change_interval_l)

# AIC profile per species (optional)
if (mostrar_graficos_AIC) {
  for (especie in unique(df_AICsABS$species)) {
    df_especie <- df_AICsABS[df_AICsABS$species == especie, ]
    melhor <- df_especie[which.min(df_especie$AIC), ]
    p <- ggplot(df_especie, aes(x = radius_m, y = AIC, color = factor(change_interval_l),
                                linetype = factor(static_lag_k),
                                group = interaction(static_lag_k, change_interval_l))) +
      geom_line(linewidth = .6) + scale_x_log10() + theme_classic() +
      labs(title = paste(especie, "- k", melhor$static_lag_k, "/ l", melhor$change_interval_l,
                         "/ radius", melhor$radius_m),
           color = "Change interval l", linetype = "Static lag k")
    print(p)
    readline(prompt = "Press [Enter] to see the next species...")
  }
}

# =============================================================================
# SELECTED (LOWEST-AIC) MODEL PER SPECIES
# =============================================================================
melhores_modelos <- df_AICsABS %>%
  group_by(species) %>%
  slice_min(order_by = AIC, n = 1, with_ties = FALSE) %>%
  ungroup()

# how well does the selected model stand out from the others?
destaque_AIC <- df_AICsABS %>%
  group_by(species) %>%
  arrange(AIC, .by_group = TRUE) %>%
  summarise(dAIC_second_best = AIC[2] - AIC[1],
            n_models_dAIC_below_2 = sum(AIC - min(AIC) < 2))
melhores_modelos <- left_join(melhores_modelos, destaque_AIC, by = "species")

# visual diagnostics (optional)
if (mostrar_diagnosticos) {
  for (i in 1:nrow(melhores_modelos)) {
    m <- melhores_modelos[i, ]
    temp_df <- padronizar_preditores(montar_dados(m$species, m$static_lag_k, m$change_interval_l, m$radius_m))
    mod <- ajustar_glm(temp_df)
    sim_res <- simulateResiduals(fittedModel = mod, n = 1000)
    layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE), heights = c(.8, 1))
    par(mar = c(4, 4, 3, 2) + 0.2)
    visreg(mod, "preditor.current", main = paste(m$species, "- static, k =", m$static_lag_k, "-", m$radius_m, "m"))
    visreg(mod, "preditor.delta",   main = paste(m$species, "- change, l =", m$change_interval_l, "-", m$radius_m, "m"))
    testDispersion(sim_res)
    layout(1)
    readline(prompt = " ")
  }
}

# confidence intervals and residual (dispersion) test of each selected model
nova_tabela <- melhores_modelos
nova_tabela$CI_low_change <- NA; nova_tabela$CI_high_change <- NA
nova_tabela$CI_low_static <- NA; nova_tabela$CI_high_static <- NA
nova_tabela$P_dispersion <- NA; nova_tabela$CI_method <- NA

for (i in 1:nrow(melhores_modelos)) {
  m <- melhores_modelos[i, ]
  temp_df <- padronizar_preditores(montar_dados(m$species, m$static_lag_k, m$change_interval_l, m$radius_m))
  mod1 <- ajustar_glm(temp_df)
  coefs <- coef(summary(mod1))

  # profile-likelihood CI; Wald CI if the profile fails
  conf_int <- tryCatch(Confint(mod1, exponentiate = FALSE, type = "LR"), error = function(e) NULL)
  nova_tabela$CI_method[i] <- "LR"
  if (is.null(conf_int) || any(is.na(conf_int[c("preditor.current", "preditor.delta"), 2:3]))) {
    t_crit <- qt(0.975, df.residual(mod1))
    conf_int <- cbind(Estimate = coefs[, "Estimate"],
                      coefs[, "Estimate"] - t_crit * coefs[, "Std. Error"],
                      coefs[, "Estimate"] + t_crit * coefs[, "Std. Error"])
    nova_tabela$CI_method[i] <- "Wald"
  }
  nova_tabela$coef_change[i] <- conf_int["preditor.delta", 1]
  nova_tabela$coef_static[i] <- conf_int["preditor.current", 1]
  nova_tabela$CI_low_change[i]  <- conf_int["preditor.delta", 2]
  nova_tabela$CI_high_change[i] <- conf_int["preditor.delta", 3]
  nova_tabela$CI_low_static[i]  <- conf_int["preditor.current", 2]
  nova_tabela$CI_high_static[i] <- conf_int["preditor.current", 3]
  nova_tabela$P_dispersion[i] <- testDispersion(simulateResiduals(fittedModel = mod1, n = 1000), plot = FALSE)$p.value
  # P-values: the type III LR tests from the main loop are kept (the same used in the selection)
  print(paste(i, "of", nrow(melhores_modelos), "-", m$species))
}
nova_tabela$estimate_within_CI <- with(nova_tabela,
  coef_change >= CI_low_change & coef_change <= CI_high_change &
  coef_static >= CI_low_static & coef_static <= CI_high_static)
print(table(nova_tabela$estimate_within_CI, useNA = "ifany"))   # should be TRUE only

melhores_modelos_todos <- nova_tabela                             # all species
melhores_modelos <- filter(nova_tabela, P_dispersion >= 0.05)     # adequate residuals only
write.xlsx(list(residuals_ok = melhores_modelos, all_species = melhores_modelos_todos),
           file.path(pasta_base, "selected_models.xlsx"))
# file cited in the manuscript: all combinations + selected models
write.xlsx(list(all_models = df_AICsABS, selected_models = melhores_modelos_todos),
           file.path(pasta_base, "model_parameters.xlsx"))

# =============================================================================
# DESCRIPTIVE RESULTS (nominal P < 0.05, before the permutation test)
# =============================================================================
N <- nrow(melhores_modelos)
sd_ <- melhores_modelos$P_change < 0.05
sc_ <- melhores_modelos$P_static < 0.05
cat("\nN species (adequate residuals):", N, "\n")
cat("Change  positive:", sum(sd_ & melhores_modelos$coef_change > 0),
    "| negative:", sum(sd_ & melhores_modelos$coef_change < 0), "\n")
cat("Static  positive:", sum(sc_ & melhores_modelos$coef_static > 0),
    "| negative:", sum(sc_ & melhores_modelos$coef_static < 0), "\n")
cat("No effect:", sum(!sd_ & !sc_), "| both:", sum(sd_ & sc_), "\n")
print(table(radius_m = melhores_modelos$radius_m))
print(table(change_interval_l = melhores_modelos$change_interval_l))
print(table(static_lag_k = melhores_modelos$static_lag_k))
cat("Species whose second-best model is within dAIC < 2:", sum(melhores_modelos$dAIC_second_best < 2), "of", N, "\n")

# =============================================================================
# FIGURE 3 — standardized coefficients of the selected models
# =============================================================================
dados_plot <- bind_rows(
  melhores_modelos %>%
    select(species, Coef = coef_static, P.value = P_static,
           IC_inf = CI_low_static, IC_sup = CI_high_static, radius_m, Ano = static_lag_k) %>%
    mutate(Tipo = "static"),
  melhores_modelos %>%
    select(species, Coef = coef_change, P.value = P_change,
           IC_inf = CI_low_change, IC_sup = CI_high_change, radius_m, Ano = change_interval_l) %>%
    mutate(Tipo = "delta")
) %>%
  mutate(Tipo = factor(Tipo, levels = c("static", "delta")))
dados_plot$species <- gsub("_", " ", dados_plot$species)

# points and bars whose CI crosses zero are drawn lighter
dados_plot$IC_cruza_zero <- dados_plot$IC_inf <= 0 & dados_plot$IC_sup >= 0
alpha_sig <- 1; alpha_nsig <- 0.45; alpha_static <- 0.5
dados_plot$alpha_ponto <- ifelse(dados_plot$IC_cruza_zero, alpha_nsig, alpha_sig)
dados_plot$alpha_barra <- dados_plot$alpha_ponto * ifelse(dados_plot$Tipo == "static", alpha_static, 1)

# colour = temporal scale of the predictor itself (static: lag k; change: interval l)
plot1 <- ggplot(dados_plot, aes(x = Coef, y = reorder(species, Coef), shape = Tipo, color = Ano)) +
  geom_errorbarh(aes(xmin = IC_inf, xmax = IC_sup, linewidth = Tipo, alpha = alpha_barra),
                 height = 0.4, position = position_dodge(width = 0.5)) +
  geom_point(aes(alpha = alpha_ponto), size = 2.5, position = position_dodge(width = 0.5)) +
  scale_alpha_identity() +
  scale_linewidth_manual(values = c(0.35, 0.7), name = "") +
  scale_color_viridis_c(name = "Temporal scale\n(years)", option = "magma", end = 0.85) +
  geom_vline(xintercept = 0, linetype = "solid", color = "red", linewidth = .5) +
  scale_shape_manual(values = c(1, 19), name = "") +
  labs(x = "Standardized coefficients", y = "") +
  scale_y_discrete(labels = function(x) parse(text = paste0("italic('", x, "')"))) +
  theme_classic() +
  theme(legend.position = c(0.85, 0.25),
        legend.justification = c(0.3, 0.4),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        axis.text.y = element_text(size = 8, face = "italic"),
        axis.text.x = element_text(size = 10)) +
  guides(color = guide_legend(order = 2),
         shape = guide_legend(order = 1),
         linewidth = guide_legend(order = 1))
print(plot1)
# "Figure 3" label in the top-left corner, in a strip above the plot (as in Figures 2 and 4)
fig3 <- cowplot::ggdraw() +
  cowplot::draw_plot(plot1, x = 0, y = 0, width = 1, height = 0.965) +
  cowplot::draw_label("Figure 3", x = 0.005, y = 0.995, hjust = 0, vjust = 1, fontface = "bold", size = 10)
ggsave(file.path(pasta_base, "Figure3.tif"), fig3, width = 17, height = 19, units = "cm",
       dpi = 600, compression = "lzw", bg = "white")

# =============================================================================
# PERMUTATION TEST OF THE WHOLE SELECTION PROCEDURE
# =============================================================================
# Shuffles the rows of the habitat predictors within each species (the same
# permutation for all lags, intervals and radii), keeping abundance and effort
# paired, and repeats exactly the same AIC-based selection.
n_perm <- 999   # value used in the manuscript (use 9 for a quick test)
set.seed(123)

melhor_modelo_especie <- function(especie, perm = NULL) {
  melhor.aic <- Inf; melhor.mod <- NULL
  for (ano.est in anos_estaticos) {
    for (ano.del in anos_delta) {
      for (escala in escalas_espaciais) {
        temp_df <- montar_dados(especie, ano.est, ano.del, escala, perm)
        if (nrow(temp_df) < 5) next
        if (!isTRUE(var(temp_df$preditor.current) > 0) || !isTRUE(var(temp_df$preditor.delta) > 0)) next
        temp_df <- padronizar_preditores(temp_df)
        mod <- tryCatch(ajustar_glm(temp_df), error = function(e) NULL)
        if (is.null(mod)) next
        if (AIC(mod) < melhor.aic) { melhor.aic <- AIC(mod); melhor.mod <- mod }
      }
    }
  }
  teste <- tryCatch(Anova(melhor.mod, type = "III"), error = function(e) NULL)
  if (is.null(teste)) return(c(p.curr = NA, p.delta = NA))
  c(p.curr  = teste["preditor.current", "Pr(>Chisq)"],
    p.delta = teste["preditor.delta", "Pr(>Chisq)"])
}

n.linhas <- sapply(especies, function(sp) sum(base$Species == sp))

# observed P-values of the selected model of each species
obs <- t(sapply(especies, function(sp) melhor_modelo_especie(sp, NULL)))
obs.sig <- obs[, "p.curr"] < 0.05 | obs[, "p.delta"] < 0.05
n.sig.obs <- sum(obs.sig, na.rm = TRUE)

# P-values of the selected model in each permutation
pcurr.perm <- matrix(NA, nrow = n_perm, ncol = length(especies), dimnames = list(NULL, especies))
pdelta.perm <- pcurr.perm
for (p in 1:n_perm) {
  for (sp in especies) {
    pv <- melhor_modelo_especie(sp, sample(n.linhas[sp]))
    pcurr.perm[p, sp]  <- pv["p.curr"]
    pdelta.perm[p, sp] <- pv["p.delta"]
  }
  print(paste("Permutation", p, "of", n_perm))
}
saveRDS(list(obs = obs, pcurr.perm = pcurr.perm, pdelta.perm = pdelta.perm),
        file.path(pasta_base, "permutation_p_values.rds"))

# assemblage level: number of species with P < 0.05 under the null
pmin.perm <- pmin(pcurr.perm, pdelta.perm)
sig.perm  <- pmin.perm < 0.05
taxa_falso_positivo <- colMeans(sig.perm, na.rm = TRUE)
n.sig.nulo <- rowSums(sig.perm, na.rm = TRUE)
p.global <- (sum(n.sig.nulo >= n.sig.obs) + 1) / (n_perm + 1)

# species level: permutation P-value = proportion of null P-values <= observed
p_permutacao <- function(p.obs, p.nulo) {
  if (is.na(p.obs)) return(NA)
  p.nulo <- p.nulo[!is.na(p.nulo)]
  (sum(p.nulo <= p.obs) + 1) / (length(p.nulo) + 1)
}
pmin.obs <- pmin(obs[, "p.curr"], obs[, "p.delta"])
P.perm.habitat <- sapply(especies, function(sp) p_permutacao(pmin.obs[sp], pmin.perm[, sp]))
P.perm.curr    <- sapply(especies, function(sp) p_permutacao(obs[sp, "p.curr"],  pcurr.perm[, sp]))
P.perm.delta   <- sapply(especies, function(sp) p_permutacao(obs[sp, "p.delta"], pdelta.perm[, sp]))

# "any" = either habitat predictor (smallest of the two P-values)
resul_perm <- data.frame(
  species = especies,
  P_static_obs = obs[, "p.curr"], P_change_obs = obs[, "p.delta"],
  any_P_below_0.05_obs = obs.sig,
  false_positive_rate = taxa_falso_positivo,
  P_perm_any = P.perm.habitat,   P_FDR_any = p.adjust(P.perm.habitat, method = "BH"),
  P_perm_static = P.perm.curr,   P_FDR_static = p.adjust(P.perm.curr, method = "BH"),
  P_perm_change = P.perm.delta,  P_FDR_change = p.adjust(P.perm.delta, method = "BH")
) %>%
  mutate(sig_any = P_FDR_any < 0.05, sig_static = P_FDR_static < 0.05, sig_change = P_FDR_change < 0.05) %>%
  left_join(melhores_modelos_todos %>%
              select(species, static_lag_k, change_interval_l, radius_m, coef_static, coef_change, P_dispersion),
            by = "species") %>%
  arrange(P_FDR_any)

cat("\nSpecies with an effect (nominal, P < 0.05):", n.sig.obs, "\n",
    "Mean expected by chance:", round(mean(n.sig.nulo), 2),
    "(", min(n.sig.nulo), "-", max(n.sig.nulo), ")\n",
    "Assemblage-level P:", round(p.global, 4), "\n",
    "After permutation + FDR (73 species):", sum(resul_perm$sig_any, na.rm = TRUE), "\n",
    "  static:", sum(resul_perm$sig_static, na.rm = TRUE),
    "(+", sum(resul_perm$sig_static & resul_perm$coef_static > 0, na.rm = TRUE),
    "/ -", sum(resul_perm$sig_static & resul_perm$coef_static < 0, na.rm = TRUE), ")\n",
    "  change:", sum(resul_perm$sig_change, na.rm = TRUE),
    "(+", sum(resul_perm$sig_change & resul_perm$coef_change > 0, na.rm = TRUE),
    "/ -", sum(resul_perm$sig_change & resul_perm$coef_change < 0, na.rm = TRUE), ")\n")

hist(n.sig.nulo, breaks = seq(-0.5, length(especies) + 0.5, 1),
     main = "Null distribution", xlab = "Number of species with a significant selected model",
     col = "grey80", border = "white")
abline(v = n.sig.obs, col = "red", lwd = 2)

# -----------------------------------------------------------------------------
# FDR among the 65 species with adequate residuals (values reported in the manuscript)
# -----------------------------------------------------------------------------
r65 <- subset(resul_perm, P_dispersion >= 0.05)
r65$P_FDR_any    <- p.adjust(r65$P_perm_any, "BH")
r65$P_FDR_static <- p.adjust(r65$P_perm_static, "BH")
r65$P_FDR_change <- p.adjust(r65$P_perm_change, "BH")
r65$sig_any    <- r65$P_FDR_any < 0.05
r65$sig_static <- r65$P_FDR_static < 0.05
r65$sig_change <- r65$P_FDR_change < 0.05
cat("\n65 species (adequate residuals): any habitat", sum(r65$sig_any),
    "| static", sum(r65$sig_static), "(+", sum(r65$sig_static & r65$coef_static > 0),
    "/ -", sum(r65$sig_static & r65$coef_static < 0), ")",
    "| change", sum(r65$sig_change), "(+", sum(r65$sig_change & r65$coef_change > 0),
    "/ -", sum(r65$sig_change & r65$coef_change < 0), ")\n")
# expected: any habitat 37 | static 7 (+1 / -6) | change 27 (+23 / -4)
write.xlsx(list(species_65 = r65, species_73 = resul_perm,
                null_distribution = data.frame(permutation = 1:n_perm, n_species_P_below_0.05 = n.sig.nulo)),
           file.path(pasta_base, "permutation_results.xlsx"))
