# =============================================================================
# 03 — PHYLOGENETIC SIGNAL AND PGLS (traits x scales of effect)
# =============================================================================
# Run from the repository root (paths are relative):
#   data/phylogeny/  output.nex                           BirdTree sample of trees
#   data/traits/     AVONET Supplementary dataset 1.xlsx  Tobias et al. 2022
#   results/         selected_models.xlsx (script 01) and outputs
#
#  - Uses the selected models of the 65 species with adequate residuals
#    (sheet "residuals_ok").
#  - The consensus tree is pruned automatically to the species in the data.
#  - Phylogenetic signal (Pagel's lambda) of the standardized coefficients of
#    static habitat and habitat change.
#  - PGLS of the three scales of effect (static lag k, change interval l and
#    radius) against species traits. Radius enters on a log10 scale because
#    the radii are geometrically spaced (200-20,000 m).
#
# Outputs (results/):
#   pgls_results.txt    full printed output
#   pgls_summary.xlsx   pgls_summary, pgls_coefficients, phylogenetic_signal,
#                       pca_importance
# =============================================================================

library(ape)
library(phytools)
library(caper)
library(openxlsx)
library(dplyr)
library(stringr)

# ---- PATHS --------------------------------------------------------------------
pasta_base   <- "results"
pasta_arvore <- file.path("data", "phylogeny")   # output.nex (BirdTree)
arquivo_avonet <- file.path("data", "traits", "AVONET Supplementary dataset 1.xlsx")  # Tobias et al. 2022
arquivo_melhores <- "selected_models.xlsx"
aba_melhores <- "residuals_ok"

# ---- TREE ---------------------------------------------------------------------
arvores <- read.nexus(file.path(pasta_arvore, "output.nex"))
tree <- consensus.edges(arvores)
tree$node.label <- 1:tree$Nnode

# BirdTree names -> names used in the data (current taxonomy)
nomes_antigos <- c("Basileuterus_flaveolus", "Basileuterus_leucoblepharus", "Buteo_magnirostris",
                   "Parula_pitiayumi", "Thraupis_sayaca", "Troglodytes_aedon", "Vireo_olivaceus")
nomes_novos   <- c("Myiothlypis_flaveola", "Myiothlypis_leucoblephara", "Rupornis_magnirostris",
                   "Setophaga_pitiayumi", "Tangara_sayaca", "Troglodytes_musculus", "Vireo_chivi")
indices <- match(tree$tip.label, nomes_antigos)
tree$tip.label[!is.na(indices)] <- nomes_novos[indices[!is.na(indices)]]

# ---- DATA ---------------------------------------------------------------------
melhores <- read.xlsx(file.path(pasta_base, arquivo_melhores), sheet = aba_melhores)

# match tree and data
fora_da_arvore <- setdiff(melhores$species, tree$tip.label)
cat("Species in the data that are NOT in the tree (left out):", if (length(fora_da_arvore)) fora_da_arvore else "none", "\n")
tree <- drop.tip(tree, setdiff(tree$tip.label, melhores$species))
melhores <- melhores[match(tree$tip.label, melhores$species), ]
rownames(melhores) <- melhores$species
cat("Species in the analysis:", nrow(melhores), "\n")

# ---- PHYLOGENETIC SIGNAL OF THE COEFFICIENTS -----------------------------------
set.seed(123)
sinal_delta    <- phylosig(tree, setNames(melhores$coef_change, melhores$species), method = "lambda", test = TRUE)
sinal_estatico <- phylosig(tree, setNames(melhores$coef_static, melhores$species), method = "lambda", test = TRUE)
print(sinal_delta); print(sinal_estatico)

# ---- TRAITS (AVONET, eBird taxonomy) --------------------------------------------
avonet <- read.xlsx(arquivo_avonet, sheet = "AVONET2_eBird")
correcoes_taxonomicas <- c("Tangara sayaca" = "Thraupis sayaca",
                           "Troglodytes musculus" = "Troglodytes aedon",
                           "Veniliornis spilogaster" = "Dryobates spilogaster")

dados <- melhores %>%
  mutate(species_underscore = species,
         species_avonet = str_replace(species, "_", " "),
         species_avonet = dplyr::recode(species_avonet, !!!correcoes_taxonomicas)) %>%
  left_join(avonet %>% dplyr::select(Species2, Mass, Beak.Length_Culmen, Beak.Length_Nares,
                                     Beak.Width, Beak.Depth, Tarsus.Length, Wing.Length,
                                     Kipps.Distance, Hand_Wing.Index, Tail.Length,
                                     Migration, Trophic.Level, Trophic.Niche, Primary.Lifestyle),
            by = c("species_avonet" = "Species2"))

sem_trait <- dados$species[is.na(dados$Mass)]
cat("Species without AVONET traits (check the name):", if (length(sem_trait)) sem_trait else "none", "\n")
dados <- dados[!is.na(dados$Mass), ]

# PCA of the continuous traits (+ migration)
variaveis <- dados %>% dplyr::select(Mass, Beak.Length_Culmen, Beak.Length_Nares, Beak.Width, Beak.Depth,
                                     Tarsus.Length, Wing.Length, Kipps.Distance, Hand_Wing.Index,
                                     Tail.Length, Migration)
pca <- prcomp(variaveis, scale. = TRUE)
print(summary(pca))
for (k in 1:6) dados[[paste0("PC", k)]] <- pca$x[, k]

# radius on a log10 scale
dados$log10_radius <- log10(dados$radius_m)

comparative_data <- comparative.data(tree, dados, names.col = species_underscore)

# ---- PGLS -----------------------------------------------------------------------
formula_traits <- "~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + Trophic.Level + Trophic.Niche + Primary.Lifestyle"
respostas <- c(static = "static_lag_k", change = "change_interval_l", spatial = "log10_radius")

modelos <- list(); resumo <- list(); coefs <- list()
for (nm in names(respostas)) {
  f <- as.formula(paste(respostas[nm], formula_traits))
  m <- pgls(f, data = comparative_data, lambda = "ML")
  s <- summary(m)
  fs <- s$fstatistic
  modelos[[nm]] <- m
  resumo[[nm]] <- data.frame(
    response = respostas[nm],
    n = nrow(comparative_data$data),
    lambda = m$param["lambda"],
    R2 = s$r.squared, R2_adj = s$adj.r.squared,
    F = unname(fs["value"]), df1 = unname(fs["numdf"]), df2 = unname(fs["dendf"]),
    P_model = unname(pf(fs["value"], fs["numdf"], fs["dendf"], lower.tail = FALSE)),
    n_terms_P_below_0.05 = sum(s$coefficients[-1, 4] < 0.05)
  )
  co <- as.data.frame(s$coefficients); co$term <- rownames(co); co$response <- respostas[nm]
  coefs[[nm]] <- co
  print(s)
}
resumo_pgls <- do.call(rbind, resumo)
coefs_pgls  <- do.call(rbind, coefs)
print(resumo_pgls)

# ---- SAVE -------------------------------------------------------------------------
sink(file.path(pasta_base, "pgls_results.txt"))
cat("Species not in the tree:", fora_da_arvore, "\n")
cat("Species without traits:", sem_trait, "\n\n")
cat("== Phylogenetic signal: habitat change coefficient ==\n"); print(sinal_delta)
cat("\n== Phylogenetic signal: static habitat coefficient ==\n"); print(sinal_estatico)
cat("\n== PCA ==\n"); print(summary(pca))
for (nm in names(modelos)) { cat("\n== PGLS:", respostas[nm], "==\n"); print(summary(modelos[[nm]])) }
cat("\n== Summary ==\n"); print(resumo_pgls)
sink()

importancia <- as.data.frame(summary(pca)$importance)
importancia <- cbind(statistic = rownames(importancia), importancia)
write.xlsx(list(
  pgls_summary = resumo_pgls,
  pgls_coefficients = coefs_pgls,
  phylogenetic_signal = data.frame(
    coefficient = c("change", "static"),
    lambda = c(sinal_delta$lambda, sinal_estatico$lambda),
    logL = c(sinal_delta$logL, sinal_estatico$logL),
    P = c(sinal_delta$P, sinal_estatico$P)),
  pca_importance = importancia
), file.path(pasta_base, "pgls_summary.xlsx"))
