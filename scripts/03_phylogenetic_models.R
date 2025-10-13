{
library(lme4)
library(phyr)
library(ape)
library(phylobase)
library(adephylo)
library(yaml)
library(phytools)
library(caper)
library(glmmTMB)
library(cowplot)
library(AICcmodavg)
library(bbmle)
library(terra)
library(sf)
library(landscapemetrics)
library(openxlsx)
library(units)
library(dplyr)
library(stringr)
library(sjmisc)
library(ggplot2)
library(rlist)
library(nlme)
library(visreg)
library(car)
library(truncnorm)
library(spdep)
library(vegan)
library(fields)
library(MASS)
library(pscl)
library(rsq)
library(ggrepel)
library(DHARMa)
library(SuppDists)
library(multcomp)
library(emmeans)
}
setwd("C:\\Users\\User\\Desktop\\ArtigoCap1\\birdtree")
arvores = read.nexus("output.nex")
info = read_yaml("config.yaml")
info[["bad_names"]][["EricsonStage2Full_0004.tre"]]
tree = consensus.edges(arvores)
tree$node.label = 1:tree$Nnode
especies_remover <- c("Coereba_flaveola", "Pitangus_sulphuratus")
tree <- drop.tip(tree, especies_remover)

plotTree(tree)
nomes_antigos <- c("Basileuterus_flaveolus", "Basileuterus_leucoblepharus", "Buteo_magnirostris", 
                   "Parula_pitiayumi", "Thraupis_sayaca", "Troglodytes_aedon", "Vireo_olivaceus")
nomes_novos <- c("Myiothlypis_flaveola", "Myiothlypis_leucoblephara", "Rupornis_magnirostris", 
                 "Setophaga_pitiayumi", "Tangara_sayaca", "Troglodytes_musculus", "Vireo_chivi")

indices <- match(tree$tip.label, nomes_antigos)
tree$tip.label[!is.na(indices)] <- nomes_novos[indices[!is.na(indices)]]

setwd('C:/Users/User/Desktop/ArtigoCap1/GEE_planilhas')
melhores_modelos = read.xlsx("melhores_modelos.xlsx")

melhores_modelos = `rownames<-`(melhores_modelos, melhores_modelos$Especie)
melhores_modelos = melhores_modelos[match(tree$tip.label, melhores_modelos$Especie),]
correlacao_filogenetica = phylosig(tree, melhores_modelos$Coef.delta.hab, method = "lambda", test = TRUE)
phylosig(tree, melhores_modelos$Coef.delta.hab, method = "lambda", test = TRUE)
phylosig(tree, melhores_modelos$Coef.curr.hab, method = "lambda", test = TRUE)
plot(correlacao_filogenetica)

dados_long <- data.frame(
  coeficiente = c(melhores_modelos$Coef.delta.hab, melhores_modelos$Coef.curr.hab), 
  tipo_habitat = rep(c("Variação do Habitat", "Quantidade de Habitat"), each = nrow(melhores_modelos)),  # Repete os rótulos
  especies = rep(rownames(melhores_modelos), 2) 
)
dados_long$tipo_habitat <- as.factor(dados_long$tipo_habitat)
model = gls(coeficiente ~ tipo_habitat, correlation = corPagel(value = correlacao_filogenetica$lambda, phy = tree, form = ~especies), data = dados_long)
summary(model)
anova(model)
visreg(model)

visreg_output <- visreg(model)
visreg_data <- visreg_output$fit
visreg_res <- visreg_output$res

plot3 = ggplot() +
  geom_ribbon(data = visreg_data, aes(x = tipo_habitat, ymin = visregLwr, ymax = visregUpr), 
              alpha = 0.2, fill = "#4E84C4") +
  
  geom_point(data = visreg_res, aes(x = tipo_habitat, y = visregRes), 
             color = "#5A5A5A", size = 1.5, shape = 16, alpha = 0.4,
             position = position_jitter(width = 0.1, height = 0)) +
  
  geom_point(data = visreg_data, aes(x = tipo_habitat, y = visregFit), 
             color = "black", fill = "#D55E00", size = 3.5, shape = 21, stroke = 0.8) +

  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
  
  labs(x = "", y = "Coefficients") +
  
  scale_y_continuous(breaks = seq(-1.5, 1.5, by = 0.5)) +
  scale_x_discrete(labels = c("Quantidade de Habitat" = "Static habitat",
                              "Variação do Habitat" = "Delta habitat")) +
  theme_classic(base_size = 12) +
  theme(
    text = element_text(family = "sans", color = "black"),
    axis.text = element_text(color = "black", size = 10),
    axis.line = element_line(color = "black", linewidth = 0.4),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.2),
    plot.margin = unit(c(1, 1, 1, 1.5), "cm")
  )


plotx = grid.arrange(plot2, plot3, ncol=2)
library(gridExtra)
grid.arrange(plot1, plot3, ncol = 1)

ggsave("fig1.tif", width = 19, height = 27, units = "cm")
ggsave("fig2.tif", plot = plot2, width = 6, height = 9, units = "cm")
ggsave("fig3.tif", plot = plot3, width = 12, height = 10, units = "cm")

setwd("C:\\Users\\User\\Dropbox\\doutorado_cristiano\\dados\\AVONET")
avonet = read.xlsx("AVONET Supplementary dataset 1.xlsx", sheet = "AVONET2_eBird")

correcoes_taxonomicas = list(
  "Tangara sayaca" = "Thraupis sayaca",
  "Troglodytes musculus" = "Troglodytes aedon",
  "Veniliornis spilogaster" = "Dryobates spilogaster"
)

melhores_modelos = melhores_modelos %>%
  mutate(Especie = str_replace(Especie, "_", " "))

melhores_modelos <- melhores_modelos %>%
  mutate(Especie = dplyr::recode(Especie, !!!correcoes_taxonomicas))

melhores_modelos <- melhores_modelos %>%
  left_join(avonet %>% dplyr::select(Species2, Mass, Beak.Length_Culmen, Beak.Length_Nares,
                              Beak.Width, Beak.Depth, Tarsus.Length, Wing.Length,
                              Kipps.Distance, Hand_Wing.Index, Tail.Length,
                              Migration, Trophic.Level, Trophic.Niche, Primary.Lifestyle), by = c("Especie" = "Species2"))

melhores_modelos <- melhores_modelos %>%
  mutate(Especie_underline = gsub(" ", "_", Especie))

variaveis <- melhores_modelos %>%
  dplyr::select(Mass, Beak.Length_Culmen, Beak.Length_Nares,
         Beak.Width, Beak.Depth, Tarsus.Length, Wing.Length,
         Kipps.Distance, Hand_Wing.Index, Tail.Length,
         Migration)

pca <- prcomp(variaveis, scale. = TRUE)
summary(pca)

melhores_modelos$PC1 <- pca$x[,1]
melhores_modelos$PC2 <- pca$x[,2]
melhores_modelos$PC3 <- pca$x[,3]
melhores_modelos$PC4 <- pca$x[,4]
melhores_modelos$PC5 <- pca$x[,5]
melhores_modelos$PC6 <- pca$x[,6]

comparative_data <- comparative.data(tree, melhores_modelos, names.col = Especie_underline)

modelo_pgls_estatico <- pgls(Ano.estatico ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + Trophic.Level + Trophic.Niche + Primary.Lifestyle,
                    data = comparative_data,
                    lambda = "ML")

summary.pgls(modelo_pgls_estatico)

modelo_pgls_delta <- pgls(Ano_Delta ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + Trophic.Level + Trophic.Niche + Primary.Lifestyle,
                             data = comparative_data,
                             lambda = "ML")

summary.pgls(modelo_pgls_delta)

modelo_pgls_spatial <- pgls(Escala ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + Trophic.Level + Trophic.Niche + Primary.Lifestyle,
                          data = comparative_data,
                          lambda = "ML")

summary.pgls(modelo_pgls_spatial)
