library(openxlsx)
library(MASS)
library(stringr)
library(ggplot2)
library(dplyr)
library(DHARMa)
library(car)
library(visreg)
library(tidyr)
library(cowplot)

setwd("C:\\Users\\User\\Desktop\\ArtigoCap1\\GEE_planilhas\\valores_abs_porAno")
tabelas.current <- list.files()[1:11]
setwd("C:\\Users\\User\\Desktop\\ArtigoCap1\\GEE_planilhas\\valores_delta")
tabelas.delta <- list.files()[1:11]

lista_AICs <- list()

contador <- 1

for (i in 1:length(tabelas.current)) {  
  setwd("C:\\Users\\User\\Desktop\\ArtigoCap1\\GEE_planilhas\\valores_abs_porAno")
  tabela.current.i <- read.xlsx(tabelas.current[i])  
  
  especies <- unique(tabela.current.i$Species) 
  
  for (k in 1:length(tabelas.delta)) { 
    setwd("C:\\Users\\User\\Desktop\\ArtigoCap1\\GEE_planilhas\\valores_delta")
    tabela.delta.k <- read.xlsx(tabelas.delta[k])
    
    for (especie in especies) {  
      tabela.current <- filter(tabela.current.i, Species == especie)
      tabela.delta <- filter(tabela.delta.k, Species == especie)
      
      for (j in 10:ncol(tabela.current)) { 
        temp_df <- data.frame(
          PE_PAI = tabela.current$PE_PAI,
          preditor.current = tabela.current[, j],
          preditor.delta = tabela.delta[, j],  
          PE_total_effort = tabela.current$PE_total_effort
        )
        temp_df <- na.omit(temp_df)
        temp_df$PE_total_effort[temp_df$PE_total_effort == 0] <- 0.001
        
       
        if (var(temp_df$preditor.current, na.rm = TRUE) == 0) next
        
        mod.i <- tryCatch(
          glm(PE_PAI ~ preditor.current + preditor.delta + PE_total_effort, 
              family = Gamma(link = "log"), data = temp_df),
          error = function(e) {
            message(paste("Erro ao ajustar modelo para:", tabelas.current[i], 
                          "especie:", especie, "coluna:", colnames(tabela.current)[j], 
                          "delta:", tabelas.delta[k]))
            return(NULL)
          }
        )
        
        if (is.null(mod.i)) next
        
        resul <- tryCatch(summary(mod.i), error = function(e) return(NULL))
        teste <- tryCatch(Anova(mod.i, type = 'III'), error = function(e) return(NULL))
        
        if (is.null(resul) || is.null(teste)) next
        
      
        lista_AICs[[contador]] <- data.frame(
          Ano.estatico = as.numeric(str_extract(tabelas.current[i], "[0-9]+(?=anos)")),
          Ano.delta = as.numeric(str_extract(tabelas.delta[i], "[0-9]+(?=anos)")),
          Escala = as.numeric(str_extract(str_extract(paste0(tabelas.current[i], colnames(tabela.current)[j]), "buffer_size[0-9]+"), "[0-9]+")),
          Ano_Delta = as.numeric(str_extract(tabelas.delta[k], "[0-9]+(?=anos)")),
          AIC = AIC(mod.i),
          Especie = especie,
          Coef.delta.hab = resul$coefficients['preditor.delta', 'Estimate'],
          Coef.curr.hab = resul$coefficients['preditor.current', 'Estimate'],
          P.value.delta.hab = teste$`Pr(>Chisq)`[2],
          P.value.curr.hab = teste$`Pr(>Chisq)`[1],
          pseudoR2 = (resul$null.deviance-resul$deviance)/resul$null.deviance
        )
        
        contador <- contador + 1
        print(c(i, k, j))
      }
    }
  }
}



df_AICsABS = do.call(rbind, lista_AICs)
df_AICsABS = filter(df_AICsABS, P.value.delta.hab < 0.05 | P.value.curr.hab < 0.05)


head(df_AICsABS)
getwd()
write.xlsx(df_AICsABS, "resul_coefs_allCombinations.xlsx")

setwd("C:\\Users\\User\\Desktop\\ArtigoCap1\\GEE_planilhas")
df_AICsABS = read.xlsx('resul_coefs_allCombinations.xlsx')
df_AICsABS$Anos = paste0(df_AICsABS$Ano.estatico, ".", df_AICsABS$Ano_Delta)


for(especie in unique(df_AICsABS$Especie)) {
  df_especie = df_AICsABS[df_AICsABS$Especie == especie, ]
  menor_aic_info = df_especie %>%
    arrange(AIC) %>%
    slice(1) %>%
    summarise(ano_menor_aic = Anos, escala_menor_aic = Escala)
  
  ano_menor_aic = menor_aic_info$ano_menor_aic
  escala_menor_aic = menor_aic_info$escala_menor_aic
  
  p = ggplot(df_especie, aes(x = as.numeric(Escala), y = as.numeric(AIC), color = as.numeric(Anos), group = Anos)) +
    geom_line(size = .6) +  
    theme_classic() +
    labs(
      title = paste(especie, "- Ano =", ano_menor_aic, ", Escala =", escala_menor_aic),
      x = "Escala",
      y = "AIC",
      color = "Ano"
    ) +
    scale_color_gradient(low = "black", high = "lightgray", name = "Ano")
  
  print(p)
  
  readline(prompt = "Pressione [Enter] para ver a próxima espécie...")
}

melhores_modelos <- df_AICsABS %>%
  group_by(Especie) %>%
  slice_min(order_by = AIC, n = 1) %>%
  ungroup()

lista_modelos <- list()
modelos.resid.ok = data.frame(Especie = NA, Ano.Est = NA, Ano.Del = NA, Escala = NA, P.resid = NA)

for (i in 1:nrow(melhores_modelos)) {
  especie <- melhores_modelos$Especie[i]
  ano.est <- melhores_modelos$Ano.estatico[i]
  ano.del <- melhores_modelos$Ano_Delta[i]
  escala <- melhores_modelos$Escala[i]
  

  setwd("C:\\Users\\User\\Desktop\\ArtigoCap1\\GEE_planilhas\\valores_abs_porAno")
  tabela_nome.current <- tabelas.current[str_detect(tabelas.current, paste0("dados_absolutos_", ano.est, "anos_arrumado"))]
  tabela.current <- read.xlsx(tabela_nome.current) %>% filter(Species == especie)
  setwd("C:\\Users\\User\\Desktop\\ArtigoCap1\\GEE_planilhas\\valores_delta")
  tabela_nome.delta <- tabelas.delta[str_detect(tabelas.delta, paste0("dados_delta", ano.del, "anos"))]
  tabela.delta <- read.xlsx(tabela_nome.delta) %>% filter(Species == especie)
  
  tabela.current <- tabela.current %>% arrange(PE_PAI)
  tabela.delta <- tabela.delta %>% arrange(PE_PAI)
  

  temp_df <- data.frame(
    PE_PAI = tabela.current$PE_PAI,
    preditor.current = tabela.current[, which(colnames(tabela.current) == paste0("buffer_size", escala, "_absolute_Year_finish2_minus_", ano.est, "years"))],
    preditor.delta = tabela.delta[, which(colnames(tabela.delta) == paste0("buffer_size", escala, "_delta_Year_finish2-year"))],
    PE_total_effort = tabela.current$PE_total_effort
  )
  temp_df <- na.omit(temp_df)
  temp_df$PE_total_effort[temp_df$PE_total_effort == 0] <- 0.001
  
  mod <- glm(PE_PAI ~ preditor.current + preditor.delta + PE_total_effort, family = Gamma(link = "log"), data = temp_df)

  sim_res <- simulateResiduals(fittedModel = mod, n = 1000)

  layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE), heights = c(.8, 1))
  par(mar=c(4, 4, 3, 2) + 0.2)
  visreg(mod, 'preditor.current', main = c(especie, str_c("Ano", " ", ano.est), str_c("Escala", " ", escala)))
  visreg(mod, 'preditor.delta', main = c(especie, str_c("Ano", " ", ano.del), str_c("Escala", " ", escala)))
  testDispersion(sim_res)

  layout(1)
  

  lista_modelos[[especie]] <- list(
    modelo = mod,
    residuos = sim_res
  )
  readline(prompt = ' ')
}

nova_tabela <- melhores_modelos

for (i in 1:nrow(melhores_modelos)) {
  especie <- melhores_modelos$Especie[i]
  ano.est <- melhores_modelos$Ano.estatico[i]
  ano.del <- melhores_modelos$Ano_Delta[i]
  escala <- melhores_modelos$Escala[i]
  
  setwd("C:\\Users\\User\\Desktop\\ArtigoCap1\\GEE_planilhas\\valores_abs_porAno")
  tabela_nome.current <- tabelas.current[str_detect(tabelas.current, paste0("dados_absolutos_", ano.est))]
  tabela.current <- read.xlsx(tabela_nome.current) %>% filter(Species == especie)
  
  setwd("C:\\Users\\User\\Desktop\\ArtigoCap1\\GEE_planilhas\\valores_delta")
  tabela_nome.delta <- tabelas.delta[str_detect(tabelas.delta, paste0("dados_delta", ano.del, "anos"))]
  tabela.delta <- read.xlsx(tabela_nome.delta) %>% filter(Species == especie)
  
  tabela.current <- tabela.current %>% arrange(PE_PAI)
  tabela.delta <- tabela.delta %>% arrange(PE_PAI)

  temp_df <- data.frame(
    PE_PAI = tabela.current$PE_PAI,
    preditor.current = tabela.current[, which(colnames(tabela.current) == paste0("buffer_size", escala, "_absolute_Year_finish2_minus_", ano.est, "years"))],
    preditor.delta = tabela.delta[, which(colnames(tabela.delta) == paste0("buffer_size", escala, "_delta_Year_finish2-year"))],
    PE_total_effort = tabela.current$PE_total_effort
  )
  temp_df <- na.omit(temp_df)
  temp_df$PE_total_effort[temp_df$PE_total_effort == 0] <- 0.001
  
  mod <- glm(PE_PAI ~ preditor.current + preditor.delta + PE_total_effort, family = Gamma(link = "log"), data = temp_df)

  sim_res <- simulateResiduals(fittedModel = mod, n = 1000)
  p_valor_disp <- testDispersion(sim_res)$p.value
 
  if (p_valor_disp > 0.05) {
    mod1 <- glm(PE_PAI ~ preditor.current + preditor.delta + PE_total_effort, family = Gamma(link = "log"), data = temp_df)
  }
  
  coefs <- coef(summary(mod1))
  conf_int <- Confint(mod1, exponentiate = F, type = "LR")
  
  p_value_col <- grep("Pr\\(>\\|[zt]\\|\\)", colnames(coefs), value = TRUE)
  
  nova_tabela$Coef.delta.hab[i] <- conf_int["preditor.delta", "Estimate"]
  nova_tabela$Coef.curr.hab[i] <- conf_int["preditor.current", "Estimate"]
  nova_tabela$P.value.delta.hab[i] <- coefs["preditor.delta", p_value_col]
  nova_tabela$P.value.curr.hab[i] <- coefs["preditor.current", p_value_col]
  nova_tabela$IC_inf.delta[i] <- conf_int["preditor.delta", 2]
  nova_tabela$IC_sup.delta[i] <- conf_int["preditor.delta", 3]
  nova_tabela$IC_inf.curr[i] <- conf_int["preditor.current", 2]
  nova_tabela$IC_sup.curr[i] <- conf_int["preditor.current", 3]
  nova_tabela$P.resid[i] = testDispersion(simulateResiduals(fittedModel = mod1, n = 1000))$p.value
}

melhores_modelos = nova_tabela
melhores_modelos = filter(melhores_modelos, P.resid >= 0.05)
write.xlsx(melhores_modelos, "melhores_modelos.xlsx")

setwd('C:/Users/User/Desktop/ArtigoCap1/GEE_planilhas')

melhores_modelos = read.xlsx("melhores_modelos.xlsx")

dados_plot <- bind_rows(
  melhores_modelos %>%
    select(Especie,
           Coef = Coef.curr.hab,
           P.value = P.value.curr.hab,
           IC_inf = IC_inf.curr,
           IC_sup = IC_sup.curr,
           Escala,
           Anos) %>%
    mutate(Tipo = "static"),
  
  melhores_modelos %>%
    select(Especie,
           Coef = Coef.delta.hab,
           P.value = P.value.delta.hab,
           IC_inf = IC_inf.delta,
           IC_sup = IC_sup.delta,
           Escala,
           Anos) %>%
    mutate(Tipo = "delta")
) %>%
  mutate(
    Tipo = factor(Tipo, levels = c("static", "delta")),
    Significativo = ifelse(P.value < 0.05, "Significativo (p < 0.05)", "Não significativo")
  )

{
ggplot(dados_plot, aes(x = Coef, y = Especie, shape = Tipo, color = as.numeric(Anos))) +
  geom_errorbarh(aes(xmin = IC_inf, xmax = IC_sup), 
                 height = 0.5, position = position_dodge(width = 0.5)) +
  geom_point(aes(size=as.numeric(Escala)), position = position_dodge(width = 0.5)) +
  #scale_color_manual(values = c("#E69F00", "#56B4E9")) +
  scale_color_viridis_c(name = "Riqueza", option = "plasma", na.value = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  scale_shape_manual(values = c(1, 19)) +
  labs(x = "Coeficientes") +
  theme_classic() +
  theme(legend.position = "top",
        legend.box = "horizontal",
        legend.title = element_blank(),
        legend.text = element_text(size = 9),
        axis.text.y = element_text(size = 9))
  }

dados_plot$Especie <- gsub("_", " ", dados_plot$Especie)


plot1 = ggplot(dados_plot, aes(x = Coef, y = reorder(Especie, Coef), shape = Tipo, color = as.numeric(Anos))) +
  geom_errorbarh(aes(xmin = IC_inf, xmax = IC_sup), 
                 height = 0.4, position = position_dodge(width = 0.5), linewidth = 0.6) +
  geom_point(aes(size = as.numeric(Escala)), position = position_dodge(width = 0.5)) +
  scale_color_viridis_c(name = "Temporal scale", option = "magma") +
  geom_vline(xintercept = 0, linetype = "solid", color = "red", linewidth = .5) +
  annotate("rect", xmin = -0.2, xmax = 0.2, ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = "gray10")+
  scale_shape_manual(values = c(1, 19), name = "") +
  scale_size_continuous(name = "Spatial scale") +
  labs(x = "Coefficients", y = "") +
  scale_y_discrete(labels = function(x) parse(text = paste0("italic('", x, "')"))) +
  theme_classic() +
  theme(legend.position = c(0.85, 0.25),
        legend.justification = c(0.3, 0.4),
        legend.direction = "vertical",
        legend.box = "vertical",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        axis.text.y = element_text(size = 8, face = "italic"),
        axis.text.x = element_text(size = 10)) +
  guides(
    color = guide_legend(order = 2),
    size = guide_legend(order = 3),
    shape = guide_legend(order = 1))

ggsave('fig1.tif', plot1, width = 17, height = 19, units = 'cm')

n_spp_positive_delta = filter(melhores_modelos, Coef.delta.hab > 0)
n_spp_positive_delta = filter(n_spp_positive_delta, P.value.delta.hab < 0.05)
length(n_spp_positive_delta$Especie)

n_spp_negative_delta = filter(melhores_modelos, Coef.delta.hab < 0)
n_spp_negative_delta = filter(n_spp_negative_delta, P.value.delta.hab < 0.05)
length(n_spp_negative_delta$Especie)

n_spp_positive_curr = filter(melhores_modelos, Coef.curr.hab > 0)
n_spp_positive_curr = filter(n_spp_positive_curr, P.value.curr.hab < 0.05)
length(n_spp_positive_curr$Especie)

n_spp_negative_delta = filter(melhores_modelos, Coef.delta.hab < 0)
n_spp_negative_delta = filter(n_spp_negative_delta, P.value.delta.hab < 0.05)
length(n_spp_negative_delta$Especie)

n_spp_sem_efeito = filter(melhores_modelos, P.value.delta.hab > 0.05 & P.value.curr.hab > 0.05)
length(n_spp_sem_efeito$Especie)

nrow(melhores_modelos[which(melhores_modelos$Escala == min(melhores_modelos$Escala)),])/72
mean(melhores_modelos$Escala)
sd(melhores_modelos$Escala)

nrow(melhores_modelos[which(melhores_modelos$Ano_Delta == min(melhores_modelos$Ano_Delta)),])/72
summary(melhores_modelos$Ano_Delta)
mean(melhores_modelos$Ano_Delta)
sd(melhores_modelos$Ano_Delta)

nrow(melhores_modelos[which(melhores_modelos$Ano.estatico == min(melhores_modelos$Ano.estatico)),])/72
summary(melhores_modelos$Ano.estatico)
mean(melhores_modelos$Ano.estatico)
sd(melhores_modelos$Ano.estatico)

library(ggalluvial)
dados <- data.frame(
  categoria = rep(c("∆hab", "Habitat atual"), each = 72),
  efeito = c(ifelse(1:72 %in% c(1:36), "Yes", "No"),
             ifelse(1:72 %in% c(1:2, 15:28), "Yes", "No")),
  especie = rep(1:72, 2)
)

plot2 = ggplot(dados, aes(x = categoria, stratum = efeito, alluvium = especie, fill = efeito)) +
  geom_flow() +
  geom_stratum() +
  labs(fill = "Variable effect") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "top"
  )

plot2 = ggplot(dados_plot, aes(x = Tipo, y = Coef)) +
  geom_jitter(width = 0.1, color = "grey70", size = 1.5, alpha = 0.7) +
  
  stat_summary(fun = mean, geom = "point", color = "orange", size = 4) +
  #stat_summary(fun.data = mean_cl_boot, geom = "errorbar",
  #             width = 0.1, color = "orange", size = 1) +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "", y = "Coefficients") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"))

final_plot = ggdraw() +
  draw_plot(plot1) +
  draw_plot(plot2, x = 0.65, y = 0.60, width = 0.3, height = 0.3)

ggsave('fig1.tif', final_plot, width = 17, height = 19, units = 'cm')
