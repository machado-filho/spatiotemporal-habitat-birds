# =============================================================================
# Figure 2 — analytical workflow (schematic)
# =============================================================================
# Four panels read left-to-right, top-to-bottom:
#   a  spatial scales: seven buffer radii around a sampling site
#   b  temporal predictors: static habitat amount H(t0 - k) and habitat change
#      dH(l) = H(t0) - H(t0 - l)
#   c  model selection: 175 candidate GLMs per species, lowest AIC selected
#   d  inference: the whole selection procedure repeated on 999 permutations
# Panels a and c are schematic (simulated landscape and AIC surface).
# Panel d is also schematic (the real null distribution is a result and is
# shown in Appendix S4, Figure S1).
# Output: results/Figure2.tif (600 dpi, 18 x 13 cm) and results/Figure2.png
# All text is >= 7 pt at final size (ggplot text size in mm x 2.845 = pt).
# Requires only ggplot2 and cowplot.
# =============================================================================

library(ggplot2)
library(cowplot)

pasta_saida <- "results"
dir.create(pasta_saida, showWarnings = FALSE)

# ---- visual style -------------------------------------------------------------
tinta   <- "grey15"
cinza   <- "grey55"
floresta <- "#8DB38B"
matriz  <- "grey95"
destaque <- "#B5452F"
tam_txt <- 2.6          # main annotation text size (mm; ~7.4 pt)
tam_menor <- 2.5        # smaller labels (mm; ~7.1 pt)
mostrar_equacao <- FALSE  # TRUE = model equation centred below panel c (it is already in the Methods)

tema_base <- theme_void(base_size = 8) +
  theme(plot.title = element_text(size = 8.5, face = "bold", hjust = 0, margin = margin(b = 3)),
        plot.subtitle = element_text(size = 7, colour = "grey30", hjust = 0, margin = margin(b = 4)),
        plot.margin = margin(4, 6, 4, 6))

circulo <- function(r, n = 400) {
  a <- seq(0, 2 * pi, length.out = n)
  data.frame(x = r * cos(a), y = r * sin(a))
}

# ---- (a) spatial scales --------------------------------------------------------
set.seed(7)
n <- 180
paisagem <- expand.grid(x = seq(-1.05, 1.05, length.out = n), y = seq(-1.05, 1.05, length.out = n))
campo <- with(paisagem,
              sin(3.1 * x + 0.4) * cos(2.7 * y - 0.3) +
              0.6 * sin(5.3 * x * y + 1) +
              0.5 * cos(4.1 * x - 3.2 * y) +
              0.3 * sin(9 * x + 2) * cos(8 * y))
paisagem$floresta <- campo > 0.1

raios_m <- c(200, 600, 1000, 2000, 5000, 10000, 20000)
# radii drawn on a log scale so that all seven are visible (not to scale)
raios_plot <- 0.1 + (log10(raios_m) - log10(200)) / (log10(20000) - log10(200)) * 0.88
circs <- do.call(rbind, lapply(seq_along(raios_plot), function(i)
  cbind(circulo(raios_plot[i]), id = i)))
# labels alternate between the upper-right and lower-right diagonals so that
# neighbouring radii do not overlap at 7 pt
angulo_rot <- ifelse(seq_along(raios_plot) %% 2 == 1, pi / 4, -pi / 4)
rotulos_raio <- data.frame(
  x = raios_plot * cos(angulo_rot), y = raios_plot * sin(angulo_rot),
  lab = c("0.2", "0.6", "1", "2", "5", "10", "20 km"))

pa <- ggplot() +
  geom_raster(data = paisagem, aes(x, y, fill = floresta)) +
  scale_fill_manual(values = c(`TRUE` = floresta, `FALSE` = matriz), guide = "none") +
  geom_path(data = circs, aes(x, y, group = id), colour = tinta, linewidth = 0.3) +
  geom_label(data = rotulos_raio, aes(x, y, label = lab), size = tam_menor, label.size = 0,
             label.padding = unit(0.06, "lines"), fill = "white", colour = tinta) +
  annotate("point", x = 0, y = 0, shape = 21, size = 1.8, fill = "white", colour = tinta, stroke = 0.6) +
  coord_equal(xlim = c(-1.05, 1.05), ylim = c(-1.05, 1.05), expand = FALSE) +
  labs(title = "a   Spatial scales",
       subtitle = "Native forest cover (H) within seven buffer radii (s)") +
  tema_base

# ---- (b) temporal predictors ---------------------------------------------------
anos <- -5:0
rotulos_anos <- parse(text = c("t[0]-5", "t[0]-4", "t[0]-3", "t[0]-2", "t[0]-1", "t[0]"))
estaticos <- data.frame(x = -c(0, 2, 3, 4, 5), y = 1)
mudancas <- data.frame(l = 1:5, x0 = -(1:5), x1 = 0, y = -0.15 - 0.22 * (0:4))

pb <- ggplot() +
  # time axis
  annotate("segment", x = -5.3, xend = 0.3, y = 0.55, yend = 0.55, colour = tinta, linewidth = 0.35) +
  annotate("segment", x = anos, xend = anos, y = 0.52, yend = 0.58, colour = tinta, linewidth = 0.35) +
  annotate("text", x = anos, y = 0.42, label = as.character(rotulos_anos), parse = TRUE,
           size = tam_txt, colour = tinta) +
  annotate("segment", x = 0, xend = 0, y = -1.2, yend = 1.25, colour = cinza,
           linewidth = 0.3, linetype = "22") +
  annotate("text", x = 0.12, y = 1.25, label = "sampling year", hjust = 0, size = tam_menor, colour = cinza) +
  # static habitat amount
  geom_point(data = estaticos, aes(x, y), shape = 21, size = 2.2, fill = floresta, colour = tinta, stroke = 0.4) +
  annotate("text", x = -5.3, y = 1.25, hjust = 0, size = tam_txt, colour = tinta, parse = TRUE,
           label = '"Static habitat amount  "*H(t[0]-k)*",  "*k %in% "{0, 2, 3, 4, 5}"') +
  # habitat change
  geom_segment(data = mudancas, aes(x = x0, xend = x1, y = y, yend = y),
               colour = destaque, linewidth = 0.45) +
  geom_segment(data = mudancas, aes(x = x0, xend = x0, y = y - 0.05, yend = y + 0.05),
               colour = destaque, linewidth = 0.45) +
  geom_segment(data = mudancas, aes(x = x1, xend = x1, y = y - 0.05, yend = y + 0.05),
               colour = destaque, linewidth = 0.45) +
  geom_text(data = mudancas, aes(x = x0 - 0.12, y = y, label = paste0("l = ", l)),
            hjust = 1, size = tam_menor, colour = destaque) +
  annotate("text", x = -5.3, y = 0.08, hjust = 0, size = tam_txt, colour = tinta, parse = TRUE,
           label = '"Habitat change  "*Delta*H(l)==H(t[0])-H(t[0]-l)') +
  coord_cartesian(xlim = c(-5.9, 1.6), ylim = c(-1.25, 1.35), clip = "off") +
  labs(title = "b   Temporal predictors",
       subtitle = "Measured at each radius, up to five years before sampling") +
  tema_base

# ---- (c) model selection ---------------------------------------------------------
set.seed(3)
grade <- expand.grid(raio = 1:7, comb = 1:25)
grade$dAIC <- with(grade, 0.9 * (raio - 6)^2 + 0.35 * abs(comb - 8) + runif(nrow(grade), 0, 1.6))
grade$dAIC <- grade$dAIC - min(grade$dAIC)
sel <- grade[which.min(grade$dAIC), ]

eq <- 'log*" E("*A[i]*")"=="  "*beta[0]+beta[1]*" z["*H[s*","*i](t[0]-k)*"]"+beta[2]*" z["*Delta*H[s*","*i](l)*"]"+beta[3]*E[i]'

pc <- ggplot(grade, aes(comb, raio, fill = pmin(dAIC, 12))) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_tile(data = sel, fill = NA, colour = destaque, linewidth = 0.7) +
  scale_fill_gradient(low = "grey20", high = "grey95", guide = "none") +
  scale_y_continuous(breaks = 1:7, labels = c("0.2", "0.6", "1", "2", "5", "10", "20"), expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = "c   Model selection (per species)",
       subtitle = "175 gamma GLMs: 7 radii x 25 combinations of k and l; lowest AIC selected",
       x = "combination of static lag (k) and change interval (l)", y = "radius (km)",
       caption = if (mostrar_equacao) parse(text = eq) else NULL) +
  tema_base +
  theme(axis.text.y = element_text(size = 7, colour = tinta, margin = margin(r = 2)),
        axis.title.y = element_text(size = 7, angle = 90, colour = tinta, margin = margin(r = 3)),
        axis.title.x = element_text(size = 7, colour = tinta, margin = margin(t = 3)),
        plot.caption = element_text(size = 7.5, colour = tinta, hjust = 0.5, margin = margin(t = 6)),
        plot.caption.position = "panel")

# ---- (d) inference (schematic: no results shown) ------------------------------
# Illustrates the procedure only; the real null distribution is a result and is
# shown in Appendix S4 (Figure S1).
x_nulo <- seq(-4, 4, length.out = 400)
curva <- data.frame(x = x_nulo, y = dnorm(x_nulo))
pd <- ggplot(curva, aes(x, y)) +
  geom_area(fill = "grey82", colour = NA) +
  geom_line(colour = "grey45", linewidth = 0.35) +
  annotate("segment", x = 5.2, xend = 5.2, y = 0, yend = 0.36, colour = destaque, linewidth = 0.6) +
  annotate("text", x = 5.0, y = 0.39, hjust = 1, size = tam_menor, colour = destaque, label = "observed") +
  annotate("text", x = 0, y = 0.45, size = tam_menor, colour = "grey30",
           label = "null distribution\n(999 permutations)") +
  annotate("segment", x = -4.6, xend = 6.2, y = 0, yend = 0, colour = tinta, linewidth = 0.3) +
  scale_x_continuous(limits = c(-4.6, 6.4), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.52), expand = c(0, 0)) +
  labs(title = "d   Inference",
       subtitle = "Whole procedure repeated on permuted data;\npermutation P-values, FDR across species",
       x = "number of species with a\nsignificant habitat predictor", y = NULL) +
  tema_base +
  theme(axis.title.x = element_text(size = 7, colour = tinta, margin = margin(t = 3), lineheight = 0.95))

# ---- assemble --------------------------------------------------------------------
linha1 <- plot_grid(pa, pb, nrow = 1, rel_widths = c(0.42, 0.58))
linha2 <- plot_grid(pc, pd, nrow = 1, rel_widths = c(0.62, 0.38))
figura <- ggdraw() +
  draw_plot(plot_grid(linha1, linha2, ncol = 1, rel_heights = c(1, 0.95)),
            x = 0, y = 0, width = 1, height = 0.955) +
  draw_label("Figure 2", x = 0.005, y = 0.995, hjust = 0, vjust = 1, fontface = "bold", size = 10)

ggsave(file.path(pasta_saida, "Figure2.tif"), figura, width = 18, height = 13, units = "cm",
       dpi = 600, compression = "lzw", bg = "white")
ggsave(file.path(pasta_saida, "Figure2.png"), figura, width = 18, height = 13, units = "cm",
       dpi = 200, bg = "white")
print(figura)
