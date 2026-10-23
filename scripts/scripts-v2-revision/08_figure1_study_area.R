# =============================================================================
# 08 — FIGURE 1: study area and sampling sites
# =============================================================================
# Main map: Atlantic Forest (IBGE, 1:5,000,000) and sampling sites, coloured by
# habitat change over the longest interval analysed (l = 5 years) within the
# radius most often selected (20 km): dH = H(t0) - H(t0 - 5), in percentage
# points. Inset: location of the Atlantic Forest in Brazil.
# Changes from the previous version:
#  - axis labels use the standard notation (e.g. 55°W, 20°S) instead of "-55°W";
#  - colour shows 5-year change (the analyses use intervals of 1-5 years; the
#    old map showed 10-year change) on a diverging scale centred on zero
#    (red = loss, blue = gain), with a continuous colour bar and units;
#  - all points have the same size (size was mapped to the signed change, so
#    the largest losses looked like the smallest points);
#  - one point per site and sampling year; points are no longer cut at the map edge;
#  - all text >= 7 pt at final size; "Figure 1" label in bold.
# Output: Figure1.tif (600 dpi, 17 x 20 cm) and Figure1.png (preview).
# =============================================================================

library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(openxlsx)
library(stringr)
library(dplyr)
library(cowplot)

# Run from the repository root (paths are relative).
pasta_base  <- "results"
arquivo_mata <- file.path("data", "shapes", "shape_mata_atlantica_IBGE_5milhoes_policonica_sirgas2000.shp")
pasta_delta  <- file.path("data", "habitat_change")

intervalo <- 5        # change interval l (years)
raio      <- 20000    # buffer radius (m)

# ---- data -------------------------------------------------------------------------
mata_sf <- st_read(arquivo_mata, quiet = TRUE)
brasil  <- ne_countries(scale = "medium", country = "Brazil", returnclass = "sf")

# habitat change file for l = 5 years (same file-name pattern used in the analyses)
arq <- list.files(pasta_delta, pattern = "\\.xlsx$")
arq <- arq[str_detect(arq, "dados_delta") & as.numeric(str_extract(arq, "[0-9]+(?=anos)")) == intervalo]
stopifnot(length(arq) == 1)
dados <- read.xlsx(file.path(pasta_delta, arq))
coluna <- paste0("buffer_size", raio, "_delta_Year_finish2-year")

# one row per site and sampling year (the files have one row per species x site)
locais <- dados %>%
  transmute(Longitude_x, Latitude_y, Year_finish2, dH = .data[[coluna]]) %>%
  distinct(Longitude_x, Latitude_y, Year_finish2, .keep_all = TRUE) %>%
  arrange(abs(dH))                         # largest changes drawn on top
cat("Site-years:", nrow(locais), "| unique sites:",
    nrow(distinct(locais, Longitude_x, Latitude_y)),
    "| without 5-year change (NA):", sum(is.na(locais$dH)), "\n")
cat("Range of change (percentage points):", round(range(locais$dH, na.rm = TRUE), 2), "\n")
pontos <- st_as_sf(locais, coords = c("Longitude_x", "Latitude_y"), crs = 4326)

lim_cor <- max(abs(locais$dH), na.rm = TRUE)   # symmetric colour scale around zero

# ---- inset: Brazil -------------------------------------------------------------------
mapa_brasil <- ggplot() +
  geom_sf(data = brasil, fill = "grey90", colour = "grey60", linewidth = 0.2) +
  geom_sf(data = mata_sf, fill = "#A6CEE3", colour = NA) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.5))

# ---- main map ------------------------------------------------------------------------
mapa_mata <- ggplot() +
  geom_sf(data = mata_sf, fill = "#A6CEE3", colour = NA, alpha = 0.7) +
  geom_sf(data = pontos, aes(fill = dH), shape = 21, size = 1.8, colour = "grey15", stroke = 0.25) +
  scale_fill_gradient2(low = "#B2182B", mid = "#F7F7F7", high = "#2166AC", midpoint = 0,
                       limits = c(-lim_cor, lim_cor), na.value = "white",
                       name = paste0("Habitat change over ", intervalo, " years, \u0394H(", intervalo,
                                     ")\n(percentage points, ", raio / 1000, "-km radius)")) +
  annotation_north_arrow(location = "bl", which_north = "true",
                         style = north_arrow_fancy_orienteering(text_size = 7),
                         height = unit(0.9, "cm"), width = unit(0.9, "cm"),
                         pad_x = unit(0.4, "cm"), pad_y = unit(0.4, "cm")) +
  annotation_scale(location = "br", width_hint = 0.25, text_cex = 0.6,
                   pad_x = unit(0.4, "cm"), pad_y = unit(0.4, "cm")) +
  coord_sf(expand = TRUE) +                 # default graticule labels: 55°W, 20°S
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 8) +
  theme(axis.text = element_text(size = 7, colour = "grey20"),
        panel.grid = element_line(colour = "grey90", linewidth = 0.25),
        legend.position = "bottom",
        legend.title = element_text(size = 7.5, hjust = 0.5),
        legend.text = element_text(size = 7),
        legend.key.width = unit(1.6, "cm"), legend.key.height = unit(0.25, "cm")) +
  guides(fill = guide_colourbar(title.position = "top", frame.colour = "grey40", ticks.colour = "grey40"))

# ---- combine: inset + figure label ----------------------------------------------------
final_plot <- ggdraw() +
  draw_plot(mapa_mata, x = 0, y = 0, width = 1, height = 0.97) +
  draw_plot(mapa_brasil, x = 0.17, y = 0.70, width = 0.24, height = 0.24) +
  draw_label("Figure 1", x = 0.01, y = 0.995, hjust = 0, vjust = 1, fontface = "bold", size = 10)

ggsave(file.path(pasta_base, "Figure1.tif"), final_plot, width = 17, height = 20, units = "cm",
       dpi = 600, compression = "lzw", bg = "white")
ggsave(file.path(pasta_base, "Figure1.png"), final_plot, width = 17, height = 20, units = "cm", dpi = 200, bg = "white")
print(final_plot)
