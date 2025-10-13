library(readxl)
library(dplyr)
library(openxlsx)

setwd("C:\\Users\\User\\Desktop\\ArtigoCap1\\GEE_planilhas")
dados <- read_excel("data.xlsx")
forest_data <- read.csv("ForestPercentageTemporalSerie.csv")
forest_data = forest_data[,-ncol(forest_data)]
forest_data = forest_data[,-1]

buffer_sizes <- unique(forest_data$buffer_size)

for (buffer in buffer_sizes) {
  col_name <- paste0("buffer_size", buffer, "_delta_Year_finish2-year")
  dados[[col_name]] <- NA
}

dados$coord_id <- paste0(round(dados$Longitude_x, 5), "_", round(dados$Latitude_y, 5))
forest_data$coord_id <- paste0(round(forest_data$longitude, 5), "_", round(forest_data$latitude, 5))
forest_data <- forest_data %>%
  group_by(coord_id, year, buffer_size) %>%
  summarise(forest_percentage = mean(forest_percentage, na.rm = TRUE)) %>%
  ungroup()

for (buffer in buffer_sizes) {
  col_name <- paste0("buffer_size", buffer, "_delta_Year_finish2-year")
  
  dados[[col_name]] <- NA

  buffer_data <- forest_data %>%
    filter(buffer_size == buffer)

  buffer_data <- buffer_data %>%
    arrange(coord_id)
  dados <- dados %>%
    arrange(coord_id)

  for (i in 1:nrow(dados)) {
    coord_id <- dados$coord_id[i]
    year_finish <- dados$Year_finish2[i]
    year_previous <- year_finish - 10 # Change the number to vary the year delta

    current_data <- buffer_data %>%
      filter(year == year_finish)
    previous_data <- buffer_data %>%
      filter(year == year_previous)
    
    current_index <- which(current_data$coord_id == coord_id)
    previous_index <- which(previous_data$coord_id == coord_id)
    
    if (length(current_index) > 0 && length(previous_index) > 0) {

      current_forest <- current_data$forest_percentage[current_index]
      previous_forest <- previous_data$forest_percentage[previous_index]
 
      delta <- current_forest - previous_forest
      dados[[col_name]][i] <- delta
    }
  }
  print(data.frame(buf = buffer, ano_inicio = year_previous, ano_fim = year_finish,
                   cd.id = coord_id, lat = dados$Latitude_y[i], long = dados$Longitude_x[i],
                   ano_dados = dados$Year_finish2[i], f.atual = current_forest,
                   f.passado = previous_forest, perda = delta))
}  

openxlsx::write.xlsx(dados,'dados_delta0anos.xlsx')

