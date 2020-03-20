# hexagon map animation for CHIKV

# load libraries
library(tidyverse)
library(viridis)
library(hrbrthemes)
library(mapdata)
library(raster)
library(sp)
library(geojsonio)
library(RColorBrewer)
library(rgdal)
library(ggplot2)
library(tigris)
library(geogrid)
library(sf)
library(tmap)
library(gganimate)
library(transformr)

# load disease data
chik <- readRDS("animations/data/chik_monthly_inc_for_animation_incl2016.RDS")

new_df <- chik

# shapefile for Colombia
adm1data <- getData('GADM', country = 'COL', level = 1)

# admin 1
new_cells_hex <- calculate_grid(shape = adm1data, grid_type = "hexagonal", seed = 39)
resulthex <- assign_polygons(adm1data, new_cells_hex)

hexplot <- tm_shape(resulthex) + 
  tm_polygons() +
  tm_text("NAME_1")

#hexplot

centers <- cbind.data.frame(coordinates(resulthex), as.character(resulthex$NAME_1))
names(centers) <- c("x", "y", "id")

# fortify for use in ggplot
chik_map  <- fortify(resulthex, region = "NAME_1")

# San Andres does not match spatial polygons
new_df$admin1_name[new_df$admin1_name == "San Andrés" ] <- "San Andrés y Providencia"
new_df$id <- new_df$admin1_name

# fake data to avoid ugly gray in log transformation of color scale on plot
new_df_for_log <- new_df
new_df_for_log$inc_rate_x100 <- new_df_for_log$inc_rate * 100
new_df_for_log$inc_rate_x100[new_df_for_log$inc_rate_x100  <1] <- 1

hex_map_data_month <- as.tbl(chik_map) %>%
  full_join(new_df_for_log) %>%
  left_join(centers) %>%
  filter(!is.na(long))

# animation
p <- ggplot(hex_map_data_month, aes(long, lat)) +
  geom_polygon(aes(group = group, fill = inc_rate_x100),
               color = "white") +
  geom_text(aes(label = dep_short, x = x, y = y), color = "white", size = 4) + # transition_states(hex_map_data_month$month_no, transition_length = 1, states = hex_map_data_month$month_no) 
  transition_time(hex_map_data_month$month_no) + # states is name of column holding state levels in data
  labs(title = "Chikungunya incidence per 100,000 in 
       {hex_map_data_month$month_name[as.integer(frame_time)]} {hex_map_data_month$year[as.integer(frame_time)]}",
       subtitle = "Colombia") +
  scale_fill_viridis_c(trans = "log",
                       breaks = c(1, 30, 450, 6000, 97000, 98000),
                       labels = c(0, 0.3, 45, 60, 970, 980)) + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(2, "cm"),
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 20),
        plot.subtitle = element_text(size = 14))

chik_hexagon_map <- animate(p, nframes = 100, fps = 10, height = 600, width = 400)  

anim_save("chik_hexagon_map.gif", chik_hexagon_map)

