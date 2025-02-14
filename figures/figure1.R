##### Figure 1 Mapping invasion times 
# all cities with epidemics excluding San Andres and Providencia

## CHIKV

# Load these packages
library(raster)
library(sp)
library(ggplot2)   # for general plotting
library(ggpubr)    # for arranging combined ggplots
library(dplyr)
library(rnaturalearth)

# Reading shapefile data into R
spdf_world <- ne_countries()

# crop the extent
e <- extent(-80,-66,-6,12)

america <- crop(spdf_world, e)

plot(america)

admin0 <- fortify(america)
admin1 <- filter(admin0, group == 22.1)
admin2 <- filter(admin0, group == 35.1) # Colombia
admin3 <- filter(admin0, group == 46.1)
admin4 <- filter(admin0, group == 123.1)
admin5 <- filter(admin0, group == 124.1)
admin6 <- filter(admin0, group == 170.1)

# base map

map <- ggplot() +
  geom_polygon(data = admin1, aes(x = long, y = lat, group = group),
               color = "black", fill = 'gray', size = 0.25) +
  geom_polygon(data = admin2, aes(x = long, y = lat, group = group),
               color = "black", fill = 'black', size = 0.25) + 
  geom_polygon(data = admin3, aes(x = long, y = lat, group = group),
               color = "black", fill = 'gray', size = 0.25) +
  geom_polygon(data = admin4, aes(x = long, y = lat, group = group),
               color = "black", fill = 'gray', size = 0.25) +
  geom_polygon(data = admin5, aes(x = long, y = lat, group = group),
               color = "black", fill = 'gray', size = 0.25) +
  geom_polygon(data = admin6, aes(x = long, y = lat, group = group),
               color = "black", fill = 'gray', size = 0.25) +
  coord_fixed(xlim =c(-79, -67), ylim=c(-5, 13), ratio = 1) +
  xlab ('') + ylab('') +
  #remove coordinate lines on water
  theme(panel.background = element_rect(fill = "lightblue"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  #remove margins
  theme(plot.margin=unit(c(0,0,0,0),"mm"))  

map

############# DATA

# Import lat/long of cities
lat_lon <- read.csv("figures/data/lat_lon_complete.txt", 
                    sep='\t', 
                    stringsAsFactors = FALSE,
                    header = TRUE)

# format and pad admin2 with zeros
lat_lon$admin2 <- sprintf("%05d", lat_lon$admin2)
lat_lon$admin2 <- as.character(lat_lon$admin2)

# load disease data
chik <- read.csv("figures/data/chik_first_reported_cases_338.txt", # invasion times start at week 1
                 sep='\t', 
                 stringsAsFactors = FALSE,
                 header = TRUE)

# format and pad admin2 with zeros
chik$admin2_code <- sprintf("%05d", chik$admin2_code)
chik$admin2_code <- as.character(chik$admin2_code) 

chik <- rename(chik, admin2 = admin2_code)

# Merge latitudes and longitudes with city infection week by admin 2 
pop_map2 <- merge(lat_lon, chik, by = "admin2", all.y = TRUE)

# Sort by pop to reveal small circles
pop_map2 <- pop_map2[with(pop_map2, order( -1 * pop)), ] 

########## PUTTING EVERYTHING TOGETHER (panels are each 12 weeks)

data_1 <- filter(pop_map2, first_report <13)
data_2 <- filter(pop_map2, first_report >= 13 & first_report < 25)
data_3 <- filter(pop_map2, first_report >= 25 & first_report < 37)
data_4 <- filter(pop_map2, first_report >= 37 & first_report < 49)
data_5 <- filter(pop_map2, first_report >= 49 & first_report < 61)
data_6 <- filter(pop_map2, first_report >= 61)

# fake data to fix legend key labels
data_7 <- data_6
data_7$admin2_longitude <- 55
data_7$admin2_latitude <- 55

# cities with earliest infection onsets
map1 <- map +
  geom_point(data = data_1, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi)), color = "green") +
  geom_blank(data = data_2, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_blank(data = data_3, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_blank(data = data_4, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_blank(data = data_5, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_blank(data = data_6, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_point(data = data_7, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi)), color = "black") +
  scale_size(range = c(1,8), breaks = c(200, 400, 600, 800), 
             labels = c("population = 100,000", "population = 500,000", 
                        "population = 1,000,000", "population = 2,000,000")) +
  theme(plot.margin=unit(c(0,0,0,0),"mm")) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 16)) +
  theme(legend.key.size = unit(1.5, 'lines'))

#map1

map2 <- map +
  geom_blank(data = data_1, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_point(data = data_2, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi)), color = "yellow") +
  geom_blank(data = data_3, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_blank(data = data_4, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_blank(data = data_5, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_blank(data = data_6, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  scale_size(range = c(1,8), breaks = c(200, 400, 600, 800), 
             labels = c("population = 100,000", "population = 500,000", 
                        "population = 1,000,000", "population = 2,000,000")) +
  theme(plot.margin=unit(c(0,0,0,0),"mm")) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 16)) +
  theme(legend.key.size = unit(1.5, 'lines'))

#map2

map3 <- map +
  geom_blank(data = data_1, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_blank(data = data_2, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_point(data = data_3, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi)), color = "red") +
  geom_blank(data = data_4, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_blank(data = data_5, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_blank(data = data_6, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  scale_size(range = c(1,8), breaks = c(200, 400, 600, 800), 
             labels = c("population = 100,000", "population = 500,000", 
                        "population = 1,000,000", "population = 2,000,000")) +
  theme(plot.margin=unit(c(0,0,0,0),"mm")) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 16)) +
  theme(legend.key.size = unit(1.5, 'lines'))

#map3

map4 <- map +
  geom_blank(data = data_1, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_blank(data = data_2, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_blank(data = data_3, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_point(data = data_4, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi)), color = "magenta") +
  geom_blank(data = data_5, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_blank(data = data_6, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  scale_size(range = c(1,8), breaks = c(200, 400, 600, 800), 
             labels = c("population = 100,000", "population = 500,000", 
                        "population = 1,000,000", "population = 2,000,000")) +
  theme(plot.margin=unit(c(0,0,0,0),"mm")) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 16)) +
  theme(legend.key.size = unit(1.5, 'lines'))

#map4

map5 <- map +
  geom_blank(data = data_1, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_blank(data = data_2, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_blank(data = data_3, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_blank(data = data_4, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_point(data = data_5, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi)), color = "blue") +
  geom_blank(data = data_6, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  scale_size(range = c(1,8), breaks = c(200, 400, 600, 800), 
             labels = c("population = 100,000", "population = 500,000", 
                        "population = 1,000,000", "population = 2,000,000")) +
  theme(plot.margin=unit(c(0,0,0,0),"mm")) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 16)) +
  theme(legend.key.size = unit(1.5, 'lines'))

#map5

map6 <- map +
  geom_blank(data = data_1, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_blank(data = data_2, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_blank(data = data_3, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_blank(data = data_4, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_blank(data = data_5, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_point(data = data_6, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi)), color = "cyan") +
  scale_size(range = c(1,8), breaks = c(200, 400, 600, 800), 
             labels = c("population = 100,000", "population = 500,000", 
                        "population = 1,000,000", "population = 2,000,000")) +
  theme(plot.margin=unit(c(0,0,0,0),"mm")) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 16)) +
  theme(legend.key.size = unit(1.5, 'lines'))


#map6

chik_map <- ggarrange(map1, map2, map3, map4, map5, map6, ncol = 3, nrow=2, legend="none", 
                      labels = c("Weeks\n   0-11", "Weeks\n  12-23", "Weeks\n 24-35", "Weeks\n 36-47", "Weeks\n 48-59", "Weeks\n 60-71"),
                      # labels = c("A", "B", "C", "D", "E","F"),
                      font.label = list(size = 14))

saveRDS(chik_map, "chik_map_for_figure1_first_reports.RDS")

###############################################

## ZIKV

# Load these packages
library(raster)
library(sp)
library(ggplot2)   # for general plotting
library(ggpubr)    # for arranging combined ggplots
library(dplyr)
library(rnaturalearth)

# Reading shapefile data into R
spdf_world <- ne_countries()

# crop the extent
e <- extent(-80,-66,-6,12)

america <- crop(spdf_world, e)

plot(america)

admin0 <- fortify(america)
admin1 <- filter(admin0, group == 22.1)
admin2 <- filter(admin0, group == 35.1) # Colombia
admin3 <- filter(admin0, group == 46.1)
admin4 <- filter(admin0, group == 123.1)
admin5 <- filter(admin0, group == 124.1)
admin6 <- filter(admin0, group == 170.1)

# base map

map <- ggplot() +
  geom_polygon(data = admin1, aes(x = long, y = lat, group = group),
               color = "black", fill = 'gray', size = 0.25) +
  geom_polygon(data = admin2, aes(x = long, y = lat, group = group),
               color = "black", fill = 'black', size = 0.25) + 
  geom_polygon(data = admin3, aes(x = long, y = lat, group = group),
               color = "black", fill = 'gray', size = 0.25) +
  geom_polygon(data = admin4, aes(x = long, y = lat, group = group),
               color = "black", fill = 'gray', size = 0.25) +
  geom_polygon(data = admin5, aes(x = long, y = lat, group = group),
               color = "black", fill = 'gray', size = 0.25) +
  geom_polygon(data = admin6, aes(x = long, y = lat, group = group),
               color = "black", fill = 'gray', size = 0.25) +
  coord_fixed(xlim =c(-79, -67), ylim=c(-5, 13), ratio = 1) +
  xlab ('') + ylab('') +
  #remove coordinate lines on water
  theme(panel.background = element_rect(fill = "lightblue"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  #remove margins
  theme(plot.margin=unit(c(0,0,0,0),"mm"))  

map


############# DATA

# Import lat/long of cities
lat_lon <- read.csv("figures/data/lat_lon_complete.txt", 
                    sep='\t', 
                    stringsAsFactors = FALSE,
                    header = TRUE)

# format and pad admin2 with zeros
lat_lon$admin2 <- sprintf("%05d", lat_lon$admin2)
lat_lon$admin2 <- as.character(lat_lon$admin2)

# load disease data
zika <- read.csv("figures/data/zika_first_reported_cases_288.txt", # weeks start at 1
                 sep='\t', 
                 stringsAsFactors = FALSE,
                 header = TRUE)

# format and pad admin2 with zeros
zika$admin2_code <- sprintf("%05d", zika$admin2_code)
zika$admin2_code <- as.character(zika$admin2_code) 

zika <- rename(zika, admin2 = admin2_code)

# Merge latitudes and longitudes with city infection week by admin 2 
pop_map2 <- merge(lat_lon, zika, by = "admin2", all.y = TRUE)

# Sort by pop to reveal small circles
pop_map2 <- pop_map2[with(pop_map2, order( -1 * pop)), ] 

########## PUTTING EVERYTHING TOGETHER (each panel is 6 weeks)

# break up cities by invasion week
data_1 <- filter(pop_map2, first_report <7)
data_2 <- filter(pop_map2, first_report >= 7 & first_report < 13)
data_3 <- filter(pop_map2, first_report >= 13 & first_report < 19)
data_4 <- filter(pop_map2, first_report >= 19 & first_report < 25)
data_5 <- filter(pop_map2, first_report >= 25 & first_report < 31)
data_6 <- filter(pop_map2, first_report >= 31)

# fake data to fix legend key labels
data_7 <- data_6
data_7$admin2_longitude <- 55
data_7$admin2_latitude <- 55

map1 <- map +
  geom_point(data = data_1, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi)), color = "green") +
  geom_blank(data = data_2, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_blank(data = data_3, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_blank(data = data_4, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_blank(data = data_5, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_blank(data = data_6, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_point(data = data_7, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi)), color = "black") +
  scale_size(range = c(1,8), breaks = c(200, 400, 600, 800), 
             labels = c("population = 100,000", "population = 500,000", 
                        "population = 1,000,000", "population = 2,000,000")) +
  theme(plot.margin=unit(c(0,0,0,0),"mm")) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 16)) +
  theme(legend.key.size = unit(1.5, 'lines'))

#map1

map2 <- map +
  geom_blank(data = data_1, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_point(data = data_2, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi)), color = "yellow") +
  geom_blank(data = data_3, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_blank(data = data_4, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_blank(data = data_5, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_blank(data = data_6, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  scale_size(range = c(1,8), breaks = c(200, 400, 600, 800), 
             labels = c("population = 100,000", "population = 500,000", 
                        "population = 1,000,000", "population = 2,000,000")) +
  theme(plot.margin=unit(c(0,0,0,0),"mm")) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 16)) +
  theme(legend.key.size = unit(1.5, 'lines'))

#map2

map3 <- map +
  geom_blank(data = data_1, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_blank(data = data_2, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_point(data = data_3, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi)), color = "red") +
  geom_blank(data = data_4, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_blank(data = data_5, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_blank(data = data_6, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  scale_size(range = c(1,8), breaks = c(200, 400, 600, 800), 
             labels = c("population = 100,000", "population = 500,000", 
                        "population = 1,000,000", "population = 2,000,000")) +
  theme(plot.margin=unit(c(0,0,0,0),"mm")) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 16)) +
  theme(legend.key.size = unit(1.5, 'lines'))

#map3

map4 <- map +
  geom_blank(data = data_1, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_blank(data = data_2, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_blank(data = data_3, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_point(data = data_4, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi)), color = "magenta") +
  geom_blank(data = data_5, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_blank(data = data_6, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  scale_size(range = c(1,8), breaks = c(200, 400, 600, 800), 
             labels = c("population = 100,000", "population = 500,000", 
                        "population = 1,000,000", "population = 2,000,000")) +
  theme(plot.margin=unit(c(0,0,0,0),"mm")) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 16)) +
  theme(legend.key.size = unit(1.5, 'lines'))

#map4

map5 <- map +
  geom_blank(data = data_1, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_blank(data = data_2, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_blank(data = data_3, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_blank(data = data_4, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_point(data = data_5, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi)), color = "blue") +
  geom_blank(data = data_6, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  scale_size(range = c(1,8), breaks = c(200, 400, 600, 800), 
             labels = c("population = 100,000", "population = 500,000", 
                        "population = 1,000,000", "population = 2,000,000")) +
  theme(plot.margin=unit(c(0,0,0,0),"mm")) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 16)) +
  theme(legend.key.size = unit(1.5, 'lines'))

#map5

map6 <- map +
  geom_blank(data = data_1, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_blank(data = data_2, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_blank(data = data_3, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_blank(data = data_4, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_blank(data = data_5, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi))) +
  geom_point(data = data_6, aes(x = admin2_longitude, y = admin2_latitude, 
                                size = sqrt(pop/pi)), color = "cyan") +
  scale_size(range = c(1,8), breaks = c(200, 400, 600, 800), 
             labels = c("population = 100,000", "population = 500,000", 
                        "population = 1,000,000", "population = 2,000,000")) +
  theme(plot.margin=unit(c(0,0,0,0),"mm")) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 16)) +
  theme(legend.key.size = unit(1.5, 'lines'))


#map6

zika_map <- ggarrange(map1, map2, map3, map4, map5, map6, ncol = 3, nrow=2, common.legend = TRUE, legend="bottom",
                      labels = c("Weeks\n   0-5", "Weeks\n  6-11", "Weeks\n 12-17", "Weeks\n 18-23", "Weeks\n 24-29", "Weeks\n 30-35"),
                      # labels = c("A", "B", "C", "D", "E","F"),
                      font.label = list(size = 14))

saveRDS(zika_map, "zika_map_for_figure1_first_reports.RDS")

##############################################

# put CHIKV map and ZIKV maps together

library(ggplot2)   # for general plotting
library(ggpubr)    # for arranging combined ggplots

chik_map <- readRDS("chik_map_for_figure1_first_reports.RDS")
zika_map <- readRDS("zika_map_for_figure1_first_reports.RDS")

library("cowplot")
both <- plot_grid(chik_map, zika_map, 
                  labels = c("A", "B"),
                  label_size = 20,
                  ncol = 1, nrow = 2)

both # save as jpeg or png, width 1100, height 1400


