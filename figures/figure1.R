##### Figure 1 Mapping invasion times 
# all cities with epidemics excluding San Andres and Providencia

# updated to fix CHIKV dates reported week 1 of 2014

## CHIKV

# Load these packages
library(raster)
library(sp)
library(rgdal)     # R wrapper around GDAL/OGR
library(ggplot2)   # for general plotting
library(ggmap)    # for fortifying shapefiles
library(ggpubr)    # for arranging combined ggplots
library(dplyr)

# Reading areal polygon data into R

par(mar=c(0,0,0,0))

# Colombia shapefiles
adm1data <- getData('GADM', country = 'COL', level = 1)
adm0data <- getData('GADM', country = 'COL', level = 0)

# Bordering countries
adm0data_VEN <- getData('GADM', country = 'VEN', level = 0)
adm0data_BRA <- getData('GADM', country = 'BRA', level = 0)
adm0data_PER <- getData('GADM', country = 'PER', level = 0)
adm0data_ECU <- getData('GADM', country = 'ECU', level = 0)
adm0data_PAN <- getData('GADM', country = 'PAN', level = 0)

# Prepare shapefile for use in ggplot
COL_admin0 <- fortify(adm0data)
COL_admin1 <- fortify(adm1data)

VEN <- fortify(adm0data_VEN)
BRA <- fortify(adm0data_BRA)
PER <- fortify(adm0data_PER)
ECU <- fortify(adm0data_ECU)
PAN <- fortify(adm0data_PAN)

# base map

map <- ggplot() +
  geom_polygon(data = COL_admin0, aes(x = long, y = lat, group = group),
               color = "black", fill = 'gray', size = 0.25) + 
  geom_polygon(data = VEN, aes(x = long, y = lat, group = group),
               color = "black", fill = 'gray', size = 0.25) + 
  geom_polygon(data = BRA, aes(x = long, y = lat, group = group),
               color = "black", fill = 'gray', size = 0.25) + 
  geom_polygon(data = PER, aes(x = long, y = lat, group = group),
               color = "black", fill = 'gray', size = 0.25) + 
  geom_polygon(data = ECU, aes(x = long, y = lat, group = group),
               color = "black", fill = 'gray', size = 0.25) + 
  geom_polygon(data = PAN, aes(x = long, y = lat, group = group),
               color = "black", fill = 'gray', size = 0.25) + 
  
  
  geom_polygon(data = COL_admin1, aes(x = long, y = lat, group = group),
               color = "white", fill = 'black', size = 0.25) +
  theme(panel.background = element_rect(fill = "black",
                                        colour = "white",
                                        size = 0.5, linetype = "solid")) +
  coord_fixed(xlim =c(-79, -67), ylim=c(-5, 13), ratio = 1) +
  xlab ('') + ylab('') +
  #remove coordinate lines on water
  theme(panel.background = element_rect(fill = "lightblue"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  #remove margins
  theme(plot.margin=unit(c(0,0,0,0),"mm"))  

map

############# DATA

# Import lat/long of cities
lat_lon <- readRDS("figures/data/lat_lon_complete.RDS")

# load disease data
chik <- readRDS('figures/data/chik_first_reported_cases_335.RDS') # invasion times start at week 1
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

map1

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

map2

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

map3

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

map4

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

map5

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


map6

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
library(rgdal)     # R wrapper around GDAL/OGR
library(ggplot2)   # for general plotting
library(ggmap)    # for fortifying shapefiles
library(ggpubr)    # for arranging combined ggplots
library(dplyr)

# Reading areal polygon data into R

par(mar=c(0,0,0,0))

# Colombia shapefiles
adm1data <- getData('GADM', country = 'COL', level = 1)
adm0data <- getData('GADM', country = 'COL', level = 0)

# Bordering countries
adm0data_VEN <- getData('GADM', country = 'VEN', level = 0)
adm0data_BRA <- getData('GADM', country = 'BRA', level = 0)
adm0data_PER <- getData('GADM', country = 'PER', level = 0)
adm0data_ECU <- getData('GADM', country = 'ECU', level = 0)
adm0data_PAN <- getData('GADM', country = 'PAN', level = 0)

# Prepare shapefile for use in ggplot
COL_admin0 <- fortify(adm0data)
COL_admin1 <- fortify(adm1data)

VEN <- fortify(adm0data_VEN)
BRA <- fortify(adm0data_BRA)
PER <- fortify(adm0data_PER)
ECU <- fortify(adm0data_ECU)
PAN <- fortify(adm0data_PAN)

# base map

map <- ggplot() +
  geom_polygon(data = COL_admin0, aes(x = long, y = lat, group = group),
               color = "black", fill = 'gray', size = 0.25) + 
  geom_polygon(data = VEN, aes(x = long, y = lat, group = group),
               color = "black", fill = 'gray', size = 0.25) + 
  geom_polygon(data = BRA, aes(x = long, y = lat, group = group),
               color = "black", fill = 'gray', size = 0.25) + 
  geom_polygon(data = PER, aes(x = long, y = lat, group = group),
               color = "black", fill = 'gray', size = 0.25) + 
  geom_polygon(data = ECU, aes(x = long, y = lat, group = group),
               color = "black", fill = 'gray', size = 0.25) + 
  geom_polygon(data = PAN, aes(x = long, y = lat, group = group),
               color = "black", fill = 'gray', size = 0.25) + 
  
  
  geom_polygon(data = COL_admin1, aes(x = long, y = lat, group = group),
               color = "white", fill = 'black', size = 0.25) +
  theme(panel.background = element_rect(fill = "black",
                                        colour = "white",
                                        size = 0.5, linetype = "solid")) +
  coord_fixed(xlim =c(-79, -67), ylim=c(-5, 13), ratio = 1) +
  xlab ('') + ylab('') +
  #remove coordinate lines on water
  theme(panel.background = element_rect(fill = "lightblue"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  #remove margins
  theme(plot.margin=unit(c(0,0,0,0),"mm"))  

map

############# DATA

# Import lat/long of cities
lat_lon <- readRDS("figures/data/lat_lon_complete.RDS")

# load disease data
zika <- readRDS('figures/data/zika_first_reported_cases_288.RDS') # weeks start at 1
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

map1

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

map2

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

map3

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

map4

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

map5

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


map6

zika_map <- ggarrange(map1, map2, map3, map4, map5, map6, ncol = 3, nrow=2, common.legend = TRUE, legend="bottom",
                      labels = c("Weeks\n   0-5", "Weeks\n  6-11", "Weeks\n 12-17", "Weeks\n 18-23", "Weeks\n 24-29", "Weeks\n 30-35"),
                      # labels = c("A", "B", "C", "D", "E","F"),
                      font.label = list(size = 14))

saveRDS(zika_map, "zika_map_for_figure1_first_reports.RDS")

##############################################

# put CHIKV map and ZIKV maps together

library(ggplot2)   # for general plotting
library(ggmap)    # for fortifying shapefiles
library(ggpubr)    # for arranging combined ggplots

chik_map <- readRDS("chik_map_for_figure1_first_reports.RDS")
zika_map <- readRDS("zika_map_for_figure1_first_reports.RDS")

library("cowplot")
both <- plot_grid(chik_map, zika_map, 
                  labels = c("A", "B"),
                  label_size = 20,
                  ncol = 1, nrow = 2)

both # save as jpeg or png, width 1100, height 1400


