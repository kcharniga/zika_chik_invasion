##### Correlation plots

library(dplyr)

## CHIKV

# Import lat/long of cities
lat_lon <- readRDS("figures/data/lat_lon_complete.RDS")

# city invasion times
invasion_times <- readRDS('figures/data/chik_first_reported_cases_338.RDS') 
invasion_times <- rename(invasion_times, admin2 = admin2_code)

# subtract 1 from first reports to start at 0
invasion_times$first_report <- invasion_times$first_report - 1

# Merge latitudes and longitudes with invasion times by admin 2 
invasion_plot <- merge(lat_lon, invasion_times, by = "admin2", all.y = TRUE)

# import matrix of distance between cities
d_cities <- readRDS("figures/data/distance_matrix_km.RDS")

# Re-order rows and columns in distance matrix to match rows in invasion times df
row.order <- rownames(invasion_times)
d_cities <- d_cities[row.order, row.order]

# change NAs on diagonal to 0
d_cities[is.na(d_cities)] <- 0

# extract row of distances that correspond to most likely origin city (Barranquilla)  - fourth row
first_city_dist <- d_cities[4,]

# create named vector of invasion times
first_report_ordered <- as.numeric(invasion_times$first_report)
names(first_report_ordered) <- invasion_times$admin2

# save these to plot later
saveRDS(first_city_dist, "chik_first_city_dist.RDS")
saveRDS(first_report_ordered, "chik_first_report_ordered.RDS")

################## restart R

library(dplyr)

## ZIKV

# Import lat/long of cities
lat_lon <- readRDS("figures/data/lat_lon_complete.RDS")

# city invasion times
invasion_times <- readRDS('figures/data/zika_first_reported_cases_288.RDS') 
invasion_times <- rename(invasion_times, admin2 = admin2_code)

# subtract 1 from first reports to start at 0
invasion_times$first_report <- invasion_times$first_report - 1

# Merge latitudes and longitudes with invasion times by admin 2 
invasion_plot <- merge(lat_lon, invasion_times, by = "admin2", all.y = TRUE)

# import matrix of distance between cities
d_cities <- readRDS("figures/data/distance_matrix_km.RDS")

# Re-order rows and columns in distance matrix to match rows in invasion times df
row.order <- rownames(invasion_times)
d_cities <- d_cities[row.order, row.order]

# change NAs on diagonal to 0
d_cities[is.na(d_cities)] <- 0

# extract row of distances that correspond to the most likely origin of invasion (Barranquilla)
first_city_dist <- d_cities[14,]

# create named vector of invasion times
first_report_ordered <- as.numeric(invasion_times$first_report)
names(first_report_ordered) <- invasion_times$admin2

# save these to plot later
saveRDS(first_city_dist, "zika_first_city_dist.RDS")
saveRDS(first_report_ordered, "zika_first_report_ordered.RDS")

################### restart R

# load in saved files

## CHIKV
chik_first_city_dist <- readRDS("chik_first_city_dist.RDS")
chik_first_report_ordered <- readRDS("chik_first_report_ordered.RDS")

# correlations
cor.test(x = chik_first_city_dist, y = chik_first_report_ordered, method = "pearson", conf.level = 0.95)

## ZIKV
zika_first_city_dist <- readRDS("zika_first_city_dist.RDS")
zika_first_report_ordered <- readRDS("zika_first_report_ordered.RDS")

# correlations
cor.test(x = zika_first_city_dist, y = zika_first_report_ordered, method = "pearson", conf.level = 0.95)


par(mfrow=c(1,2))
par(mar=c(4,4,4,2))

# CHIKV
# plot invasion onset against straight line distance from Baranquilla for all cities with epidemic
plot(x = chik_first_city_dist, y = chik_first_report_ordered,
     xlab = "Great Circle Distance (km)", ylab = "Week of invasion", pch = 20, cex = 0.75, 
     cex.lab = 0.8, cex.axis = 0.8, xaxs="i", yaxs="i", xlim = c(0,1250), ylim = c(0,72))
abline(lm(chik_first_report_ordered~chik_first_city_dist), col="black") # regression line (y~x) 
mtext(expression(paste("Correlation: 0.37, p<0.0001")), line = 1, cex = 0.8)
mtext(text = "A", col = "black", cex = 1.2, side = 3, adj = 0, line = 2, font = (face=2))

par(mar=c(4,2,4,4))

# ZIKV

# plot invasion onset against straight line distance from San Andres for all cities with epidemic
plot(x = zika_first_city_dist, y = zika_first_report_ordered,
     xlab = "Great Circle Distance (km)", ylab = "", pch = 20, cex = 0.75, 
     cex.lab = 0.8, cex.axis = 0.8, xaxs="i", yaxs="i", xlim = c(0,1250), ylim = c(0,35))
abline(lm(zika_first_report_ordered~zika_first_city_dist), col="black") # regression line (y~x) 
mtext(expression(paste("Correlation: 0.24, p<0.0001")), line = 1, cex = 0.8)
mtext(text = "B", col = "black", cex = 1.2, side = 3, adj = 0, line = 2, font = (face=2))

