##### Heatmaps

## CHIKV

# load disease and population data
chik <- readRDS("figures/data/chik_admin1_series.RDS")

depts <- readRDS("figures/data/admin1_names.RDS")

#Creating a column with correct (accented) department names

dati <- chik

dati$admin1_map <- NA
for (i in unique(dati$admin1)){
  admin1_map <- as.character(depts$admin1_map[depts$admin1_code==i])
  if (length(admin1_map) > 0) {
    dati$admin1_map[dati$admin1==i] <- admin1_map
  } else {
    dati$admin1_map[dati$admin1==i] <- NA}
}

#Load in population weighted centroids for all countries
adm1_cod <- read.csv("figures/data/pop_wt_centroids.txt", 
                     #Load Chik and population data
                     sep='\t', 
                     stringsAsFactors = FALSE,
                     header = FALSE, encoding="latin1")

#Subset Colombia information
adm1_cod <- adm1_cod[adm1_cod$V1=="Colombia",]

adm1_cod$V2[adm1_cod$V4 == 26] <- "San Andrés"

#Sort population weighted latitude (V8) from South to North
adm1_cod <- adm1_cod[order(adm1_cod$V8, decreasing=FALSE),]

# re-name variables
adm1_cod <- dplyr::rename(adm1_cod, country = V1, admin1_name = V2, country_id = V3, admin1_id = V4,
                          lon = V5, lat = V6, pop_wt_lon = V7, pop_wt_lat = V8, pop = V9)

#Change year from character to numeric
dati$year = as.numeric(as.character(dati$year))

#Assign week numbers over course of epidemic
first_year <- min(dati$year)
first_week <- min(dati$week[dati$year==first_year])

##Cases reported in week 53 of 2014

for (i in seq_len(nrow(dati))) {
  if (dati$year[i]==first_year) {
    dati$week_no[i] <- 1+(dati$week[i]-first_week)
  } else {
    dati$week_no[i] <-
      1+(53-first_week)+((dati$year[i]-(first_year+1))*53) + dati$week[i]
  }
}

no_weeks <- max(dati$week_no)

# fix 2016 weeks
dati1 <- dati[(dati$week_no < 84),]
dati2 <- dati[(dati$week_no > 84),]

dati2$week_no <- dati2$week_no - 1

dati <- rbind(dati1, dati2)

no_weeks <- max(dati$week_no)

#Create matrix for heatmap
hm <- matrix(nrow = nrow(adm1_cod), ncol = no_weeks)
for (wk in 1:no_weeks) {
  for (ad in 1:nrow(adm1_cod)) {
    cases <- dati$cases[dati$week_no == wk &
                          dati$admin1_map == adm1_cod$admin1_name[ad]]
    if (length(cases)==1) {
      hm[ad, wk] <- cases
    } else {
      hm[ad,wk] <- 0
    }
  }
}

colnames(hm)<-1:no_weeks
rownames(hm)<-adm1_cod$admin1_name

library(ggplot2)
library(dplyr)
library(reshape2)

# re-shape data for ggplot
hm2 <- melt(hm)

p <- ggplot(hm2, aes(Var2, Var1)) +
  geom_tile(aes(fill = value), colour = "white") +
  #scale_fill_gradient("Cases", low = "white", high = "steelblue") +
  scale_fill_viridis_c("Cases", option = "viridis", trans = "log", 
                       breaks = c(1, 10, 60, 500, 6800),labels = c(1, 10, 60, 500, 6800)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Week") +
  ylab(NULL) +
  theme(axis.text=element_text(size=16), axis.title.x=element_text(size=18)) +
  theme(legend.title=element_text(size=15), legend.text=element_text(size=14))

p      

saveRDS(p, "chik_heatmap_for_manuscript.RDS")


###################################### re-start R


## ZIKV

#Load disease and population data

zika <- readRDS("figures/data/zika_admin1_series.RDS")

depts <- readRDS("figures/data/admin1_names.RDS")

#Creating a column with correct (accented) department names

dati <- zika

dati$admin1_map <- NA
for (i in unique(dati$admin1)){
  admin1_map <- as.character(depts$admin1_map[depts$admin1_code==i])
  if (length(admin1_map) > 0) {
    dati$admin1_map[dati$admin1==i] <- admin1_map
  } else {
    dati$admin1_map[dati$admin1==i] <- NA}
}

#Load in population weighted centroids for all countries
adm1_cod <- read.csv("figures/data/pop_wt_centroids.txt", 
                     sep='\t', 
                     stringsAsFactors = FALSE,
                     header = FALSE, encoding="latin1")

#Subset Colombia information
adm1_cod <- adm1_cod[adm1_cod$V1=="Colombia",]

adm1_cod$V2[adm1_cod$V4 == 26] <- "San Andrés"

#Sort population weighted latitude (V8) from South to North
adm1_cod <- adm1_cod[order(adm1_cod$V8, decreasing=FALSE),]

# re-name variables
adm1_cod <- dplyr::rename(adm1_cod, country = V1, admin1_name = V2, country_id = V3, admin1_id = V4,
                          lon = V5, lat = V6, pop_wt_lon = V7, pop_wt_lat = V8, pop = V9)

#Assign week numbers over course of epidemic
first_year <- min(dati$year)
first_week <- min(dati$week[dati$year==first_year])

for (i in seq_len(nrow(dati))) {
  if (dati$year[i]==first_year) {
    dati$week_no[i] <- 1+(dati$week[i]-first_week)
  } else {
    dati$week_no[i] <-
      1+(52-first_week)+((dati$year[i]-(first_year+1))*52) + dati$week[i]
  }
}

no_weeks <- max(dati$week_no)

#Create matrix for heatmap
hm <- matrix(nrow = nrow(adm1_cod), ncol = no_weeks)
for (wk in 1:no_weeks) {
  for (ad in 1:nrow(adm1_cod)) {
    cases <- dati$cases[dati$week_no == wk &
                          dati$admin1_map == adm1_cod$admin1_name[ad]]
    if (length(cases)==1) {
      hm[ad, wk] <- cases
    } else {
      hm[ad,wk] <- 0
    }
  }
}

colnames(hm)<-1:no_weeks
rownames(hm)<-adm1_cod$admin1_name

library(ggplot2)
library(dplyr)
library(reshape2)

# re-shape data for ggplot
hm2 <- melt(hm)

p <- ggplot(hm2, aes(Var2, Var1)) +
  geom_tile(aes(fill = value), colour = "white") +
  #scale_fill_gradient("Cases", low = "white", high = "steelblue") +
  scale_fill_viridis_c("Cases", option = "viridis", trans = "log",
                       breaks=c(1, 7, 50, 300, 1500),labels=c(1, 7, 50, 300, 1500)) +
  #scale_fill_viridis_c("log(cases)", option = "viridis") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Week") +
  ylab(NULL) +
  theme(axis.text=element_text(size=16), axis.title.x=element_text(size=18)) +
  theme(legend.title=element_text(size=15), legend.text=element_text(size=14)) 

p

saveRDS(p, "zika_heatmap_for_manuscript.RDS")

############################################### re-start R

# put them together

zika_hm <- readRDS("zika_heatmap_for_manuscript.RDS")
chik_hm <- readRDS("chik_heatmap_for_manuscript.RDS")

library("cowplot")
plot_grid(chik_hm, zika_hm, 
          labels = c("A", "B"),
          label_size = 20,
          ncol = 1, nrow = 2)

# save as pdf, resolution 1000 x 1200

