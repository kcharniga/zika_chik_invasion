##### Figure 4 Probability distribution of first reported cases 

# 12-8-20

## CHIKV


# invasion times
invasion_week <- read.csv("figures/data/chik_first_reported_cases_338.txt", # invasion times start at week 1
                          sep='\t', 
                          stringsAsFactors = FALSE,
                          header = TRUE)

# format and pad admin2 with zeros
invasion_week$admin2_code <- sprintf("%05d", invasion_week$admin2_code)
invasion_week$admin2_code <- as.character(invasion_week$admin2_code) 

invasion_week <- invasion_week[order(invasion_week$first_report), ]

# load data
N <- as.matrix(read.csv("figures/data/chik_N.txt", 
                        sep='\t', 
                        stringsAsFactors = FALSE,
                        header = TRUE))

I <- as.matrix(read.csv("figures/data/chik_I.txt", 
                        sep='\t', 
                        stringsAsFactors = FALSE,
                        header = TRUE,
                        check.names = FALSE))

infectious <- as.matrix(read.csv("figures/data/chik_infectious.txt", 
                                 sep='\t', 
                                 stringsAsFactors = FALSE,
                                 header = TRUE,
                                 check.names = FALSE))

J <- as.matrix(read.csv("figures/data/chik_J.txt", 
                        sep='\t', 
                        stringsAsFactors = FALSE,
                        header = TRUE,
                        check.names = FALSE))

d_cities <- as.matrix(read.csv("figures/data/chik_d_cities.txt", 
                               sep='\t', 
                               stringsAsFactors = FALSE,
                               header = TRUE,
                               check.names = FALSE))

H <- as.matrix(read.csv("figures/data/chik_H.txt", 
                        sep='\t', 
                        stringsAsFactors = FALSE,
                        header = TRUE,
                        check.names = FALSE))


################################################
# prepare data to plot probability 
# distributions for central parameter estimates
################################################

# median parameter estimates from MCMC
b <- 0.24
g <- 1.68
v <- 0.65
m <- 0
e <- 0.83
y <- 0.35

# Define kernel

inv_dist <- (d_cities)^g 
num_kernel <-  ((N^m)%*%t(N^v))

kernel <- num_kernel/inv_dist

w <- ( colSums(matrix(N^v,length(N),length(N))/inv_dist ,na.rm = TRUE) )^e
W <- matrix(1/w, nrow = length(N), ncol = length(N), byrow = TRUE)
kernel <- kernel * W

# Change NAs back to zero for matrix multiplication
kernel[is.na(kernel)] <- 0

# Model equation
foi = b*(H^y*infectious)%*%kernel

Ptj2 <- matrix(NA, 69, 338)
logPtj <- matrix(NA, 69, 338)
I2 <- matrix(1, 69, 338)
J2 <- matrix(NA, 69, 338)

for (i in 2:69){ # looping over all weeks
  
  
  I2[1:(i-1),] <- 0 # shift time of first invasion
  J2[1:i,] <- 1 # also shifting
  log_Ptj2 <- J2*(I2*log(1-exp(-foi)) -foi * (1-I2)) 
  temp <- colSums(log_Ptj2, na.rm = TRUE)
  
  
  Ptj2[i,] <- exp(temp)
}


# Set NAs to 0
Ptj2[is.na(Ptj2)] <- 0

# normalize Ptj
Ptj_norm2 <- LICORS::normalize(Ptj2, byrow = FALSE)

prob <- Ptj_norm2

colnames(prob) <- invasion_week$admin2_code
rownames(prob) <- 1:69

# set first column to 0
prob[,1] <- 0

dat1 <- prob

dat2 <- matrix(NA, 69, 338)

for (i in 1:338) {
  dat2[,i] <- cumsum(dat1[,i])
}

colnames(dat2) <- invasion_week$admin2_code

dat3 <- matrix(0, 69, 338)

invasion_week_v <- as.vector(invasion_week$first_report)

for (i in 1:338){
  inf_wk <- invasion_week_v[i]
  dat3[inf_wk,i] <- 1
}

colnames(dat3) <- invasion_week$admin2_code

# load packages
library(reshape2)
library(ggplot2)

# melt data into form acceptable for ggplot
datt1 <- data.frame(dat1)
colnames(datt1) <- colnames(dat1)
datt1$week <- 1: nrow(datt1)
dat1_m <- melt(datt1, id = 'week', variable.name = 'city', value.name = 'prob_ix')

datt2 <- data.frame(dat2)
colnames(datt2) <- colnames(dat2)
datt2$week <- 1: nrow(datt2)
dat2_m <- melt(datt2, id = 'week', variable.name = 'city', value.name = 'cumm_p_ix')

# determine upper and lower limits
dattl <- data.frame(dat2)
colnames(dattl) <- colnames(dat2)
lower_week <- dattl[1,]
lower_week[1, ] <- NA
for (i  in 1:NCOL(dattl)) {
  lower_week[,i] <- (which(dattl[,i] > 0.025)[1] - 0.5)
}

upper_week <- dattl[1,]
upper_week[1, ] <- NA
for (i  in 1:NCOL(dattl)) {
  upper_week[,i] <- (which(dattl[,i] > 0.975)[1] - 0.5)
}

limitsu <- melt(upper_week, variable.name = 'city', value.name = 'upper')
limitsl <- melt(lower_week, variable.name = 'city', value.name = 'lower')
limits <- merge(limitsl, limitsu, by = 'city')

# more melting data
datt3 <- data.frame(dat3)
colnames(datt3) <- colnames(dat3)
datt3$week <- 1: nrow(datt3)
dat3_m <- melt(datt3, id = 'week', variable.name = 'city', value.name = 'observed_ix')

# put all the data together to plot
datf <- merge(dat1_m, dat2_m, by = c('week', 'city'))
datf <- merge(datf, dat3_m, by = c('week', 'city'))

datf$observed_ix[datf$observed_ix == 1] <- -0.001 # assign value to invasion week so it is plotted just under plot
datf$observed_ix[datf$observed_ix == 0] <- NA # make the 0's NA so they are not plotted

saveRDS(datf, 'chik_validation_first_reports_338.RDS')

################# re-start R

## ZIKV


# invasion times
invasion_week <- read.csv("figures/data/zika_first_reported_cases_288.txt", # weeks start at 1
                          sep='\t', 
                          stringsAsFactors = FALSE,
                          header = TRUE)

# format and pad admin2 with zeros
invasion_week$admin2_code <- sprintf("%05d", invasion_week$admin2_code)
invasion_week$admin2_code <- as.character(invasion_week$admin2_code) 

invasion_week <- invasion_week[order(invasion_week$first_report), ]

# load data
N <- as.matrix(read.csv("figures/data/zika_N.txt", 
                        sep='\t', 
                        stringsAsFactors = FALSE,
                        header = TRUE))

I <- as.matrix(read.csv("figures/data/zika_I.txt", 
                        sep='\t', 
                        stringsAsFactors = FALSE,
                        header = TRUE,
                        check.names = FALSE))

infectious <- as.matrix(read.csv("figures/data/zika_infectious.txt", 
                                 sep='\t', 
                                 stringsAsFactors = FALSE,
                                 header = TRUE,
                                 check.names = FALSE))

J <- as.matrix(read.csv("figures/data/zika_J.txt", 
                        sep='\t', 
                        stringsAsFactors = FALSE,
                        header = TRUE,
                        check.names = FALSE))

d_cities <- as.matrix(read.csv("figures/data/zika_d_cities.txt", 
                               sep='\t', 
                               stringsAsFactors = FALSE,
                               header = TRUE,
                               check.names = FALSE))

H <- as.matrix(read.csv("figures/data/zika_H.txt", 
                        sep='\t', 
                        stringsAsFactors = FALSE,
                        header = TRUE,
                        check.names = FALSE))

################################################
# prepare data to plot probability 
# distributions for central parameter estimates
################################################

b <- 1.10
g <- 1.74
v <- 0.55
m <- 0
e <- 0.68
y <- 0.27

# Define kernel
inv_dist <- (d_cities)^g 
num_kernel <-  ((N^m)%*%t(N^v))

kernel <- num_kernel/inv_dist

w <- ( colSums(matrix(N^v,length(N),length(N))/inv_dist ,na.rm = TRUE) )^e
W <- matrix(1/w, nrow = length(N), ncol = length(N), byrow = TRUE)
kernel <- kernel * W

# Change NAs back to zero for matrix multiplication
kernel[is.na(kernel)] <- 0

# Model equation
foi = b*(H^y*infectious)%*%kernel

Ptj2 <- matrix(NA, 34, 288)
logPtj <- matrix(NA, 34, 288)
I2 <- matrix(1, 34, 288)
J2 <- matrix(NA, 34, 288)

for (i in 2:34){ # looping over all weeks
  
  
  I2[1:(i-1),] <- 0 # shift time of first invasion
  J2[1:i,] <- 1 # also shifting
  log_Ptj2 <- J2*(I2*log(1-exp(-foi)) -foi * (1-I2)) 
  temp <- colSums(log_Ptj2, na.rm = TRUE)
  
  
  Ptj2[i,] <- exp(temp)
}

# Set NAs to 0
Ptj2[is.na(Ptj2)] <- 0

# normalize Ptj
Ptj_norm2 <- LICORS::normalize(Ptj2, byrow = FALSE)

prob_median_est <- Ptj_norm2

prob <- prob_median_est
#prob <- prob_avg

colnames(prob) <- invasion_week$admin2_code
rownames(prob) <- 1:34

# set first five columns to 0 (invaded in week 1)
prob[,1] <- 0
prob[,2] <- 0
prob[,3] <- 0
prob[,4] <- 0
prob[,5] <- 0

dat1 <- prob

dat2 <- matrix(NA, 34, 288)

for (i in 1:288) {
  dat2[,i] <- cumsum(dat1[,i])
}

colnames(dat2) <- invasion_week$admin2_code

dat3 <- matrix(0, 34, 288)

invasion_week_v <- as.vector(invasion_week$first_report)

for (i in 1:288){
  inf_wk <- invasion_week_v[i]
  dat3[inf_wk,i] <- 1
}

colnames(dat3) <- invasion_week$admin2_code

# load packages
library(reshape2)
library(ggplot2)

# melt data into form acceptable for ggplot
datt1 <- data.frame(dat1)
colnames(datt1) <- colnames(dat1)
datt1$week <- 1: nrow(datt1)
dat1_m <- melt(datt1, id = 'week', variable.name = 'city', value.name = 'prob_ix')

datt2 <- data.frame(dat2)
colnames(datt2) <- colnames(dat2)
datt2$week <- 1: nrow(datt2)
dat2_m <- melt(datt2, id = 'week', variable.name = 'city', value.name = 'cumm_p_ix')

# determine upper and lower limits
dattl <- data.frame(dat2)
colnames(dattl) <- colnames(dat2)
lower_week <- dattl[1,]
lower_week[1, ] <- NA
for ( i  in 1: NCOL(dattl)) {
  lower_week[,i] <- (which(dattl[,i] > 0.025)[1] - 0.5)
}

upper_week <- dattl[1,]
upper_week[1, ] <- NA
for ( i  in 1: NCOL(dattl)) {
  upper_week[,i] <- (which(dattl[,i] > 0.975)[1] - 0.5)
}

limitsu <- melt(upper_week, variable.name = 'city', value.name = 'upper')
limitsl <- melt(lower_week, variable.name = 'city', value.name = 'lower')
limits <- merge(limitsl, limitsu, by = 'city')

# more melting data  
datt3 <- data.frame(dat3)
colnames(datt3) <- colnames(dat3)
datt3$week <- 1: nrow(datt3)
dat3_m <- melt(datt3, id = 'week', variable.name = 'city', value.name = 'observed_ix')

# put all the data together to plot
datf <- merge(dat1_m, dat2_m, by = c('week', 'city'))
datf <- merge(datf, dat3_m, by = c('week', 'city'))

datf$observed_ix[datf$observed_ix == 1] <- -0.001 # assign value to invasion week so it is plotted just under plot
datf$observed_ix[datf$observed_ix == 0] <- NA # make the 0's NA so they are not plotted

saveRDS(datf, 'zika_validation_first_reports_288.RDS')

################### re-start R

library(tidyverse)
library(cowplot)
library(aweek)

# 16-3-20

zika <- readRDS("zika_validation_first_reports_288.RDS") %>% mutate(disease = 'ZIKA') %>% arrange(observed_ix)

chik <- readRDS("chik_validation_first_reports_338.RDS") %>% mutate(disease = 'CHIKV') %>% arrange(observed_ix)

zika$w_intro <- zika$week - 1 # start at week 0
chik$w_intro <- chik$week - 1 # start at week 0

#--- Dates

min_date_zika_intro <- as.Date('08/08/2015', format = "%d/%m/%Y") 
more_recent_date    <- as.Date('03/04/2019', format = "%d/%m/%Y") 

zika_dates_sequence <- seq(min_date_zika_intro, more_recent_date, by = 1)
zika_weeks_sequence <- unique(date2week(zika_dates_sequence, 'Sunday', floor_day = TRUE))
dates_zika <- data.frame(w_intro = min(zika$w_intro): (length(zika_weeks_sequence)-1))  %>%
  mutate(epiw_intro = zika_weeks_sequence)


#---  Zika
zika <- zika %>% 
  arrange(w_intro) %>% 
  merge(dates_zika, by = 'w_intro') %>%
  mutate(date = week2date(epiw_intro, week_start = 'Sunday', floor_day = TRUE))


#--- Dates

min_date_chik_intro <- as.Date('31/05/2014', format = "%d/%m/%Y")  
more_recent_date    <- as.Date('03/04/2019', format = "%d/%m/%Y") 

chik_dates_sequence <- seq(min_date_chik_intro, more_recent_date, by = 1)
chik_weeks_sequence <- unique(date2week(chik_dates_sequence, 'Sunday', floor_day = TRUE))
dates_chik <- data.frame(w_intro = min(chik$w_intro): (length(chik_weeks_sequence)-1))  %>%
  mutate(epiw_intro = chik_weeks_sequence)


#---  Chik
chik <- chik %>% 
  arrange(w_intro) %>% 
  merge(dates_chik, by = 'w_intro') %>%
  mutate(date = week2date(epiw_intro, week_start = 'Sunday', floor_day = TRUE))


# ---- Zika plot
# change scale so better contrast on plots
zikap <- zika
zikap$prob_ix100 <- zikap$prob_ix * 100
zikap$prob_ix100[zikap$prob_ix100 <1] <- 1
daft_ix_zika <- filter(zikap, observed_ix == -0.001) # to add observed inf week on top of heatmap

max_value_scale <- 100
addz <- zikap[1,]
addz$prob_ix100 <- max_value_scale
zikap <- rbind(addz, zikap)

pzika <- ggplot(zikap, 
                aes(date, city)) + 
  geom_tile(aes(fill = prob_ix100), colour = "white") + 
  theme_bw(20) +
  geom_point(data = daft_ix_zika, aes(x=date, y=reorder(city, - prob_ix100)), colour = 'black', size = 0.4) +
  scale_fill_viridis_c(option = "viridis", direction = -1,  trans = "log", 
                       breaks = c(1, 10, 55), labels = c(0.01, 0.10, 0.55)) +
  scale_x_date(date_breaks = "2 months", date_labels =  "%b %Y",
               expand = c(0,0)) +
  # scale_fill_viridis_c(option = "viridis", direction = -1,  limits = c(0,50)) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y = element_blank()) +
  facet_wrap(~disease, scale='free_y', nrow = 2) +
  labs(fill = "") +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +
  theme(plot.margin = unit(c(1.2, 0, 0, 0), "cm")) +
  labs(x = "Date", y = "City")
  




# ---- Chik plot
chikp <- chik
chikp$prob_ix100 <- chikp$prob_ix * 100
chikp$prob_ix100[chikp$prob_ix100 <1] <- 1
daft_ix_chik <- filter(chikp, observed_ix == -0.001) # to add observed inf week on top of heatmap

addc <- chikp[98,]
addc$prob_ix100 <- max_value_scale
chikp <- rbind(addc, chikp)
pchik <- ggplot(chikp, 
                aes(date, city)) + 
  geom_tile(aes(fill = prob_ix100), colour = "white") + 
  theme_bw(20) +
  geom_point(data = daft_ix_chik, aes(x=date, y=reorder(city, - prob_ix100)), colour = 'black', size = 0.4) +
  scale_fill_viridis_c(option = "viridis", direction = -1,  trans = "log",
                       begin = 0, end = 1) +
  scale_x_date(date_breaks = "3 months", date_labels =  "%b %Y",
               expand = c(0,0)) +
  # scale_fill_viridis_c(option = "viridis", direction = -1,  limits = c(0,50)) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y = element_blank()) +
  facet_wrap(~disease, scale='free_y', nrow = 2) +
  theme(legend.position = 'none') +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +
  theme(plot.margin = unit(c(1.2, 0, 0, 0), "cm")) +
  labs(x = "Date", y = "City")


# ---- plot both diseases

png(filename = "figure4.png",
    height = 800, width = 600*2, bg = "white", units = "px")
plot_grid(pchik, pzika, nrow = 1, rel_widths = c(0.94, 1.12),
          labels = c("A", "B"),
          label_size = 20)

dev.off()

