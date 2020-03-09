## CHIKV

# Read in posterior distribution of parameters for best model
# best model = infectivity model with mu set to 0
# straight line distance (no elevation)
# using first reported cases in each city (335 cities)
params <- readRDS("figures/data/chik_mcmc_chains.RDS")
params <- params$chain_param1 # first chain
burnIn <- 100000
params <- params[-(1:burnIn),] # remove burnIn

# randomly sample 1000 parameter sets from posterior distribution
sampled_params <- params[sample(nrow(params), size = 1000, replace = FALSE), ]

# invasion times
invasion_week <- readRDS("figures/data/chik_first_reported_cases_335.RDS")
invasion_week <- invasion_week[order(invasion_week$first_report), ]

# load data
#data <- readRDS("figures/data/matrices_chik_first_reported_cases_335.RDS")

# load data with bigger matrices - more time steps
data <- readRDS("figures/data/chik_matrices_93_weeks.RDS")

###################################################
# Epidemic simulations for 1000 parameter estimates
###################################################

# initialize array, giving dimensions
my_I2_array <- array(0, dim=c(93, 335, 1000))

# this takes a few minutes (~10)
start.time <- Sys.time()

for (i in 1:nrow(sampled_params)) { 
  
  N <- data$N
  I <- data$I
  infectious <- data$infectious
  J <- data$J
  d_cities <- data$d_cities
  H <- data$H
  
  # model
  b <- sampled_params[i,1]
  g <- sampled_params[i,2]
  v <- sampled_params[i,3]
  m <- 0
  e <- sampled_params[i,4]
  y <- sampled_params[i,5]
  
  # Define kernel
  
  inv_dist <- (d_cities)^g 
  num_kernel <-  ((N^m)%*%t(N^v))
  
  kernel <- num_kernel/inv_dist
  
  w <- ( colSums(matrix(N^v,length(N),length(N))/inv_dist ,na.rm = TRUE) )^e
  W <- matrix(1/w, nrow = length(N), ncol = length(N), byrow = TRUE)
  kernel <- kernel * W
  
  # Change NAs back to zero for matrix multiplication
  kernel[is.na(kernel)] <- 0
  
  # define new matrices in order to simulate epidemic 
  Ptj_better <- matrix(NA, 93, 335)
  
  I2 <- matrix(0, 93, 335) 
  I2[,1] <- 1 # start with one city invaded in week 1
  
  # make infectious matrix
  infectious <- I2 
  first_row <- matrix(data = 0, nrow = 1, ncol = 335) # add row of zeros
  infectious <- rbind(first_row, infectious)
  infectious <- infectious[1:93,]  # get rid of last row
  
  J2 <- infectious
  J2[J2==1] <- NA
  J2[J2==0] <- 1
  J2[1,1] <- NA # set first city to NA for 1st time step
  
  ################################################
  
  H_sim <- matrix(0, nrow(H), ncol(H))
  H_sim[,1] <- H[,1]
  
  for (k in 2:93){ # looping over all weeks
    
    # calculate foi at time t
    foi = b*(H_sim^y*infectious)%*%kernel
    
    Ptj_better <- J2[k,]*(1-exp(-foi[k,])) * (1-I2[k,]) # remove part of likelihood where city escapes invasion
    
    dog <- runif(335) < Ptj_better
    # for all uninvaded cities, draw a random number r from uniform 0-1
    
    
    I2[k:93, dog] <-  1      # city is invaded
    f <- which(dog %in% 1)
    for (j in f) {
      if (j <= nrow(invasion_week)) {
        infectiousness <- H[invasion_week$first_report[j]:93,j]   # shift H 
        index <- k:(k+length(infectiousness) - 1)
        temp <- index <= 93
        H_sim[index[temp], j] <- infectiousness[temp]
      } 
      
    }
    
    # update matrices with new invasion times
    infectious <- I2 
    first_row <- matrix(data = 0, nrow = 1, ncol = 335) # add row of zeros
    infectious <- rbind(first_row, infectious)
    infectious <- infectious[1:93,]  # get rid of last row
    
    J2 <- infectious
    J2[J2==1] <- NA
    J2[J2==0] <- 1
    J2[1,1] <- NA # set first city to NA for 1st time step
  }  
  
  # add I2 to an array
  my_I2_array[,,i] <- I2
  
}

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken  

saveRDS(I, "chik_Ifirstreports.RDS")
saveRDS(my_I2_array, "chik_my_I2_arrayfirstreports.RDS")

################### re-start R

## ZIKV

# Read in posterior distribution of parameters for best model
# best model = infectivity model with mu set to 0
# straight line distance (no elevation)
# using first reported cases in each city (288 cities)
params <- readRDS("figures/data/zika_mcmc_chains.RDS")
params <- params$chain_param1
burnIn <- 100000
params <- params[-(1:burnIn),] # remove burnIn

# invasion times
invasion_week <- readRDS("figures/data/zika_first_reported_cases_288.RDS")
invasion_week <- invasion_week[order(invasion_week$first_report), ]

# # randomly sample 1000 parameter sets from posterior distribution
sampled_params <- params[sample(nrow(params), size = 1000, replace = FALSE), ]

# load data
#data <- readRDS("figures/data/matrices_zika_first_reported_cases_288.RDS")

# load data with bigger matrices to run for more simulations
data <- readRDS("figures/data/zika_matrices_50_weeks.RDS")

###################################################
# Epidemic simulations for 1000 parameter estimates
###################################################

# initialize array, giving dimensions
my_I2_array <- array(0, dim=c(50, 288, 1000))

# takes about 2 minutes
start.time <- Sys.time()

for (i in 1:nrow(sampled_params)) { 
  
  N <- data$N
  I <- data$I
  infectious <- data$infectious
  J <- data$J
  d_cities <- data$d_cities
  H <- data$H
  
  # model
  b <- sampled_params[i,1]
  g <- sampled_params[i,2]
  v <- sampled_params[i,3]
  m <- 0
  e <- sampled_params[i,4]
  y <- sampled_params[i,5]
  
  # Define kernel
  
  inv_dist <- (d_cities)^g 
  num_kernel <-  ((N^m)%*%t(N^v))
  
  kernel <- num_kernel/inv_dist
  
  w <- ( colSums(matrix(N^v,length(N),length(N))/inv_dist ,na.rm = TRUE) )^e
  W <- matrix(1/w, nrow = length(N), ncol = length(N), byrow = TRUE)
  kernel <- kernel * W
  
  # Change NAs back to zero for matrix multiplication
  kernel[is.na(kernel)] <- 0
  
  # define new matrices in order to simulate epidemic
  Ptj_better <- matrix(NA, 50, 288)
  
  I2 <- matrix(0, 50, 288) 
  I2[,(1:5)] <- 1 ############# start with first 5 cities invaded
  
  # make infectious matrix
  infectious <- I2 
  first_row <- matrix(data = 0, nrow = 1, ncol = 288) # add row of zeros
  infectious <- rbind(first_row, infectious)
  infectious <- infectious[1:50,]  # get rid of last row
  
  J2 <- infectious
  J2[J2==1] <- NA
  J2[J2==0] <- 1
  # set first five cities to NA for 1st time step
  J2[1,1] <- NA
  J2[1,2] <- NA
  J2[1,3] <- NA
  J2[1,4] <- NA
  J2[1,5] <- NA
  
  ################################################
  # 
  
  H_sim <- matrix(0, nrow(H), ncol(H)) ####### take infecitivity of first 5 invaded cities
  H_sim[,(1:5)] <- H[,(1:5)]
  
  for (k in 2:50){ # looping over all time points
    
    # calculate foi at time t
    foi = b*(H_sim^y*infectious)%*%kernel
    
    Ptj_better <- J2[k,]*(1-exp(-foi[k,])) * (1-I2[k,]) # remove part where city escapes invasion
    
    dog <- runif(288) < Ptj_better
    # for all uninvaded cities, draw a random number r from uniform 0-1
    
    
    I2[k:50, dog] <-  1      # city is invaded
    f <- which(dog %in% 1)
    for (j in f) {
      if (j <= nrow(invasion_week)) {
        infectiousness <- H[invasion_week$first_report[j]:50,j]   # shift H 
        index <- k:(k+length(infectiousness) - 1)
        temp <- index <= 50
        H_sim[index[temp], j] <- infectiousness[temp]
      } 
      
    }
    
    # update matrices with new invasion times
    infectious <- I2 
    first_row <- matrix(data = 0, nrow = 1, ncol = 288) # add row of zeros
    infectious <- rbind(first_row, infectious)
    infectious <- infectious[1:50,]  # get rid of last row
    
    J2 <- infectious
    J2[J2==1] <- NA
    J2[J2==0] <- 1
    # set first five cities to NA for 1st week 
    J2[1,1] <- NA
    J2[1,2] <- NA
    J2[1,3] <- NA
    J2[1,4] <- NA
    J2[1,5] <- NA
  }  
  
  # add I2 to an array
  my_I2_array[,,i] <- I2
  
}

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

saveRDS(I, "zika_Ifirstreports.RDS")
saveRDS(my_I2_array, "zika_my_I2_arrayfirstreports.RDS")

############### re-start R

## Plot

chik_I <- readRDS("chik_Ifirstreports.RDS")
chik_my_I2_array <- readRDS("chik_my_I2_arrayfirstreports.RDS")
chik_sim_avg <- apply(chik_my_I2_array, c(1, 2), mean, na.rm = TRUE) # compute average across 1,000 simulations

zika_I <- readRDS("zika_Ifirstreports.RDS")
zika_my_I2_array <- readRDS("zika_my_I2_arrayfirstreports.RDS")
zika_sim_avg <- apply(zika_my_I2_array, c(1, 2), mean, na.rm = TRUE) # compute average across 1,000 simulations

par(mfrow=c(1,2))
par(mar=c(4,4,3,2))

chik_weeks <- as.vector(0:92, mode = "integer") # start at week 0
zika_weeks <- as.vector(0:49, mode = "integer") # start at week 0

#CHIKV

plot(x = chik_weeks, y = rowSums(chik_I), main = "", xlab = "Weeks", ylab = "Number of Invaded Cities", type = "l", col = "red", 
     lwd = 1)
for (h in 1:1000) { # number of param sets
  lines(x = chik_weeks, y = rowSums(chik_my_I2_array[,,h]), col = "grey76")
}
lines(x = chik_weeks, y = rowSums(chik_sim_avg), col = "grey41", lwd = 3)
lines(x = chik_weeks, y = rowSums(chik_I), col = "red", lwd = 3)
mtext(text = "A", col = "black", cex = 1.2, side = 3, adj = 0, line = 1, font = (face=2))


par(mar=c(4,2,3,4))

# ZIKV
plot(x = zika_weeks, y = rowSums(zika_I), main = "", xlab = "Weeks", type = "l", col = "red", lwd = 1)
for (h in 1:1000) { # number of param sets
  lines(x = zika_weeks, y = rowSums(zika_my_I2_array[,,h]), col = "grey76")
}
lines(x = zika_weeks, y = rowSums(zika_sim_avg), col = "grey41", lwd = 3)
lines(x = zika_weeks, y = rowSums(zika_I), col = "red", lwd = 3)
mtext(text = "B", col = "black", cex = 1.2, side = 3, adj = 0, line = 1, font = (face=2))


