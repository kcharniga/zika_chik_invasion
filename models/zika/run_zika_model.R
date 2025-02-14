# Run Zika gravity model 

library(parallel)

# Use the detectCores() function to find the number of cores in system
no_cores <- detectCores() - 1

# source MCMC Functions
source("functions/gravity_model.R")
source("functions/get_params.R")
source("functions/calc_accept.R")
source("functions/calc_dic.R")
source("functions/calc_med_est.R")
source("functions/plot_res.R")

# burn In
burnIn = 100000

# load data
zika_N <- as.matrix(read.csv("data/zika_N.txt", 
                   sep='\t', 
                   stringsAsFactors = FALSE,
                   header = TRUE))

zika_I <- as.matrix(read.csv("data/zika_I.txt", 
                             sep='\t', 
                             stringsAsFactors = FALSE,
                             header = TRUE,
                             check.names = FALSE))

zika_infectious <- as.matrix(read.csv("data/zika_infectious.txt", 
                                      sep='\t', 
                                      stringsAsFactors = FALSE,
                                      header = TRUE,
                                      check.names = FALSE))

zika_J <- as.matrix(read.csv("data/zika_J.txt", 
                             sep='\t', 
                             stringsAsFactors = FALSE,
                             header = TRUE,
                             check.names = FALSE))

zika_d_cities <- as.matrix(read.csv("data/zika_d_cities.txt", 
                                    sep='\t', 
                                    stringsAsFactors = FALSE,
                                    header = TRUE,
                                    check.names = FALSE))

zika_H <- as.matrix(read.csv("data/zika_H.txt", 
                             sep='\t', 
                             stringsAsFactors = FALSE,
                             header = TRUE,
                             check.names = FALSE))

data1 <- list(N=zika_N, I=zika_I, infectious=zika_infectious, 
              J=zika_J, d_cities=zika_d_cities, H=zika_H)

# define start values and standard deviations of proposal distributions for 3 chains

vals <- list(
  data.frame(startvalue = c(3.0e-5, 0.5, 0.2, 0.1, 0.75),
             sd_proposal = c(0.3, 0.1, 0.5, 0.2, 1.0)),
  data.frame(startvalue = c(3.0e-4, 0.9, 0.6, 0.9, 1.2),
             sd_proposal = c(0.3, 0.1, 0.5, 0.2, 1.0)),
  data.frame(startvalue = c(6.0e-6, 1.3, 0.8, 0.7, 1.5),
             sd_proposal = c(0.3, 0.1, 0.5, 0.2, 1.0)))

# run 3 chains in parallel
start.time <- Sys.time()

result <- mclapply( vals,   # note that mcapply is for Mac. Can also use lapply
                    FUN = run_metropolis_MCMC,
                    mc.cores = no_cores,
                    iterations = 100000,
                    number_of_parameters = 5)

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

# get params
p <- get_params(result)

# calc acceptance rate for all params
ch1 <- calc_accept(chain_accept = p$chain_accept1, name_chain = 'chain1')
ch2 <- calc_accept(chain_accept = p$chain_accept2, name_chain = 'chain2')
ch3 <- calc_accept(chain_accept = p$chain_accept3, name_chain = 'chain3')

total_accept <- rbind(ch1, ch2, ch3)

# dic
dic_table <- calc_dic(burnIn, p)

cha1 <- calc_med_est(burnIn, chain_param = p$chain_param1, name_chain = 'chain1')
cha2 <- calc_med_est(burnIn, chain_param = p$chain_param2, name_chain = 'chain2')
cha3 <- calc_med_est(burnIn, chain_param = p$chain_param3, name_chain = 'chain3')

total_chain <- rbind(cha1, cha2, cha3)

# plots
plots <- plot_res(burnIn, p)

saveRDS(total_chain, "med_est_table.RDS")
saveRDS(dic_table, "dic_table.RDS")
saveRDS(total_accept, "accept_table.RDS")
saveRDS(p, "transformed_params.RDS")
