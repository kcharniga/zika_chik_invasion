###########################################################################
### Median posterior estimates of parameters and 95% credible intervals ###
### for each chain
###########################################################################

calc_med_est <- function(burnIn, chain_param, name_chain) {
  
dd <- data.frame()  
name_par   <- c('beta', 'gamma', 'nu', 'epsilon', "phi")
for (i in seq_along(name_par)) {
med <- median(chain_param[-(1:burnIn),i])
ci <- quantile(chain_param[-(1:burnIn),i], c(0.025, 0.975)) # for a 95% interval
med_est_table <- data.frame(name_chain, name_par[i], med, ci[1], ci[2])
dd <- rbind(dd, med_est_table)
}
rownames(dd) <- NULL
names(dd) <- c('n_chain', 'name_par', 'median_par', 'lowerci_par', 'upperci_par')
return(dd)
}

