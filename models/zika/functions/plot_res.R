##########
# Plots
##########

plot_res <- function(burnIn, p) {

  chain_param1 <- p$chain_param1 
  chain_param2 <- p$chain_param2
  chain_param3 <- p$chain_param3
  
  chain_loglik1 = p$chain_loglik1
  chain_loglik2 = p$chain_loglik2
  chain_loglik3 = p$chain_loglik3

  pdf("mcmc_plots.pdf")
  
  par(mfrow=c(4,3))
  par(mar=c(2, 2, 2.5, 1))
  
# Plot histograms of posterior distribution of parameters from 3 starting points

# beta
hist(chain_param1[-(1:burnIn),1], col = "darkturquoise", main = "beta - chain1", 
      xlab = "", ylab = "", cex.axis = 1.2, cex.main = 1.2)
hist(chain_param2[-(1:burnIn),1], col = "deeppink", main = 'beta - chain2')
hist(chain_param3[-(1:burnIn),1], col = "darkviolet", main = 'beta - chain3')

# gamma
hist(chain_param1[-(1:burnIn),2], col = "olivedrab4",
      main = "gamma - chain1", xlab = "", ylab = "",
      cex.axis = 1.2, cex.main = 1.2)
hist(chain_param2[-(1:burnIn),2], col = "thistle3", main = "gamma - chain2")
hist(chain_param3[-(1:burnIn),2], col = "peachpuff", main = "gamma - chain3")

# nu
hist(chain_param1[-(1:burnIn),3], col = "tomato2",
      main = "nu - chain1", xlab = "", ylab = "",
      cex.axis = 1.2, cex.main = 1.2)
hist(chain_param2[-(1:burnIn),3], col = "dodgerblue", main = "nu - chain2")
hist(chain_param3[-(1:burnIn),3], col = "darkseagreen4", main = "nu - chain3")

# Next page
par(mfrow=c(2,3))
par(mar=c(2, 2, 2.5, 1))

# epsilon
hist(chain_param1[-(1:burnIn),4], col = "grey37",
     main = "epsilon - chain1", xlab = "", ylab = "",
     cex.axis = 1.2, cex.main = 1.2)
hist(chain_param2[-(1:burnIn),4], col = "darkred", main = "epsilon - chain2")
hist(chain_param3[-(1:burnIn),4], col = "gold", main = "epsilon - chain3")

# phi
hist(chain_param1[-(1:burnIn),5], col = "slategray2",
     main = "phi - chain1", xlab = "", ylab = "",
     cex.axis = 1.2, cex.main = 1.2)
hist(chain_param2[-(1:burnIn),5], col = "tan2", main = "phi - chain2")
hist(chain_param3[-(1:burnIn),5], col = "chocolate4", main = "phi - chain3")

# Next page
par(mfrow=c(3,3))
par(mar=c(2, 2, 2.5, 1))

# Plot traces of chains and log likelihood from 3 starting points 

# Log likelihood
plot(chain_loglik1[-(1:burnIn),1], type = "l", main = "Chain values of log likelihood", 
      xlab = "", ylab = "", col = "burlywood4", cex.axis = 1.2, cex.main = 1.2, pch = 20)
lines(chain_loglik2[-(1:burnIn),1], col = "mediumaquamarine", pch = 20)
lines(chain_loglik3[-(1:burnIn),1], col = "sienna3", pch = 20)

# Chain values of beta
plot(chain_param1[-(1:burnIn),1], type = "l", main = "Chain values of beta",
      xlab = "", ylab = "", col = "darkturquoise", cex.axis = 1.2, cex.main = 1.2, pch = 20)
lines(chain_param2[-(1:burnIn),1], col = "deeppink", pch = 20)
lines(chain_param3[-(1:burnIn),1], col = "darkviolet", pch = 20)

# Chain values of gamma
plot(chain_param1[-(1:burnIn),2], type = "l", main = "Chain values of gamma",
      xlab = "", ylab = "", col = "olivedrab4", cex.axis = 1.2, cex.main = 1.2, pch = 20)
lines(chain_param2[-(1:burnIn),2], col = "thistle3", pch = 20)
lines(chain_param3[-(1:burnIn),2], col = "peachpuff", pch = 20)

# Chain values of nu
plot(chain_param1[-(1:burnIn),3], type = "l", main = "Chain values of nu",
      xlab = "", ylab = "", col = "tomato2", cex.axis = 1.2, cex.main = 1.2, pch = 20)
lines(chain_param2[-(1:burnIn),3], col = "dodgerblue", pch = 20)
lines(chain_param3[-(1:burnIn),3], col = "darkseagreen4", pch = 20)

# Chain values of epsilon
plot(chain_param1[-(1:burnIn),4], type = "l", main = "Chain values of epsilon",
     xlab = "", ylab = "", col = "grey37", cex.axis = 1.2, cex.main = 1.2, pch = 20)
lines(chain_param2[-(1:burnIn),4], col = "darkred", pch = 20)
lines(chain_param3[-(1:burnIn),4], col = "gold", pch = 20)

# Chain values of phi
plot(chain_param1[-(1:burnIn),5], type = "l", main = "Chain values of phi",
     xlab = "", ylab = "", col = "slategray2", cex.axis = 1.2, cex.main = 1.2, pch = 20)
lines(chain_param2[-(1:burnIn),5], col = "tan2", pch = 20)
lines(chain_param3[-(1:burnIn),5], col = "chocolate4", pch = 20)

#################################
# Autocorrelation of parameters?
#################################

# Next page
par(mfrow=c(3,2))
par(mar=c(4, 4, 4, 4))

plot(chain_param1[-(1:burnIn),1], chain_param1[-(1:burnIn),2], 
     col = "chartreuse2", xlab = "beta", ylab = "gamma", pch= 19,
      main = "Autocorrelation - chain1")


plot(chain_param1[-(1:burnIn),1], chain_param1[-(1:burnIn),4], 
     col = "red", xlab = "beta", ylab = "epsilon", pch= 19,
     main = "Autocorrelation - chain1")

plot(chain_param1[-(1:burnIn),2], chain_param1[-(1:burnIn),4], 
     col = "yellow", xlab = "gamma", ylab = "epsilon", pch= 19,
     main = "Autocorrelation - chain1")

plot(chain_param1[-(1:burnIn),1], chain_param1[-(1:burnIn),5], 
     col = "mediumpurple2", xlab = "beta", ylab = "phi", pch= 19,
     main = "Autocorrelation - chain1")

dev.off()



}

