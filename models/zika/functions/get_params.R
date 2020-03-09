
get_params <- function(result) {

# extract chains and log likelihood from result
chain1 <- result[[1]]
chain2 <- result[[2]]
chain3 <- result[[3]]
  
# Get parameter values from MCMC chains
chain_param1 = chain1$chain
chain_param2 = chain2$chain
chain_param3 = chain3$chain

# Get log likelihood
chain_loglik1 = chain1$LL_chain
chain_loglik2 = chain2$LL_chain
chain_loglik3 = chain3$LL_chain

# Get acceptance rates
chain_accept1 = chain1$accept
chain_accept2 = chain2$accept
chain_accept3 = chain3$accept


pp <- list(chain_param1 = chain_param1, 
          chain_param2 = chain_param2, 
          chain_param3 = chain_param3,
          
          chain_loglik1 = chain_loglik1, 
          chain_loglik2 = chain_loglik2, 
          chain_loglik3 = chain_loglik3,

          chain_accept1 = chain_accept1,
          chain_accept2 = chain_accept2,
          chain_accept3 = chain_accept3)
          

return(pp)

}



