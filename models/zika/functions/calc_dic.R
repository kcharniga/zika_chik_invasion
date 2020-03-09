

calculateDIC <- function(chain_param, 
                         chain_loglik) {
  
  #Calculate L
  theta_hat = apply(chain_param, 2, median) # apply median to each of the columns in chain param
  L = likelihood(data1, theta_hat, kernel=NA, k = 0, num_kernel=NA) 
  
  ll_med = median(chain_loglik)
  P = 2 * (L$LL - ll_med)
  
  #Calculate DIC
  DIC = -2 * (L$LL - P)
  
  #Return the results
  return(list(DIC=DIC, P=P, L=L$LL))
}


calc_dic <- function(burnIn, p) {
  
  chain_param1 <- p$chain_param1 
  chain_param2 <- p$chain_param2
  chain_param3 <- p$chain_param3
  
  chain_loglik1 = p$chain_loglik1
  chain_loglik2 = p$chain_loglik2
  chain_loglik3 = p$chain_loglik3
  
  dic1 <- calculateDIC(chain_param = chain_param1[-(1:burnIn), ], 
                       chain_loglik = chain_loglik1[-(1:burnIn),1])
  
  dic2 <- calculateDIC(chain_param = chain_param2[-(1:burnIn), ], 
                       chain_loglik = chain_loglik2[-(1:burnIn),1])
  
  dic3 <- calculateDIC(chain_param = chain_param3[-(1:burnIn), ], 
                       chain_loglik = chain_loglik3[-(1:burnIn),1])
  
  dic <- data.frame(rbind(dic1, dic2, dic3))
  
  return(dic)
  
}

