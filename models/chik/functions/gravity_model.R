########################
# Functions for MCMC
# Infectivity model with mu set to 0
# Singly update params (sequential proposal)
# 12-2-20
########################

#1 --------------------------------------  Define likelihood function

likelihood <- function(data1, param, kernel,k,num_kernel,inv_dist) {
  
  N <- data1$N
  I <- data1$I
  infectious <- data1$infectious
  J <- data1$J
  d_cities <- data1$d_cities
  H <- data1$H
  
  # parameters sampled on log scale because they are positive and can vary on large scale
  # b (beta) is transmission between cities
  # g (gamma) estimates spatial effect in transmission
  # v (nu) estimating effect of source pop size
  # m (mu) estimating effect of destination pop size
  # e (epsilon) estimating effect of spatial interaction
  # y (phi) estimating infectivity
  b <- as.numeric(param[1])
  g <- as.numeric(param[2])
  v <- as.numeric(param[3])
  m <- 0
  e <- as.numeric(param[4])
  y <- as.numeric(param[5])
  
  # Define kernel
  if (k!=1){
    if (k == 0){
      inv_dist <- (d_cities)^g 
      num_kernel <-  ((N^m)%*%t(N^v))
    }
    if (k == 2){
      inv_dist <- (d_cities)^g 
    }else{
      num_kernel <-  ((N^m)%*%t(N^v))
    }
    kernel <- num_kernel/inv_dist
    w <- ( colSums(matrix(N^v,length(N),length(N))/inv_dist ,na.rm = TRUE) )^e
    W <- matrix(1/w, nrow = length(N), ncol = length(N), byrow = TRUE)
    kernel <- kernel * W
  }
  
  # Change NAs back to zero for matrix multiplication
  kernel[is.na(kernel)] <- 0
  
  # Model equation
  foi = b*(H^y*infectious)%*%kernel
  
  log_Ptj <- J * (I*log(1-exp(-foi)) -foi * (1-I))   
  
  # Sum log likelihoods of Ptj over all susceptible cities to get conditional log likelihood
  LL <- sum(log_Ptj, na.rm = TRUE)
  
  return( list(LL = LL, kernel = kernel, num_kernel = num_kernel, inv_dist = inv_dist) )
  
}

#2 ---------------------------- Prior distribution - parameters can't be negative
prior <- function(param) {
  b <- as.numeric(param[1])
  g <- as.numeric(param[2])
  v <- as.numeric(param[3])
  e <- as.numeric(param[4])
  y <- as.numeric(param[5])
  b_prior <- dunif(b, min = 0, max = 40, log = TRUE)
  g_prior <- dunif(g, min = 0, max = 40, log = TRUE)
  v_prior <- dunif(v, min = 0, max = 40, log = TRUE)
  e_prior <- dunif(e, min = 0, max = 40, log = TRUE)
  y_prior <- dunif(y, min = 0, max = 40, log = TRUE)
  return(b_prior + g_prior + v_prior + e_prior + y_prior)
}

#3 ---------------------------- Posterior distribution 
posterior <- function(param,kernel,k,num_kernel,inv_dist) {
  temp <- likelihood(data1, param, kernel,k,num_kernel,inv_dist)
  kernel <- temp$kernel
  L <- temp$L
  pr <- prior(param)
  post <- L + pr
  
  return(list(L = L, pr = pr, post = post, kernel = kernel, num_kernel = temp$num_kernel, inv_dist = temp$inv_dist))
}

################### Metropolis algorithm ####################


#4 ---------------------------- 
proposal_function_rand <- function(param, sd_proposal) {
  
  # Sampling from proposal distribution
  
  return(rlnorm(1, meanlog = log(param), sdlog = sd_proposal)) # return one value at a time
}

#5 ---------------------------- 
run_metropolis_MCMC <- function(vals, number_of_parameters, iterations) {
  
  
  startvalue <- vals[ ,"startvalue"]
  sd_proposal <- vals[, "sd_proposal"]
  
  
  chain = array(dim = c((iterations*number_of_parameters)+1, number_of_parameters)) # Array for storing chain data
  colnames(chain) <- c("b", "g", "v", "e", "y")
  LL_chain = array(dim = c((iterations*number_of_parameters)+1, 1)) # Array for storing log likelihood
  
  accept_chain = array(dim = c(iterations, number_of_parameters)) # Array for storing acceptances.
  colnames(accept_chain) <- c("b", "g", "v", "e", "y")
  
  proposal <- startvalue
  
  # Start at a random parameter value
  chain[1, ] = startvalue 
  temp_post <- posterior(proposal, kernel=NA, k=0, num_kernel=NA, inv_dist=NA)
  LL_chain[1] <- temp_post$L
  
  
  # Choose a new parameter value close to the old value based on
  # probability density (proposal function) and calculate likelihood
  index <- 1
  for (i in 1:iterations) {
    
    for (j in 1:number_of_parameters) {
      
      proposal = chain[index,] # this line resets proposal to the chain, so rejected proposals aren't carried through. Next line updates just one param
      proposal[j] = proposal_function_rand(chain[index,j], sd_proposal[j])
      
      # Jump to new point proportional to likelihood of new value divided by likelihood of old value 
      # Note: p1/p2 = exp[log(p1) - log(p2)]
      # This is where acceptance probability is calculated
      probab1 = posterior(proposal, kernel = temp_post$kernel , k = j, 
                          num_kernel = temp_post$num_kernel, inv_dist = temp_post$inv_dist)
      
      current <- chain[index,j]
      proposed <- proposal[j]
      
      # correction for asymmetrical proposal distribution
      correction <- log(proposed) - log(current)
      
      # acceptance ratio 
      probab = probab1$post - temp_post$post + correction
      
      # If p > 1, use proposed value, if not, use the old value. Store acceptances
      if (log(runif(1)) < probab) {
        chain[index+1, ] = proposal
        accept_chain[i, j] = 1
        # Save log-likelihood
        LL_chain[index+1, ] <- probab1$L
        temp_post <- probab1
      } else {
        chain[index+1, ] = chain[index, ]
        accept_chain[i, j] = 0
        # Save log-likelihood
        LL_chain[index+1, ] <- temp_post$L
      }
      # Save log-likelihood
      #LL_chain[index+1, ] <- likelihood(data, chain[index+1, ])
      index = index + 1
    }
  }
  
  chain_matrix <- matrix(chain, ncol = number_of_parameters)
  
  both_outcomes <- list(chain = chain_matrix, LL_chain = LL_chain, accept = accept_chain)
  return(both_outcomes)
  
}



