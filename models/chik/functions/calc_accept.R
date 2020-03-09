# correct acceptance rate for single update

calc_accept <- function(chain_accept, name_chain) {
  dd <- data.frame()  
  name_par   <- c('beta', 'gamma', 'nu', 'epsilon', "phi")
  iterations <- nrow(chain_accept)
  burnIn_accept <- 0.2*iterations
  for (i in seq_along(name_par)) {
    acceptance <- sum(chain_accept[-(1:burnIn_accept),i])/(iterations - burnIn_accept)
    accept_table <- data.frame(name_chain, name_par[i], acceptance)
    dd <- rbind(dd, accept_table)
  }
  rownames(dd) <- NULL
  names(dd) <- c('n_chain', 'name_par', 'accept_rate')
  return(dd)
}

