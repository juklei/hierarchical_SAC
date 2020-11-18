## hsac_collector_raw
##
## First edit: 20201115
## Last edit: 20201115
##
## Author: Julian Klein

model{

  ## Likelihood:
  
  ## Observation model:
  for(k in 1:nrep){
    for(s in 1:nsite){
      for(i in 1:ntree[s]){
      obs_all[i,k,s] ~ dpois(lambda_obs[i,k,s])
      sim_obs[i,k,s] ~ dpois(lambda_obs[i,k,s])
      lambda_obs[i,k,s] <- (gdiv[k,s]*i)/(bdiv[k,s]*u - u + i)
  }}}
  
  ## Process model:
  for(k in 1:nrep){
    for(s in 1:nsite){
      ## Gamma diversity:
      gdiv[k,s] ~ dpois(lambda_gdiv[k])
      ## Beta diversity:
      bdiv[k,s] ~ dgamma(shape[k], rate[k])
    }
    shape[k] <- max(0.00001, mu_bdiv[k]^2/sigma_bdiv[k]^2)
    rate[k] <- max(0.00001, mu_bdiv[k]/sigma_bdiv[k]^2)
  }
  
  ## Priors:
  
  for(k in 1:nrep){
    lambda_gdiv[k] ~ dgamma(0.001, 0.001)
    mu_bdiv[k] ~ dunif(1.000001, 10)
    sigma_bdiv[k] ~ dgamma(0.001, 0.001)
  }
  
  ## Posterior calculations:
  
  ## Among replications comparison and comparison of gdiv estimates with 
  ## vegan::specpool():
  for(s in 1:nsite){
    gdiv_mean[s] <- mean(gdiv[,s])
    gdiv_se[s] <- sd(gdiv[,s])/sqrt(nrep)
    bdiv_mean[s] <- mean(bdiv[,s])
    for(k in 1:nrep){
      gdiv_diff[k,s] <- gdiv[k,s] - gdiv_mean[s]
      bdiv_diff[k,s] <- bdiv[k,s] - bdiv_mean[s]
  }}
  
  ## Alpha diversity:
  for(k in 1:nrep){
    for(s in 1:nsite){
      adiv[k,s] <- gdiv[k,s]/bdiv[k,s]
  }}
  
}

