## hsac_raw
##
## First edit: 20201115
## Last edit: 20201115
##
## Author: Julian Klein

model{
  
  ## Likelihood:
  
  ## Observation model:
  for(s in 1:nsite){
    for(k in 1:nrep){
      for(i in 1:ntree[s]){
        obs[i,k,s] ~ dpois(lambda_obs[i,k,s])
        lambda_obs[i,k,s] <- (gdiv[k,s]*i)/(bdiv[k,s]*u - u + i)
  }}}
  
  ## Process model:
  for(s in 1:nsite){
    for(k in 1:nrep){
      gdiv[k,s] ~ dpois(lambda_gdiv[k])
      bdiv[k,s] ~ dnorm(mu_bdiv[k], 1/sigma_bdiv^2) T(1,)
  }}
  
  ## Priors:
  for(k in 1:nrep){
    lambda_gdiv[k] ~ dgamma(0.001, 0.001)
    mu_bdiv[k] ~ dnorm(1, 0.001) T(1,)
  }
  sigma_bdiv ~ dgamma(0.001, 0.001)
  
}

