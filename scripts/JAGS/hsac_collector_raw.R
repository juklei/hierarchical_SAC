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
    shape[k] <- mu_bdiv[k]^2/sigma_bdiv^2
    rate[k] <- mu_bdiv[k]/sigma_bdiv^2
  }
  
  ## Priors:
  
  for(k in 1:nrep){
    lambda_gdiv[k] ~ dgamma(0.001, 0.001)
    log(mu_bdiv[k]) <- mu_bdiv_log[k] 
    mu_bdiv_log[k] ~ dgamma(0.001, 0.001)
  }
  sigma_bdiv ~ dt(0, pow(2.5,-2), 1)T(0,)
  
  ## Posterior calculations:
  
  ## Alpha diversity:
  for(k in 1:nrep){
    for(s in 1:nsite){
      adiv[k,s] <- gdiv[k,s]/bdiv[k,s]
  }}
  
  ## Among replications comparison:
  for(s in 1:nsite){
    gdiv_cv[s] <- sd(gdiv[,s])/mean(gdiv[,s])
    bdiv_cv[s] <- sd(bdiv[,s])/mean(bdiv[,s])
    adiv_cv[s] <- sd(adiv[,s])/mean(adiv[,s])
  }
  
  ## Extracting estimates for the permutation used in the covariate model parts
  ## and for comparison of gdiv estimates with vegan::specpool():
  for(s in 1:nsite){
    gdiv_sel[s] <- gdiv[select, s]
    bdiv_sel[s] <- bdiv[select, s]
    adiv_sel[s] <- adiv[select, s]
  }
  
}

