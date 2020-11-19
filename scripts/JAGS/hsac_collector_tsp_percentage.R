## hsac_collector
##
## First edit: 20201115
## Last edit: 20201115
##
## Author: Julian Klein

model{

  ## Likelihood:
  
  ## Observation model:
  for(s in 1:nsite){
    for(i in 1:ntree[s]){
      obs_red[i,s] ~ dpois(lambda_obs[i,s])
      lambda_obs[i,s] <- (gdiv[s]*i)/(bdiv[s]*u - u + i)
  }}
  
  ## Process model:
  for(s in 1:nsite){
    ## Gamma diversity:
    gdiv[s] ~ dpois(lambda_gdiv[s])
    log(lambda_gdiv[s]) <- g_icpt + 
                           g_perc*perc[s] + 
                           g_perc2*perc[s]^2 + 
                           g_dbh*dbh[s]
    ## Beta diversity:
    bdiv[s] ~ dgamma(shape[s], rate[s])
    shape[s] <- max(0.00001, mu_bdiv[s]^2/sigma_bdiv^2)
    rate[s] <- max(0.00001, mu_bdiv[s]/sigma_bdiv^2)
    log(mu_bdiv[s]) <- b_icpt + b_perc*perc[s] + b_perc2*perc[s]^2 + b_dbh*dbh[s]
  }
  
  ## Priors:
  
  ## Process model:
  ## Gamma diversity:
  g_icpt ~ dgamma(0.001, 0.001)
  g_perc ~ dnorm(0, 0.001)
  g_perc2 ~ dnorm(0, 0.001)
  g_dbh ~ dnorm(0, 0.001)
  ## Beta diversity:
  sigma_bdiv ~ dgamma(0.001, 0.001)
  b_icpt ~ dgamma(0.001, 0.001)
  b_perc ~ dnorm(0, 0.001)
  b_perc2 ~ dnorm(0, 0.001)
  b_dbh ~ dnorm(0, 0.001) 
  
  ## Predictions:

  ## Diversity metrics predictions:
  for(m in 1:length(perc_pred)){
    ## Gamma diversity:
    log(gdiv_pred[m]) <- g_icpt + g_perc*perc_pred[m] + g_perc2*perc_pred[m]^2
    ## Beta diversity:
    log(bdiv_pred[m]) <- b_icpt + b_perc*perc_pred[m] + b_perc2*perc_pred[m]^2
    ## Alpha diversity:
    adiv_pred[m] <- gdiv_pred[m]/bdiv_pred[m]
  }

  ## Maxima:
  g_perc_max <- - g_perc/(2*g_perc2)
  b_perc_max <- - b_perc/(2*b_perc2)
  
}

