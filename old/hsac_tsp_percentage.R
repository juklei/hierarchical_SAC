## hsac model without covariates
##
## First edit: 20190605
## Last edit: 20201110
##
## Author: Julian Klein

model{
  
  ## Likelihood:
  
  ## Observation model:
  for(s in 1:nsite){
    for(i in 1:ntree[s]){
      obs[i,s] ~ dnorm(mu_obs[i,s], 1/sigma_obs[s]^2) T(0,)
      mu_obs[i,s] <- (gdiv[s]*i)/(bdiv[s]*u - u + i)
  }}
  
  ## Process model:
  for(s in 1:nsite){
    ## Gamma diversity:
    gdiv[s] ~ dnorm(mu_gdiv[s], 1/sigma_gdiv^2) T(1,)
    mu_gdiv[s] <- g_icpt + g_perc*perc[s] + g_perc2*perc[s]^2 + g_dbh*dbh[s]
    ## Beta diversity:
    bdiv[s] ~ dnorm(mu_bdiv[s], 1/sigma_bdiv^2) T(1,)
    mu_bdiv[s] <- b_icpt + b_perc*perc[s] + b_perc2*perc[s]^2 + b_dbh*dbh[s]
  }
  
  ## Priors:
  
  ## Observational model:
  for(s in 1:nsite){
    sigma_obs[s] ~ dgamma(0.001, 0.001)
  }
  
  ## Process model:
  ## Gamma diversity:
  sigma_gdiv ~ dgamma(0.001, 0.001)
  g_icpt ~ dnorm(1, 0.001) T(1,)
  g_perc ~ dnorm(0, 0.001)
  g_perc2 ~ dnorm(0, 0.001)
  g_dbh ~ dnorm(0, 0.001)
  ## Beta diversity:
  sigma_bdiv ~ dgamma(0.001, 0.001)
  b_icpt ~ dnorm(1, 0.001) T(1,)
  b_perc ~ dnorm(0, 0.001)
  b_perc2 ~ dnorm(0, 0.001)
  b_dbh ~ dnorm(0, 0.001)

  ## Predictions:

  ## Diversity metrics predictions:
  for(m in 1:length(perc_pred)){
    ## Gamma diversity:
    gdiv_pred[m] <- g_icpt + g_perc*perc_pred[m] + g_perc2*perc_pred[m]^2
    ## Beta diversity:
    bdiv_pred[m] <- b_icpt + b_perc*perc_pred[m] + b_perc2*perc_pred[m]^2
    ## Alpha diversity:
    adiv_pred[m] <- gdiv_pred[m]/bdiv_pred[m]
  }

  ## Maxima:
  g_perc_max <- - g_perc/(2*g_perc2)
  b_perc_max <- - b_perc/(2*b_perc2)

}

