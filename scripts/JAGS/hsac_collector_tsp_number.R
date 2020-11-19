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
                           g_2tsp*tsp_2[s] + 
                           g_3tsp*tsp_3[s] +
                           g_4tsp*tsp_4[s] +
                           g_dbh*dbh[s]
    ## Beta diversity:
    bdiv[s] ~ dgamma(shape[s], rate[s])
    shape[s] <- max(0.00001, mu_bdiv[s]^2/sigma_bdiv^2)
    rate[s] <- max(0.00001, mu_bdiv[s]/sigma_bdiv^2)
    log(mu_bdiv[s]) <- b_icpt +                  
                       b_2tsp*tsp_2[s] + 
                       b_3tsp*tsp_3[s] +
                       b_4tsp*tsp_4[s] + 
                       b_dbh*dbh[s]
  }
  
  ## Priors:
  
  ## Process model:
  ## Gamma diversity:
  g_icpt ~ dgamma(0.001, 0.001)
  g_2tsp ~ dnorm(0, 0.001)
  g_3tsp ~ dnorm(0, 0.001)
  g_4tsp ~ dnorm(0, 0.001)
  g_dbh ~ dnorm(0, 0.001)
  ## Beta diversity:
  sigma_bdiv ~ dgamma(0.001, 0.001)
  b_icpt ~ dgamma(0.001, 0.001)
  b_2tsp ~ dnorm(0, 0.001)
  b_3tsp ~ dnorm(0, 0.001)
  b_4tsp ~ dnorm(0, 0.001)
  b_dbh ~ dnorm(0, 0.001)
  
  ## Predictions:
  
  ## Gamma diversity:
  log(gdiv_1tsp) <- g_icpt
  log(gdiv_2tsp) <- g_icpt + g_2tsp
  log(gdiv_3tsp) <- g_icpt + g_3tsp
  log(gdiv_4tsp) <- g_icpt + g_4tsp
  gdiv_diff_21 <- gdiv_2tsp - gdiv_1tsp
  gdiv_diff_31 <- gdiv_3tsp - gdiv_1tsp
  gdiv_diff_41 <- gdiv_4tsp - gdiv_1tsp
  gdiv_diff_32 <- gdiv_3tsp - gdiv_2tsp
  gdiv_diff_42 <- gdiv_4tsp - gdiv_2tsp
  gdiv_diff_43 <- gdiv_4tsp - gdiv_3tsp
  ## Beta diversity:
  log(bdiv_1tsp) <- b_icpt
  log(bdiv_2tsp) <- b_icpt + b_2tsp
  log(bdiv_3tsp) <- b_icpt + b_3tsp
  log(bdiv_4tsp) <- b_icpt + b_4tsp
  bdiv_diff_21 <- bdiv_2tsp - bdiv_1tsp
  bdiv_diff_31 <- bdiv_3tsp - bdiv_1tsp
  bdiv_diff_41 <- bdiv_4tsp - bdiv_1tsp
  bdiv_diff_32 <- bdiv_3tsp - bdiv_2tsp
  bdiv_diff_42 <- bdiv_4tsp - bdiv_2tsp
  bdiv_diff_43 <- bdiv_4tsp - bdiv_3tsp
  ## Alpha diversity:
  adiv_1tsp <- gdiv_1tsp/bdiv_1tsp
  adiv_2tsp <- gdiv_2tsp/bdiv_2tsp
  adiv_3tsp <- gdiv_3tsp/bdiv_3tsp
  adiv_4tsp <- gdiv_4tsp/bdiv_4tsp
  adiv_diff_21 <- adiv_2tsp - adiv_1tsp
  adiv_diff_31 <- adiv_3tsp - adiv_1tsp
  adiv_diff_41 <- adiv_4tsp - adiv_1tsp
  adiv_diff_32 <- adiv_3tsp - adiv_2tsp
  adiv_diff_42 <- adiv_4tsp - adiv_2tsp
  adiv_diff_43 <- adiv_4tsp - adiv_3tsp
  
}

