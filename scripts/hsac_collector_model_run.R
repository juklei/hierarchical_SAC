## The species accumulation curve is a modified michaelis-menten curve
## We put explanatory variables on both beta and gamma diversity
##
## First edit: 20190605
## Last edit: 20201117
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

rm(list = ls())

library(boot)
library(rjags)
library(runjags)
library(coda)
library(magrittr)
library(reshape2)
library(data.table)
library(parallel)
library(dclone)

## 2. Define or source functions used in this script ---------------------------

dir.create("results")
dir.create("figures")

## Print all rows for mcmc outputs
options(max.print = 10E5)

## Backscale function
backscale <- function(pred_data, model_input_data) {
  
  pred_data*attr(model_input_data, 'scaled:scale') + 
    attr(model_input_data, 'scaled:center')
  
}

## 3. Load and explore data ----------------------------------------------------

dir("clean")

load("clean/species_accumulation_data.rda")
load("clean/sad_tree_part.rda")
str(sad)
str(sad_tree)

## 4. Prepare data and inits ---------------------------------------------------

## Calculate the number of trees by site:
ntree <- apply(sad, 3, function(x) sum(!is.na(x[,2])))

sad <- sad[,1:5,]

## Create model data set:
data <- list(nrep = dim(sad)[2],
             nsite = dim(sad)[3],
             ntree = ntree,
             obs_all = sad,
             obs_red = sad[, sample(1:dim(sad)[2], 1), ],
             u = 1, ## Evaluation unit of alpha diversity
             dec = scale(sad_tree$dec),
             spruce = scale(sad_tree$spruce),
             pine = scale(sad_tree$pine),
             tsp_2 = ifelse(sad_tree$nr_tsp == 2, 1, 0),
             tsp_3 = ifelse(sad_tree$nr_tsp == 3, 1, 0),
             tsp_4 = ifelse(sad_tree$nr_tsp == 4, 1, 0),
             dbh = scale(sad_tree$dbh)[, 1])

## Look at correlations:
cor_list <- data[7:13]
cor_list$tsp1 <- ifelse(sad_tree$nr_tsp == 1, 1, 0)
cor_list$dec_quad <- data$dec^2
cor_list$pine_quad <- data$pine^2
cor_list$spruce_quad <- data$spruce^2

write.csv(cor(as.data.frame(cor_list)), "results/cor_hsac.csv")

str(data)

## Load the models:
m.raw <- "scripts/JAGS/hsac_collector_raw.R"
m.perc <- "scripts/JAGS/hsac_collector_tsp_percentage.R"
m.tsp <- "scripts/JAGS/hsac_collector_tsp_number.R"

start <- Sys.time()

n.adapt <- 5000; n.iter <- 10000; samples <- 5000; n.thin <- 10

## 5a. Run m.raw for model testing and comparison ------------------------------

## Add chain inits:
inits <- list(list(gdiv = matrix(60,data$nrep, data$nsite), 
                   bdiv = matrix(6,data$nrep, data$nsite),
                   lambda_gdiv = rep(37, data$nrep),
                   mu_bdiv = rep(4, data$nrep),
                   sigma_bdiv = rep(1.1, data$nrep)),
              list(gdiv = matrix(40,data$nrep, data$nsite), 
                   bdiv = matrix(2.1,data$nrep, data$nsite),
                   lambda_gdiv = rep(33, data$nrep),
                   mu_bdiv = rep(5, data$nrep),
                   sigma_bdiv = rep(3, data$nrep)),
              list(gdiv = matrix(80,data$nrep, data$nsite), 
                   bdiv = matrix(9,data$nrep, data$nsite),
                   lambda_gdiv = rep(15, data$nrep),
                   mu_bdiv = rep(3.1, data$nrep),
                   sigma_bdiv = rep(2, data$nrep)))

## Parallel computing:
cl <- makePSOCKcluster(3) ## On 3 cores

## Load model:
parJagsModel(cl = cl,
             name = "hsac",
             file = m.raw,
             data = data,
             inits = inits,
             n.chains = 3,
             n.adapt = n.adapt)

## Update:
parUpdate(cl = cl, object = "hsac", n.iter = n.iter)

## Extract samples from estimated parameters:
zc <- parCodaSamples(cl = cl, model = "hsac",
                     variable.names = c("lambda_gdiv", "mu_bdiv", "sigma_bdiv"),
                     n.iter = samples,
                     thin = n.thin)

## Export parameter estimates:
capture.output(summary(zc), HPDinterval(zc, prob = 0.95)) %>%
  write(., "results/parameters_hsac_collector_raw.txt")

## Test model convergence and mixing:
pdf("figures/hsac_collector_raw.pdf")
plot(zc); gelman.plot(zc)
dev.off()

## And with a table:
capture.output(raftery.diag(zc), heidel.diag(zc)) %>%
  write(., "results/diagnostics_hsac_collector.txt")

## Produce SAC from fitted model and compare with raw accumulation data:
zc_val <- parCodaSamples(cl = cl, model = "hsac",
                         variable.names = "sim_obs",
                         n.iter = samples,
                         thin = n.thin)
zc_val_comb <- combine.mcmc(zc_val)

pred <- as.data.table(melt(as.data.table(zc_val_comb), 
                           measure.vars = dimnames(zc_val_comb)[[2]]))
pred$variable <- gsub(",[[:digit:]]+,", ",", pred$variable)
quant_calc <- function(x) as.list(quantile(x, probs = c(0.025, 0.5, 0.975)))
pred <- pred[, c("lower", "median", "upper") := quant_calc(value), by = "variable"]
pred <- unique(pred[ , -2])
pred$tree <- unlist(lapply(data$ntree, function(x) rep(1:x)))
pred$site <- unlist(lapply(1:data$nsite, function(x) rep(x, ntree[x])))
pred <- as.data.frame(pred)

dev.off()

pdf("figures/hsac_collector_sim_vs_obs.pdf")
par(mfrow = c(3, 1))
for(i in 1:data$nsite) {
  x = c(0, pred[pred$site == i, "tree"])
  y = pred[pred$site == i, ]
  plot(x, c(0, y$upper),
       xlab = "tree number", ylab = "accumulation",
       lty = "dashed", col = "red", typ = "l")
  lines(x, c(0, y$median), col = "red", lwd = 4)
  lines(x, c(0, y$lower), lty = "dashed", col = "red")
  ## Real data:
  boxplot(accumulation ~ tree, 
          data.frame("tree" = as.factor(rep(1:(length(x)-1), dim(sad)[2])),
                     "accumulation" = na.omit(as.vector(sad[,,i]))),
          add = TRUE)
}
dev.off()

## How much do the estimates vary among replicated accumulation curves per plot,
## and does that vary with sample size (number of trees)?

zc_diff <- parCodaSamples(cl = cl, model = "hsac",
                          variable.names = c("gdiv_diff", "bdiv_diff"),
                          n.iter = samples,
                          thin = n.thin)
zc_diff_comb <- combine.mcmc(zc_diff)

## Create array with dim = mcmc_sample*replication*plot*diversity_metric:
diff <- array(zc_diff_comb, c(dim(zc_diff_comb)[1], dim(sad)[2:3], 2))
diff <- apply(diff, c(3,4), sd)

## Differences for bdiv:
plot(data$ntree, diff[,1])

## Differences for gdiv:
plot(data$ntree, diff[,2])

## Compare Gamma diversity from hsac_collector model with estimates from 
## vegan::specpool():

zc_comp <- parCodaSamples(cl = cl, model = "hsac",
                          variable.names = c("gdiv_mean", "gdiv_se"),
                          n.iter = samples,
                          thin = n.thin)




## 5b. Use m.raw to export diversity metrics for all sites ---------------------

zc_pred <- parCodaSamples(cl = cl, model = "hsac",
                          variable.names = c("adiv", "bdiv", "gdiv"),
                          n.iter = samples,
                          thin = n.thin)

site_estimates <- zj_pred
site_estimates$dec <- sad_tree$dec
site_estimates$spruce <- sad_tree$spruce
site_estimates$pine <- sad_tree$pine
save(site_estimates, file = "clean/site_estimates.rda")

stopCluster(cl)

## 6. Run m.dec ----------------------------------------------------------------

inits <- list(list(gdiv = rep(60, data$nsite), bdiv = rep(6, data$nsite),
                   g_icpt = log(37), g_dbh = 1, g_perc = 6, g_perc2 = -1.5,
                   sigma_bdiv = 0.9, 
                   b_icpt = 4, b_dbh = -0.3, b_perc = 0.5, b_perc2 = -0.25),
              list(gdiv = rep(40, data$nsite), bdiv = rep(2, data$nsite),
                   g_icpt = log(33), g_dbh = -2, g_perc = 0, g_perc2 = -4,
                   sigma_bdiv = 0.6, 
                   b_icpt = 3, b_dbh = -0.8, b_perc = 1.5, b_perc2 = -0.5),
              list(gdiv = rep(20, data$nsite), bdiv = rep(3, data$nsite),
                   g_icpt = log(41), g_dbh = 4, g_perc = 10, g_perc2 = 0.1,
                   sigma_bdiv = 1.2, 
                   b_icpt = 5, b_dbh = 0.2, b_perc = 0, b_perc2 = 0.1))

for(i in c("dec", "spruce", "pine")){

data$perc <- unlist(data[i])
data$perc_pred <- seq(min(data$perc), max(data$perc), 0.05)

cl <- makePSOCKcluster(3) 
parJagsModel(cl, "hsac", m.perc, data, inits, 3, n.adapt)
parUpdate(cl = cl, object = "hsac", n.iter = n.iter)

zc <- parCodaSamples(cl = cl, model = "hsac",
                     variable.names = c("g_icpt", "g_perc", "g_perc2", "g_dbh",
                                        "sigma_bdiv",
                                        "b_icpt", "b_perc", "b_perc2", "b_dbh"),
                     n.iter = samples, thin = n.thin)

capture.output(summary(zc), HPDinterval(zc, prob = 0.95)) %>% 
  write(., paste0("results/parameters_hsac_collector_", i, ".txt"))

pdf(paste0("figures/hsac_collector_", i, ".pdf"))
plot(zc); gelman.plot(zc) 
dev.off()

capture.output(raftery.diag(zc), heidel.diag(zc)) %>% 
  write(., paste0("results/diagnostics_hsac_collector_", i, ".txt"))

zj_pred <- parCodaSamples(cl = cl, model = "hsac",
                          variable.names = c("adiv_pred", "bdiv_pred", "gdiv_pred", 
                                             "b_perc_max", "g_perc_max"),
                          n.iter = samples, thin = n.thin)

ld <- length(data$perc_pred)
export <- list("adiv" = zj_pred[, 1:ld],
               "bdiv" = zj_pred[, (ld+2):(2*ld+1)],
               "gdiv" = zj_pred[, (2*ld+3):(3*ld+2)],
               "b_max" = backscale(unlist(zj_pred[, "b_perc_max"]), data[[i]]),
               "g_max" = backscale(unlist(zj_pred[, "g_perc_max"]), data[[i]]),
               "pred" = backscale(data$perc_pred, data[[i]]))
names(export) <- paste0(names(export), "_", i)

save(export, file = paste0("clean/export_collector_", i, ".rda"))

stopCluster(cl)

}

## 7. Run m.tsp ----------------------------------------------------------------

inits <- list(list(gdiv = rep(60, data$nsite), bdiv = rep(6, data$nsite),
                   g_icpt = log(37), g_dbh = 1, 
                   g_2tsp = 1, g_3tsp = 1, g_4tsp = 1,
                   sigma_bdiv = 0.9, 
                   b_icpt = 4, b_dbh = -0.3, 
                   b_2tsp = 1, b_3tsp = 1, b_4tsp = 1),
              list(gdiv = rep(40, data$nsite), bdiv = rep(2, data$nsite),
                   g_icpt = log(33), g_dbh = -2, 
                   g_2tsp = -1, g_3tsp = -1, g_4tsp = -1,
                   sigma_bdiv = 0.6, 
                   b_icpt = 3, b_dbh = -0.8, 
                   b_2tsp = -1, b_3tsp = -1, b_4tsp = -1),
              list(gdiv = rep(20, data$nsite), bdiv = rep(3, data$nsite),
                   g_icpt = log(41), g_dbh = 4, 
                   g_2tsp = 0, g_3tsp = 0, g_4tsp = 0,
                   sigma_bdiv = 1.2, 
                   b_icpt = 5, b_dbh = 0.2, 
                   b_2tsp = 0, b_3tsp = 0, b_4tsp = 0))

cl <- makePSOCKcluster(3) 
parJagsModel(cl, "hsac", m.tsp, data, inits, 3, n.adapt)
parUpdate(cl = cl, object = "hsac", n.iter = n.iter)

zc <- parCodaSamples(cl = cl, model = "hsac",
                     variable.names = c("g_icpt", "g_dbh",
                                        "g_2tsp", "g_3tsp", "g_4tsp", 
                                        "sigma_bdiv", "b_icpt", "b_dbh",
                                        "b_2tsp", "b_3tsp", "b_4tsp"),
                     n.iter = samples, thin = n.thin)

capture.output(summary(zc), HPDinterval(zc, prob = 0.95)) %>% 
  write(., "results/parameters_hsac_collector_tsp.txt")

pdf("figures/hsac_collector_tsp.pdf")
plot(zc); gelman.plot(zc) 
dev.off()

capture.output(raftery.diag(zc), heidel.diag(zc)) %>% 
  write(., "results/diagnostics_hsac_collector_tsp.txt")

zj_pred <- parCodaSamples(cl = cl, model = "hsac",
                          variable.names = c("adiv_1tsp", "adiv_2tsp", "adiv_3tsp", "adiv_4tsp", 
                                             "adiv_diff_21", "adiv_diff_31", 
                                             "adiv_diff_41", "adiv_diff_32",
                                             "adiv_diff_42", "adiv_diff_43",
                                             "bdiv_1tsp", "bdiv_2tsp", "bdiv_3tsp", "bdiv_4tsp", 
                                             "bdiv_diff_21", "bdiv_diff_31", 
                                             "bdiv_diff_41", "bdiv_diff_32",
                                             "bdiv_diff_42", "bdiv_diff_43",
                                             "gdiv_1tsp", "gdiv_2tsp", "gdiv_3tsp", "gdiv_4tsp", 
                                             "gdiv_diff_21", "gdiv_diff_31", 
                                             "gdiv_diff_41", "gdiv_diff_32",
                                             "gdiv_diff_42", "gdiv_diff_43"),
                          n.iter = samples, thin = n.thin)

## Extract probability that difference between nr_tsp is bigger than 0:
tsp_diff <- combine.mcmc(zj_pred)[, c(5:10, 15:20, 25:30)]
ANOVA_prob <- summary(tsp_diff)$quantiles[, c("50%", "2.5%", "97.5%")]
ANOVA_prob <- cbind(ANOVA_prob, 
                    "ecdf" = sapply(as.data.frame(tsp_diff), 
                                    function(x) 1-ecdf(x)(0))) 
write.csv(ANOVA_prob, "results/ANOVA_results_hsac_collector.csv")

export <- zj_pred[, c("adiv_1tsp", "adiv_2tsp", "adiv_3tsp", "adiv_4tsp",
                      "bdiv_1tsp", "bdiv_2tsp", "bdiv_3tsp", "bdiv_4tsp",
                      "gdiv_1tsp", "gdiv_2tsp", "gdiv_3tsp", "gdiv_4tsp")]
save(export, file = "clean/export_tsp_collector.rda")

stopCluster(cl)

end <- Sys.time()
end-start

## -------------------------------END-------------------------------------------
