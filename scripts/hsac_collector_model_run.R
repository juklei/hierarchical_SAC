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
library(data.table)
library(vegan)

## 2. Define or source functions used in this script ---------------------------

## Define seed to get the same sample in random processes (Reproducability):
seed <- 999

dir.create("results")
dir.create("figures")

## Print all rows for mcmc outputs
options(max.print = 10E5)

## Backscale function
backscale <- function(pred_data, model_input_data) {
  pred_data*attr(model_input_data, 'scaled:scale') + 
    attr(model_input_data, 'scaled:center')
}

## Calculate median and 95% CI for data.table:
quant_calc <- function(x) as.list(quantile(x, probs = c(0.025, 0.5, 0.975)))

## 3. Load and explore data ----------------------------------------------------

dir("clean")

load("clean/species_accumulation_data.rda")
load("clean/sad_tree_part.rda")
str(sad)
str(sad_tree)

## 4. Prepare data and inits ---------------------------------------------------

set.seed(seed)
select <- sample(1:dim(sad)[2], 1)

## Create model data set:
data <- list(nrep = dim(sad)[2],
             nsite = dim(sad)[3],
             ntree = apply(sad, 3, function(x) sum(!is.na(x[,2]))), 
             obs_all = sad,
             obs_red = sad[, select, ],
             select = select,
             u = 1, ## Evaluation unit of alpha diversity
             dec = scale(sad_tree$dec),
             spruce = scale(sad_tree$spruce),
             pine = scale(sad_tree$pine),
             tsp_2 = ifelse(sad_tree$nr_tsp == 2, 1, 0),
             tsp_3 = ifelse(sad_tree$nr_tsp == 3, 1, 0),
             tsp_4 = ifelse(sad_tree$nr_tsp == 4, 1, 0),
             dbh = scale(sad_tree$dbh)[, 1])

## Look at correlations:
cor_list <- data[8:14]
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

n.adapt <- 5000; n.iter <- 15000; samples <- 200000; n.thin <- 160

## 5a. Run m.raw for model testing and comparison ------------------------------

## Add chain inits:
inits <- list(list(gdiv = matrix(60,data$nrep, data$nsite), 
                   bdiv = matrix(6,data$nrep, data$nsite),
                   lambda_gdiv = rep(37, data$nrep),
                   mu_bdiv_log = rep(4, data$nrep),
                   sigma_bdiv = 1.1),
              list(gdiv = matrix(40,data$nrep, data$nsite), 
                   bdiv = matrix(2.1,data$nrep, data$nsite),
                   lambda_gdiv = rep(33, data$nrep),
                   mu_bdiv_log = rep(5, data$nrep),
                   sigma_bdiv = 3),
              list(gdiv = matrix(80,data$nrep, data$nsite), 
                   bdiv_log = matrix(9,data$nrep, data$nsite),
                   lambda_gdiv = rep(15, data$nrep),
                   mu_bdiv_log = rep(3.1, data$nrep),
                   sigma_bdiv = 2))

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
                         n.iter = samples/40,
                         thin = n.thin/40)

sim <- as.data.frame(summary(zc_val)$quantiles)
sim$ntree <- unlist(lapply(data$ntree, function(x) rep(1:x)))
sim$site <- unlist(lapply(1:data$nsite, function(x) rep(x, data$ntree[x])))

dev.off()

pdf("figures/hsac_collector_sim_vs_obs.pdf")
par(mfrow = c(3, 1))
for(i in 1:data$nsite) {
  x = c(0, sim[sim$site == i, "ntree"])
  y = sim[sim$site == i, ]
  plot(x, c(0, y$'97.5%'), xlab = "tree number", ylab = "species accumulation", 
       lty = "dashed", col = "red", typ = "l")
  lines(x, c(0, y$'50%'), col = "red", lwd = 4)
  lines(x, c(0, y$'2.5%'), lty = "dashed", col = "red")
  ## Real data:
  points(x = x[-1], y = na.omit(as.vector(sad[,data$select,i])))
}

dev.off()

## How much do the estimates vary among replicated accumulation curves per plot,
## and does that vary with sample size (number of trees)?

zc_cv <- parCodaSamples(cl = cl, model = "hsac",
                          variable.names = c("gdiv_cv", "bdiv_cv", "adiv_cv"),
                          n.iter = samples,
                          thin = n.thin)
cv <- as.data.frame(summary(zc_cv)$quantile)
cv$ntree <- data$ntree
cv$div_metric <- c("adiv", "bdiv", "gdiv")
cv$div_metric <- sort(cv$div_metric)

## Export for graphing:
write.csv(cv, "clean/hsac_collector_raw_replication_cv.csv")

## 5b. Use m.raw to export diversity metrics for all sites to make graphs ------
##     and compare gdiv with the estimates of vegan::specpool

zc_pred <- parCodaSamples(cl = cl, model = "hsac",
                          variable.names = c("adiv_sel", "bdiv_sel", "gdiv_sel"),
                          n.iter = samples,
                          thin = n.thin)
pred <- as.data.frame(summary(zc_pred)$quantiles)
pred$mean <- summary(zc_pred)$statistics[, "Mean"]
pred$div_metric <- gsub("_sel[[[:digit:]]+]", "", rownames(pred))
pred$site <- 1:data$nsite
pred <- cbind(pred, sad_tree[, c("dec", "spruce", "pine")])

## Export for graphing:
write.csv(pred, "clean/hsac_collector_raw_site_pred.csv")

stopCluster(cl)

## 6. Run m.perc ---------------------------------------------------------------

inits <- list(list(gdiv = rep(60, data$nsite), bdiv = rep(6, data$nsite),
                   g_icpt = 1, g_dbh = 1, g_perc = 1, g_perc2 = 1,
                   sigma_bdiv = 1, 
                   b_icpt = 1.1, b_dbh = 1, b_perc = 1, b_perc2 = 1),
              list(gdiv = rep(40, data$nsite), bdiv = rep(2, data$nsite),
                   g_icpt = 3.5, g_dbh = 0, g_perc = 0, g_perc2 = -0.1,
                   sigma_bdiv = 0.5, 
                   b_icpt = 1.15, b_dbh = -0.1, b_perc = 0, b_perc2 = -0.15),
              list(gdiv = rep(20, data$nsite), bdiv = rep(3, data$nsite),
                   g_icpt = 3.7, g_dbh = 0.15, g_perc = 0.25, g_perc2 = 0.01,
                   sigma_bdiv = 1.5, 
                   b_icpt = 1.5, b_dbh = 0.2, b_perc = 0.45, b_perc2 = 0.1))

perc_pred_export <- vector("list", 3)
names(perc_pred_export) <- c("dec", "spruce", "pine")
perc_max_export <- vector("list", 3)
names(perc_max_export) <- c("dec", "spruce", "pine")

for(i in c("dec", "spruce", "pine")){
  
  data$perc <- unlist(data[i])
  data$perc_pred <- seq(min(data$perc), max(data$perc), 0.05)

  cl <- makePSOCKcluster(3) 
  parJagsModel(cl, "hsac", m.perc, data, inits, 3, n.adapt)
  parUpdate(cl = cl, object = "hsac", n.iter = n.iter)

  zc <- parCodaSamples(cl = cl, model = "hsac",
                       variable.names = c("g_icpt", "g_dbh",
                                          "g_perc", "g_perc2", 
                                          "sigma_bdiv", "b_icpt", "b_dbh",
                                          "b_perc", "b_perc2"),
                       n.iter = samples*10, thin = n.thin*10)

  capture.output(summary(zc), HPDinterval(zc, prob = 0.95)) %>% 
    write(., paste0("results/parameters_hsac_collector_", i, ".txt"))

  pdf(paste0("figures/hsac_collector_", i, ".pdf"))
  plot(zc); gelman.plot(zc) 
  dev.off()

  capture.output(raftery.diag(zc), heidel.diag(zc)) %>% 
    write(., paste0("results/diagnostics_hsac_collector_", i, ".txt"))

  ## Export predicted diversity metrics:
  zc_pred_perc <- parCodaSamples(cl = cl, model = "hsac",
                                 variable.names = c("adiv_pred", 
                                                    "bdiv_pred",
                                                    "gdiv_pred"),
                                 n.iter = samples*10, thin = n.thin*10)

  pred_perc <- as.data.frame(summary(zc_pred_perc)$quantile)
  pred_perc$perc_pred <- backscale(data$perc_pred, data[[i]])
  pred_perc$div_metric <- c("adiv", "bdiv", "gdiv")
  pred_perc$div_metric <- sort(pred_perc$div_metric)
  
  perc_pred_export[[i]] <- pred_perc

  ## Export maximum percentage values for the quadratic curves of gdiv and bdiv:
  zc_pp_max <- parCodaSamples(cl = cl, model = "hsac",
                              variable.names = c("b_perc_max", "g_perc_max"),
                              n.iter = samples*10, thin = n.thin*10)

  pp_max <- data.frame("b_max" = backscale(unlist(zc_pp_max[, "b_perc_max"]), 
                                           data[[i]]),
                       "g_max" = backscale(unlist(zc_pp_max[, "g_perc_max"]), 
                                           data[[i]]))
  
  perc_max_export[[i]] <- pp_max

  stopCluster(cl)

}

## Export for graphing:
pred_perc <- reshape2::melt(perc_pred_export,
                            measure.vars = colnames(perc_pred_export))
write.csv(pred_perc, paste0("clean/hsac_collector_pred_perc.csv"))

## Export for graphing:
pp_max <- reshape2::melt(perc_max_export, 
                         measure.vars = colnames(perc_max_export))
write.csv(pp_max, "clean/hsac_collector_max_perc.csv")

## 7. Run m.tsp ----------------------------------------------------------------

inits <- list(list(gdiv = rep(60, data$nsite), bdiv = rep(6, data$nsite),
                   g_icpt = log(30), g_dbh = 1, 
                   g_2tsp = 1, g_3tsp = 1, g_4tsp = 1,
                   sigma_bdiv = 1, 
                   b_icpt = 1, b_dbh = 1, 
                   b_2tsp = 1, b_3tsp = 1, b_4tsp = 1),
              list(gdiv = rep(40, data$nsite), bdiv = rep(2, data$nsite),
                   g_icpt = log(20), g_dbh = 0, 
                   g_2tsp = 0, g_3tsp = 0.1, g_4tsp = 0.2,
                   sigma_bdiv = 0.1, 
                   b_icpt = 0.6, b_dbh = -0.08, 
                   b_2tsp = -0.4, b_3tsp = 0.8, b_4tsp = -0.2),
              list(gdiv = rep(20, data$nsite), bdiv = rep(3, data$nsite),
                   g_icpt = log(35), g_dbh = 0.15, 
                   g_2tsp = 0.6, g_3tsp = 0.6, g_4tsp = 0.8,
                   sigma_bdiv = 0.8, 
                   b_icpt = 1.4, b_dbh = 0.2, 
                   b_2tsp = 0.5, b_3tsp = 0, b_4tsp = 0.8))

cl <- makePSOCKcluster(3) 
parJagsModel(cl, "hsac", m.tsp, data, inits, 3, n.adapt)
parUpdate(cl = cl, object = "hsac", n.iter = n.iter*2)

zc <- parCodaSamples(cl = cl, model = "hsac",
                     variable.names = c("g_icpt", "g_dbh",
                                        "g_2tsp", "g_3tsp", "g_4tsp", 
                                        "sigma_bdiv", "b_icpt", "b_dbh",
                                        "b_2tsp", "b_3tsp", "b_4tsp"),
                     n.iter = samples*4, thin = n.thin*4)

capture.output(summary(zc), HPDinterval(zc, prob = 0.95)) %>% 
  write(., "results/parameters_hsac_collector_tsp.txt")

pdf("figures/hsac_collector_tsp.pdf")
plot(zc); gelman.plot(zc) 
dev.off()

capture.output(raftery.diag(zc), heidel.diag(zc)) %>% 
  write(., "results/diagnostics_hsac_collector_tsp.txt")

zc_pred_tsp <- parCodaSamples(cl = cl, model = "hsac",
                              variable.names = c("adiv_1tsp", "adiv_2tsp", 
                                                 "adiv_3tsp", "adiv_4tsp",
                                                 "bdiv_1tsp", "bdiv_2tsp", 
                                                 "bdiv_3tsp", "bdiv_4tsp",
                                                 "gdiv_1tsp", "gdiv_2tsp", 
                                                 "gdiv_3tsp", "gdiv_4tsp"),
                              n.iter = samples*4, thin = n.thin*4)

pred_tsp <- as.data.frame(summary(zc_pred_tsp)$quantiles)
pred_tsp$div_metric <- c("adiv", "bdiv", "gdiv")
pred_tsp$div_metric <- sort(pred_tsp$div_metric)
pred_tsp$nr_tsp <- rep(1:4, 3)

write.csv(pred_tsp, paste0("clean/hsac_collector_pred_tsp.csv"))

zc_tsp_diff <- parCodaSamples(cl = cl, model = "hsac",
                              variable.names = c("adiv_diff_21", "adiv_diff_31", 
                                                 "adiv_diff_41", "adiv_diff_32",
                                                 "adiv_diff_42", "adiv_diff_43",
                                                 "bdiv_diff_21", "bdiv_diff_31", 
                                                 "bdiv_diff_41", "bdiv_diff_32",
                                                 "bdiv_diff_42", "bdiv_diff_43",
                                                 "gdiv_diff_21", "gdiv_diff_31", 
                                                 "gdiv_diff_41", "gdiv_diff_32",
                                                 "gdiv_diff_42", "gdiv_diff_43"),
                              n.iter = samples*4, thin = n.thin*4)

## Extract probability that difference between nr_tsp is bigger than 0:
ANOVA_prob <- summary(zc_tsp_diff)$quantiles
ANOVA_prob <- cbind(ANOVA_prob, 
                    "ecdf" = sapply(as.data.frame(combine.mcmc(zc_tsp_diff)), 
                                    function(x) 1-ecdf(x)(0))) 
write.csv(ANOVA_prob, "results/hsac_collector_tsp_ANOVA.csv")

stopCluster(cl)

end <- Sys.time()
end-start

## -------------------------------END-------------------------------------------
