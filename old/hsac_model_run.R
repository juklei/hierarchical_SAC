## model lichen richness per stand with species accumulation curves
## The species accumulation curve is a michaelis-menten saturation curve
## We put explanatory variables on both the half saturation and the asymptote 
##
## First edit: 20190605
## Last edit: 20201113
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

obs <- apply(sad, c(1,3), mean)

## Create model data set:
data <- list(nrep = dim(sad)[2],
             nsite = dim(sad)[3],
             ntree = ntree,
             obs = obs,
             u = 1, ## Evaluation unit of alpha diversity
             dec = scale(sad_tree$dec),
             spruce = scale(sad_tree$spruce),
             pine = scale(sad_tree$pine),
             tsp_2 = ifelse(sad_tree$nr_tsp == 2, 1, 0),
             tsp_3 = ifelse(sad_tree$nr_tsp == 3, 1, 0),
             tsp_4 = ifelse(sad_tree$nr_tsp == 4, 1, 0),
             dbh = scale(sad_tree$dbh)[, 1])

## Look at correlations:
cor_list <- data[6:12]
cor_list$tsp1 <- ifelse(sad_tree$nr_tsp == 1, 1, 0)
cor_list$dec_quad <- data$dec^2
cor_list$pine_quad <- data$pine^2
cor_list$spruce_quad <- data$spruce^2

capture.output(cor(as.data.frame(cor_list))) %>% 
  write(., "results/cor_hsac.txt")

str(data)

## Add chain inits:
inits <- list(list(rep_select = data$nrep/2,
                   sigma_obs = rep(1, data$nsite),
                   sigma_gdiv = 8, g_icpt = 37, g_dbh = 1, 
                   g_perc = 6, g_perc2 = -1.5,
                   g_2tsp = 1, g_3tsp = 1, g_4tsp = 1,
                   sigma_bdiv = 0.9, b_icpt = 4, b_dbh = -0.3, 
                   b_perc = 0.5, b_perc2 = -0.25,
                   b_2tsp = 1, b_3tsp = 1, b_4tsp = 1),
              list(rep_select = data$nrep/2,
                   sigma_obs = rep(1, data$nsite),
                   sigma_gdiv = 6, g_icpt = 33, g_dbh = -2, 
                   g_perc = 0, g_perc2 = -4,
                   g_2tsp = -1, g_3tsp = -1, g_4tsp = -1,
                   sigma_bdiv = 0.6, b_icpt = 3, b_dbh = -0.8, 
                   b_perc = 1.5, b_perc2 = -0.5,
                   b_2tsp = -1, b_3tsp = -1, b_4tsp = -1),
              list(rep_select = data$nrep/2,
                   sigma_obs = rep(1, data$nsite),
                   sigma_gdiv = 11, g_icpt = 41, g_dbh = 4, 
                   g_perc = 10, g_perc2 = 0.1,
                   g_2tsp = 0, g_3tsp = 0, g_4tsp = 0,
                   sigma_bdiv = 1.2, b_icpt = 5, b_dbh = 0.2, 
                   b_perc = 0, b_perc2 = 0.1,
                   b_2tsp = 0, b_3tsp = 0, b_4tsp = 0))

## Load the models:
m.raw <- "scripts/JAGS/hsac_tsp_percentage.R"
m.perc <- "scripts/JAGS/hsac_tsp_percentage.R"
m.tsp <- "scripts/JAGS/hsac_tsp_number.R"

start <- Sys.time()

n.adapt <- 1000; n.iter <- 1000; samples <- 1000; n.thin <- 2

## 5. Run m.raw ----------------------------------------------------------------

# MODEL TESTING; SIMULATIONS; AND PLOT LEVEL CALCULATIONS...
# ## Parallel computing:
# cl <- makePSOCKcluster(3) ## On 3 cores
# 
# ## Load model:
# parJagsModel(cl = cl, 
#              name = "lpsac", 
#              file = m.raw,
#              data = data,
#              inits = inits,
#              n.chains = 3,
#              n.adapt = n.adapt) 
# 
# ## Update:
# parUpdate(cl = cl, object = "lpsac", n.iter = n.iter)
# 
# ## Extract samples from estimated parameters:
# zc <- parCodaSamples(cl = cl, model = "lpsac",
#                      variable.names = c("mu_rich", "sigma_rich",
#                                         "mu_sat", "sigma_sat"),
#                      n.iter = samples, 
#                      thin = n.thin)
# 
# ## Export parameter estimates:
# capture.output(summary(zc), HPDinterval(zc, prob = 0.95)) %>% 
#   write(., "results/parameters_lpsac_raw.txt")
# 
# ## Test model convergence and mixing:
# pdf("figures/lpsac_raw.pdf")
# plot(zc); gelman.plot(zc) 
# dev.off()
# 
# ## And with a table:
# capture.output(raftery.diag(zc), heidel.diag(zc)) %>% 
#   write(., "results/diagnostics_lpsac_raw.txt")
# 
# ## Produce SAC from fitted model and compare with raw data:
# 
# zc_val <- parCodaSamples(cl = cl, model = "lpsac",
#                          variable.names = "obs_pred",
#                          n.iter = samples,
#                          thin = n.thin)
# 
# pred <- summary(zc_val, quantiles = c(0.025, 0.5, 0.975))$quantiles
# pred <- cbind(pred, "ntree" = unlist(lapply(data$ntree, function(x) 1:x)))
# plot <- list() ; for(i in 1:data$nplot){plot[[i]] <- rep(i, ntree[i])}
# pred <- as.data.frame(cbind(pred, "plot" = unlist(plot)))
# 
# dev.off()
# 
# pdf("figures/sim_vs_obs_new.pdf")
# par(mfrow = c(3, 2))
# for(i in 1:data$nplot) {
#   x = c(0, pred[pred$plot == i, "ntree"])
#   y = pred[pred$plot == i, ]
#   plot(x, c(0, y$`97.5%`), 
#        xlab = "tree nr", ylab = "richness",
#        lty = "dashed", col = "blue", typ = "l")
#   lines(x, c(0, y$`50%`), col = "blue")
#   lines(x, c(0, y$`2.5%`), lty = "dashed", col = "blue")
#   ## Real data:
#   points(1:(length(x)-1), na.omit(as.vector(obs[,i])))
# }
# dev.off()
# 
# ## Export estimated plot_richness and sat_speed for ternary plots:
# 
# zj_pred <- parCodaSamples(cl = cl, model = "lpsac",
#                           variable.names = c("plot_richness", "sat_speed"),
#                           n.iter = samples, 
#                           thin = n.thin)
# 
# plot_estimates <- zj_pred
# plot_estimates$dec <- sad_tree$dec
# plot_estimates$spruce <- sad_tree$spruce
# plot_estimates$pine <- sad_tree$pine
# save(plot_estimates, file = "clean/plot_estimates.rda")
# 
# stopCluster(cl)

## 6. Run m.dec ----------------------------------------------------------------

for(i in c("dec", "spruce", "pine")){

data$perc <- unlist(data[i])
data$perc_pred <- seq(min(data$perc), max(data$perc), 0.05)

cl <- makePSOCKcluster(3) 
parJagsModel(cl, "hsac", m.perc, data, inits, 3, n.adapt)
parUpdate(cl = cl, object = "hsac", n.iter = n.iter)

zc <- parCodaSamples(cl = cl, model = "hsac",
                     variable.names = c("sigma_obs",
                                        "sigma_gdiv",
                                        "g_icpt", "g_perc", "g_perc2", "g_dbh",
                                        "sigma_bdiv",
                                        "b_icpt", "b_perc", "b_perc2", "b_dbh"),
                     n.iter = samples, thin = n.thin)

capture.output(summary(zc), HPDinterval(zc, prob = 0.95)) %>% 
  write(., paste0("results/parameters_hsac_", i, ".txt"))

pdf(paste0("figures/hsac_", i, ".pdf"))
plot(zc); gelman.plot(zc) 
dev.off()

capture.output(raftery.diag(zc), heidel.diag(zc)) %>% 
  write(., paste0("results/diagnostics_hsac_", i, ".txt"))

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

save(export, file = paste0("clean/export_", i, ".rda"))

stopCluster(cl)

}

## 7. Run m.tsp ----------------------------------------------------------------

cl <- makePSOCKcluster(3) 
parJagsModel(cl, "hsac", m.tsp, data, inits, 3, n.adapt)
parUpdate(cl = cl, object = "hsac", n.iter = n.iter)

zc <- parCodaSamples(cl = cl, model = "hsac",
                     variable.names = c("sigma_obs", 
                                        "sigma_gdiv", "g_icpt", "g_dbh",
                                        "g_2tsp", "g_3tsp", "g_4tsp", 
                                        "sigma_bdiv", "b_icpt", "b_dbh",
                                        "b_2tsp", "b_3tsp", "b_4tsp"),
                     n.iter = samples, thin = n.thin)

capture.output(summary(zc), HPDinterval(zc, prob = 0.95)) %>% 
  write(., "results/parameters_hsac_tsp.txt")

pdf("figures/hsac_tsp.pdf")
plot(zc); gelman.plot(zc) 
dev.off()

capture.output(raftery.diag(zc), heidel.diag(zc)) %>% 
  write(., "results/diagnostics_hsac_tsp.txt")

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
write.csv(ANOVA_prob, "results/ANOVA_results_hsac.csv")

export <- zj_pred[, c("adiv_1tsp", "adiv_2tsp", "adiv_3tsp", "adiv_4tsp",
                      "bdiv_1tsp", "bdiv_2tsp", "bdiv_3tsp", "bdiv_4tsp",
                      "gdiv_1tsp", "gdiv_2tsp", "gdiv_3tsp", "gdiv_4tsp")]
save(export, file = "clean/export_tsp.rda")

stopCluster(cl)

end <- Sys.time()
end-start

## -------------------------------END-------------------------------------------
