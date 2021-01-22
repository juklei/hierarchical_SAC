## Create species accumulation data:
##
## First edit: 20201116
## Last edit: 20201119
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

rm(list = ls())

library(data.table)
library(vegan)

## Define seed to get the same sample in random processes (Reproducability):
seed <- 999

source("scripts/hsac_data_create_function.r")

## 2. Load and explore data ----------------------------------------------------

## From 10.5281/zenodo.3899847:
l_obs <- as.data.table(read.csv("data/epiphytic_lichen_data.csv"))

## Remove plot_119, because spruce plantation:
l_obs <- l_obs[l_obs$Plot.no. != 119, ]

## 3. Create list of matrices per plot that goes into the shuffeling ----------- 

## Change names of plot and circles for consistency:
colnames(l_obs)[2] <- "plot"
l_obs$plot <- as.factor(paste0("plot_", l_obs$plot))

## Reduce l_obs to only uncut trees:
l_obs <- l_obs[!is.na(l_obs$Tree.no), c(2:7, 10:14, 17:137)]

## Reduce to trees of known species:
l_obs <- droplevels(l_obs[l_obs$Tree.species %in% c("Ag", "Bp", "Pa", "Ps", "Pt", "Qr"), ])

## Merge Physcia adscendens, P. tenella och P. adscendens/P. tenella to
## P. adscendens/P. tenella"
T1 <- l_obs[, c("Physcia.adscendens", 
                "Physcia.tenella",
                "P..adscendens.tenella")]
l_obs$P..adscendens.tenella <- ifelse(rowSums(T1, na.rm = TRUE) > 0, 1, NA)

## Reduce data set to response columns:
lo_red <- l_obs[, c(1, 7, 12:24, 26:105, 107, 109:127, 129:132)]

## Sum species obs to total per tree and lichen species:
lo_red <- lo_red[, lapply(.SD, sum, na.rm = TRUE), by = c("plot", "Tree.no")]
lo_red[, 3:ncol(lo_red)][lo_red[, 3:ncol(lo_red)] > 1] <- 1

## 4. Make a table with species observation numbers structured by tree species:

# Get the host tree species identity back into the data set.
lo_freq <- merge(lo_red, l_obs[, c(1, 7:8)], by = c("plot", "Tree.no"))

## Calculate % of trees on which lichen species X was found per tree species:
lo_freq <- lo_freq[, lapply(.SD[, c(-1, -2)], function(x) sum(x)/nrow(.SD)*100), 
                   by = "Tree.species"]

## And the same for all trees:
lo_freq <- rbind(lo_freq, 
                 data.frame("Tree.species" = "All", 
                            t(colSums(lo_red[, c(-1, -2)])/nrow(lo_red)*100)))

## The total number of lichen species observed per tree species:
lo_freq[, "Total" := sum(.SD > 0), by = "Tree.species"]

## Prepare for output:
lo_freq_out <- as.data.frame(t(lo_freq[, -1]))
colnames(lo_freq_out) <- as.character(lo_freq$Tree.species)

## 5. Calculate estimated gamma diversity for all plots based on 
##    vegan::specpool():

lo_spec <- lo_red[, - 2] ## Remove Tree.no!
gamma_est <- lo_spec[, specpool(.SD), by = "plot"]

## 6. Create the data set which is needed to fit the accumultation curve and 
##    produce forest data per plot level with the according order

## Split into a list of data frames per plot:
lo_list <- split(lo_spec, by = "plot", keep.by = FALSE)

## Create the species accumulation data set for all plots:
sad <- sac.create(data = lo_list, n_permutation = 10, seed = seed)

## Store order of the list elements for ordering forest data:
plot_order <- dimnames(sad[[3]])

## Make a data frame with the tree data with the same plot order as "sad":
sad_tree <- as.data.table(l_obs[, c(1, 8:9)])
sad_tree <- sad_tree[, list("pine" = sum(Tree.species == "Ps")/nrow(.SD),
                            "spruce" = sum(Tree.species == "Pa")/nrow(.SD),
                            "dbh" = mean(Tree.diameter.130.cm.above.ground,
                                         na.rm = TRUE),
                            "nr_tsp" = length(unique(Tree.species))),
                     by = "plot"]
sad_tree$dec <- 1-(sad_tree$pine + sad_tree$spruce)

## Check if the order of l_sac and sad_tree are the same:
all(sad_tree$plot == names(sad))

## 7. Store the data used in the jags models and in the figures script: --------

dir.create("clean")

## Export:
write.csv(lo_freq_out, "results/obs_frequencies.csv")
write.csv(gamma_est, "clean/specpool_gamma.csv", row.names = FALSE)
save(sad, file = "clean/species_accumulation_data.rda")
save(sad_tree, file = "clean/sad_tree_part.rda")

## -------------------------------END-------------------------------------------
