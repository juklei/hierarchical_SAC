## Create species accumulation data:
##
## First edit: 20201116
## Last edit: 20201116
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

rm(list = ls())

library(data.table)

## 2. Load and explore data ----------------------------------------------------

l_obs <- as.data.table(read.csv("data/Data_lavar_Almunge15_March_2019.csv"))

## Remove plot_119, because spruce plantation:
l_obs <- l_obs[l_obs$Plot.no. != 119, ]



## 3. Create list of matrices per plot that goes into the shuffeling ----------- 

## Change names of plot and circles for consistency:
colnames(l_obs)[2] <- "plot"
l_obs$plot <- as.factor(paste0("plot_", l_obs$plot))

## Reduce l_obs to only uncut trees:
l_obs <- l_obs[!is.na(l_obs$Tree.no), c(2:7, 10:14, 17:137)]

## Reduce to trees of known species:
l_obs <- droplevels(l_obs[!is.na(l_obs$Tree.species) | 
                            l_obs$Tree.species != "check", ])

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




specaccum(lo_red[lo_red$plot == "plot_1", 3:length(lo_red)], "random")


## lo_red will be used below at 5. for species accumulation curves

## 4. Create the data set which is needed to fit the tree level richness models
##    and add tree data again from above

head(lo_red)

ltr <- lo_red[, list("richness" = sum(as.vector(.SD))), 
              by = c("plot", "Tree.no")]

ltr_tree <- merge(ltr, unique(l_obs[, c(1, 7:9)]), by = c("plot", "Tree.no"))

## 5. Create the data set which is needed to fit the accumultation curve and 
##    produce forest data per plot level with the according order

## Split into a list of data frames per plot:
lo_list <- split(lo_red, by = "plot", keep.by = FALSE)

## Turn the elements into matrices without the tree number:
lo_list <- lapply(lo_list, function(x) as.matrix(x[, -1]))

## Delete all species columns which have never been seen on a plot:
lo_list <- lapply(lo_list, function(x) x[, colSums(x) != 0])

## How many shuffelings?
S <- 100

## Define the function:
sac_create <- function(x) {
  
  out <- matrix(NA, max(l_obs$Tree.no), S)
  for(i in 1:S) {
    ## Select random rows in x:
    D <- x[sample(nrow(x)), ]
    for(j in 1:ncol(D)) {D[min(which(D[, j] == 1)):nrow(D), j] <- 1}
    out[1:nrow(D), i] <- rowSums(D)
  }
  return(out)

}

## Apply the function:
l_sac <- lapply(lo_list, sac_create)

## Store order of the list elements for ordering forest data:
plot_order <- names(l_sac)

## Turn l_sac into an array:
sad <- array(unlist(l_sac), dim = c(max(lengths(l_sac)/S), S, length(l_sac)))

## Make a data frame with the tree data with the same plot order as "sad":

sad_tree <- as.data.table(l_obs[, c(1, 8:9)])
sad_tree <- sad_tree[, list("pine" = sum(Tree.species == "Ps")/nrow(.SD),
                            "spruce" = sum(Tree.species == "Pa")/nrow(.SD),
                            "dbh" = mean(Tree.diameter.130.cm.above.ground,
                                         na.rm = TRUE),
                            "nr_tsp" = length(unique(Tree.species)),
                            "nr_dec" = sum(!unique(Tree.species) %in% c("Pa", 
                                                                        "Ps"))),
                     by = "plot"]
sad_tree$dec <- 1-(sad_tree$pine + sad_tree$spruce)

## Check if the order of l_sac and sad_tree are the same:
all(sad_tree$plot == names(l_sac))

## 6. Store the data sets used in the jags models here: ------------------------ 

dir.create("clean")

## Export:
save(ltr_tree, file = "clean/ltr_with_tree_data.rda")
save(sad, file = "clean/species_accumulation_data.rda")
save(sad_tree, file = "clean/sad_tree_part.rda")

## -------------------------------END-------------------------------------------
