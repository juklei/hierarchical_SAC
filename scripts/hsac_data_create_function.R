## This function creates species accumulation data, to which species 
## accumulation curves can be fitted.

## Author: Julian Klein
## Created: 20210122

## Input:

## data => A list of data frames. One data frame for every sampling plot, for 
##         which alpha, beta, and gamma diversity should be estimated. Each data
##         frame consists consists of the species identity as the rows and sample 
##         number as columns. The values in each data frame is the occurrence 
##         (0 or 1) of a certain species in a certain sample on a certain plot.
##         The format is the same as required for vegan::specpool().

## n_permutation => How many random orders of the samples, e.g. how many species
##                  acccumulation curves should be created?

## seed => For reproducability, create always same random collection of species
##         accumulation curves.

## Ouput:

## An array with dimensions sample*permutation*plot including plot names.

sac.create <- function(data, n_permutation, seed){

  ## Turn the elements into matrices:
  data <- lapply(data, function(x) as.matrix(x))

  ## Delete all species columns which have never been seen on a plot:
  data <- lapply(data, function(x) x[, colSums(x) != 0])

  ## Extract the maximum number of samples per plot to make sure all matrices 
  ## that will be part of the output array have the same dimension:
  max_sample <- max(sapply(lo_list, nrow))
  
  ## Define the function:
  sac.create.per.plot <- function(x) {
    out1 <- matrix(NA, max_sample, n_permutation)
    ## Create matrix with one random order of samples for every n_permutation:
    set.seed(seed)
    selection <- t(apply(matrix(nrow(x), nrow = n_permutation), 1, sample))
    ## Calculate number of unique species added, from the first to the last 
    ## sample in selection[i = 1, ] to selection[i = n_permutation, ]:
    for(i in 1:n_permutation){
      D <- x[selection[i,], ]
      ## For each species j in D, j's occurrence becomes 1 in all samples, 
      ## starting  from the first sample with occurence of j onward:  
      for(j in 1:ncol(D)) {D[min(which(D[, j] == 1)):nrow(D), j] <- 1}
      ## The row sums accross all species, then gives the accumulayed species 
      ## for every sample: 
      out1[1:nrow(D), i] <- rowSums(D)
    }
    return(out1)
  }

  ## Apply the function to all plots:
  out2 <- lapply(data, sac.create.per.plot)
  out2 <- array(unlist(out2), 
                dim = c(max_sample, n_permutation, length(out2)),
                dimnames = list(NULL, NULL, names(out2)))
  return(out2)

}