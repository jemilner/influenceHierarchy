#setwd("U:/Documents/R-code/practice2/spat-hom-hierarchy")
setwd("C:/Users/Owner/Documents/PhD/R-code/practice2/spat-hom-hierarchy")

library(parallel)

rm(list = ls())

num.cores <- 2
# num.cores <- detectCores() - 5

main <- function(i)
{
    param.set <<- 1
    seed <<- i
    
    source('gfit-setup-sim-mass.R')
}

#setup cluster
clust <- makeCluster(num.cores)
#use parLapply to run function across the cluster
parLapply(clust, 1:4, main)

#free up cluster once finished
stopCluster(clust)