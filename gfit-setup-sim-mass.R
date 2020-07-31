#for the dirichlet distribution
library(MCMCpack)

source('gfit-move-calc.R')
source('slim-mvnorm.R')

source('gfit-sim-move.R')
source('gfit-move-like-slim.R')

source('gfit-update-rate-het.R')
source('gfit-state-tab.R')

source('gfit-update-params.R')
source('gfit-utilities.R')

##########################
####load observed data####
##########################
#dat will be part of the output filenames
set.seed(1107)
dat <- paste("set", param.set, "-seed", seed, sep = "")
    
load(paste("data/mass-testing/linear-set-", param.set, "-seed-", seed, "-data.RData", sep = ""))
raw.data <- linear.data

num.cols <- 7
col.x <- 1; col.y <- 2; col.time <- 3; col.id <- 4; col.type <- 5; col.state <- 6; col.jump <- 7
type.obs <- -1; type.partial <- 0

#order the data by animal id and time
raw.data <- raw.data[order(raw.data[, col.id],raw.data[, col.time]), ]

animal.id <- unique(raw.data[, col.id])
unique.times <- unique(raw.data[, col.time])
ord <- order(unique.times)
unique.times <- unique.times[ord]

total.animals <- length(animal.id)
num.pairs <- dim(combn(total.animals, 2))[2]
num.obs <- length(unique.times)

#set up data matrix
animals <- array(NA, c(num.obs, num.cols, total.animals))

for(i in 1:total.animals)
{
    which.obs <- which(raw.data[, col.id] == i)
    tmp.times <- raw.data[which.obs, col.time]
    #which.times puts animal is data in the correct rows in the animals matrix
    which.times <- which(unique.times %in% tmp.times)
    which.missing <- c(1:num.obs)[-which.times]
    
    animals[which.times, col.x, i] <- raw.data[which.obs, col.x]
    animals[which.times, col.y, i] <- raw.data[which.obs, col.y]
    animals[, col.time, i] <- unique.times
    animals[, col.id, i] <- i
    animals[which.times, col.type, i] <- type.obs
    animals[which.missing, col.type, i] <- type.partial
    animals[, col.jump, i] <- 0
    
    #sample locations to 'complete' partial observations
    num.missing <- length(which.missing)
    
    if(num.missing > 0)
    {
        #if first obs is partial, fill it with the next available data
        if(any(which.missing == 1))
        {
            next.obs <- which.times[1]
            animals[1, c(col.x, col.y), i] <- animals[next.obs, c(col.x, col.y), i]
        }
        
        #if last obs is partial, fill it by last  available obs
        if(any(which.missing == num.obs))
        {
            prev.obs <- which.times[length(which.times)]
            animals[num.obs, c(col.x, col.y), i] <- animals[prev.obs, c(col.x, col.y), i]
        }
        
        #for all others, use an average of the prev and next available obs
        for(j in 1:num.missing)
        {
            if(which.missing[j] != 1 & which.missing[j] != num.obs)
            {
                prev.obs <- which.missing[j] - 1 #will already by filled in, obs or not
                next.obs <- which(animals[, col.time, i] > animals[which.missing[j], col.time, i] &
                                      animals[, col.type, i] == type.obs)[1]
                
                dt <- (animals[which.missing[j], col.time, i] - animals[prev.obs, col.time, i]) /
                    (animals[next.obs, col.time, i] - animals[prev.obs, col.time, i])
                
                animals[which.missing[j], col.x, i] <- (1 - dt) * animals[prev.obs, col.x, i] +
                    dt * animals[next.obs, col.x, i]
                animals[which.missing[j], col.y, i] <- (1 - dt) * animals[prev.obs, col.y, i] +
                    dt * animals[next.obs, col.y, i]
            }
        }
    }
}

colnames(animals) <- c("X", "Y", "Time", "ID", "Obs", "State", "Jump")

rm(raw.data, linear.data)

#############################################
####set up MCMC initial values and priors####
#############################################
#number of states
num.bm <- 2
num.states <- total.animals + num.bm

follower.states <- c(1:total.animals)
bm.states <- c((total.animals + 1):(total.animals + num.bm))

#initialise model parameters
if(param.set == 1) {
    alpha <- 0.7
    sigma <- 1
    rho <- c(0.9, 1.1)
    
    prop.alpha <- 0.015
    prop.sigma <- 0.015
    prop.slow <- 0.015
    prop.fast <- 0.08
} else if(param.set == 2) {
    alpha <- 0.1
    sigma <- 1
    rho <- c(0.2, 1)
    
    prop.alpha <- 0.0025
    prop.sigma <- 0.015
    prop.slow <- 0.0025
    prop.fast <- 0.02
} else if(param.set == 3) {
    alpha <- 0.7
    sigma <- 1
    rho <- c(0.9, 1.1)
    
    prop.alpha <- 0.01
    prop.sigma <- 0.015
    prop.slow <- 0.015
    prop.fast <- 0.06 
} else if(param.set == 4) {
    alpha <- 0.1 #start weak to encourage false positives
    sigma <- 1
    rho <- c(0.9, 1.1)
    
    prop.alpha <- 0.025
    prop.sigma <- 0.05
    prop.slow <- 0.01
    prop.fast <- 0.05
}

#variance for if the first obs is partial
partial.var <- diag(1) * 0.1

#max rate for the switching state poisson process
max.rate <- 0.2

kappa <- max.rate * total.animals
#priors for transition rates
shape1 <- shape2 <- 1 / num.states

#initial lambda values
init.lambda <- max.rate / num.states
lambda <- matrix(0, nrow = num.states, ncol = num.states)
lambda[c(upper.tri(lambda) | lower.tri(lambda))] <- init.lambda

#common matrices
diag.a <- diag(total.animals)
diag.o <- diag(1)
mat.a0 <- matrix(0, total.animals, total.animals)
mat.s0 <- matrix(0, num.states, num.states)

########################
####model and params####
########################
#types of state - leader & follower
move.code <- c(rep('O', total.animals), rep('B', num.bm))

#num of params for each state
num.params <- rep(0, num.states)
num.params[move.code == 'B'] <- 3
move.offset <- c(0, cumsum(num.params)[-num.states])
total.params <- sum(num.params)

para.move <- rep(NA, total.params)
para.prop <- rep(NA, total.params)

para.move <- set.params(para.move, alpha, sigma, rho)

#init states
animals[, col.state, ] <- sample(bm.states, size = num.obs * total.animals, replace = T)
#randomly add in subordinates, default to BM if it creates a cyclic structure
for(i in 2:num.obs)
{
    rnd.animals <- sample(1:total.animals)
    
    for(j in rnd.animals)
    {
        animals[i, col.state, j] <- sample(1:num.states, size = 1)
        
        if(is.cyclic(animals[i, col.state, ]))
            animals[i, col.state, j] <- sample(bm.states, size = 1)
    }
}
##############
####output####
##############
file.alpha <- paste("chains/group-alpha-", dat, ".txt", sep = "")
file.sigma <- paste("chains/group-sigma-",  dat, ".txt", sep = "")
file.rho <- paste("chains/group-rho-",  dat, ".txt", sep = "")
file.rates <- paste("chains/group-rates-",  dat, ".txt", sep = "")
file.states1 <- paste("chains/group-1states-", dat, ".txt", sep = "")
file.states2 <- paste("chains/group-2states-", dat, ".txt", sep = "")
file.states3 <- paste("chains/group-3states-", dat, ".txt", sep = "")
file.states4 <- paste("chains/group-4states-", dat, ".txt", sep = "")
file.states5 <- paste("chains/group-5states-", dat, ".txt", sep = "")
file.obs <- paste("chains/group-obs-", dat, ".txt", sep = "")
file.times <- paste("chains/group-times-", dat, ".txt", sep = "")
file.acc.rates <- paste("chains/group-acc-", dat, ".txt", sep = "")

#length of mcmc iteration
iter <- 10000

traj.min <- 3
traj.max <- 12
traj.iter <- ceiling(num.obs / ((traj.min + traj.max) / 2)) * total.animals

#control of output
thin.move <- 2
thin.behav <- 20

#Met-Has acceptance rates
acc.rates <- matrix(0, nrow = 2, ncol = 2)

#################################
####set up switching matrices####
#################################
num.accepted <- num.obs - 1
num.start.switches <- num.accepted * total.animals
a.switches <- array(NA, c(num.start.switches, num.cols, total.animals))

for(i in 1:num.accepted)
{
    #sample switching times
    time.sim <- sort(runif(total.animals))
    interval <- (((i - 1) * total.animals) + 1):(i * total.animals)
    
    a.switches[interval, col.time, ] <- (1 - time.sim)*animals[i, col.time, 1] + time.sim*animals[i + 1, col.time, 1]
    
    cyclic <- TRUE
    
    while((cyclic))
    {
        cyclic <- FALSE
        #order of switches
        rnd.animals <- sample(1:total.animals, size = total.animals)
        
        a.switches[interval, col.type, ] <- rnd.animals
        a.switches[interval, col.jump, ] <- 0
        
        for(j in 1:length(interval))
        {
            a.switches[interval[j], col.id, ] <- 1:total.animals
            
            #state of the switch animal is taken from the next obs
            #class as a jump if it's a switch
            a.switches[interval[j], col.state, rnd.animals[j]] <- animals[i + 1, col.state, rnd.animals[j]]
            if(animals[i + 1, col.state, rnd.animals[j]] != animals[i, col.state, rnd.animals[j]])
                a.switches[interval[j], col.jump, rnd.animals[j]] <- 1
            
            #when j = 1, the states of the non-switch animals are taken from the previous obs
            #otherwise, it's taken from the previous switching time
            if(j == 1)
            {
                a.switches[interval[j], col.state, -rnd.animals[j]] <- animals[i, col.state, -rnd.animals[j]]
            }
            else
            {
                a.switches[interval[j], col.state, -rnd.animals[j]] <- a.switches[interval[j - 1], col.state, -rnd.animals[j]]
            }
            
        }
        
        #if the result is cyclic, start again
        for(j in 1:length(interval))
        {
            if(is.cyclic(a.switches[interval[j], col.state, ]))
                cyclic <- TRUE
        }
    }
}

#a blank switches matrix to avoid having to create a new one for every sim - just needs to large enough row-wise
#used to stop a memory leak
switch.rows <- 1000
switches <- array(0, c(switch.rows, num.cols, total.animals))
switches[, col.x, ] <- NA
switches[, col.y, ] <- NA
switches[, col.time, ] <- NA
switches[, col.id, ] <- rep(1:total.animals, rep(switch.rows, total.animals))
switches[, col.state, ] <- NA
switches[, col.jump, ] <- 0

source('gfit-iter.R')