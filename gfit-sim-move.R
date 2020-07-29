#' Simulate a new behaviour trajectory for \code{sim.animal}
#'
#' @param sim.animal integer of the animal that is being sampled
#' @param sub.obs array of observations for the trajectory for all animals
#' @param sub.switches array of switches for all animals except sim.animal
#' @param first.obs boolean of whether the trajectory contains the first obs
#' @param last.obs boolean of whether the trajectory contains the last obs
#' @return a list of \code{sim.data} - array of the new trajectory for all animals;
#' \code{ind.obs}- vector of indices of observations in \code{sim.data};
#' \code{ind.switch}- vector of indices of switching times in \code{sim.data};
#' \code{bk} - boolean on whether the new \code{sim.data} meets all conditions
sim.move <- function(sim.animal, sub.obs, sub.switches, first.obs, last.obs)
{
    #################################
    ## Simulate potential switches ##
    #################################
    bk <- FALSE
    
    #number of potential switches
    num.switches <- rpois(1, (time.end - time.beg) * max.rate)
    
    #if first obs, sample initial state for sim.animal from initial state vector
    if(first.obs)
    {
        sub.obs[1, col.state, sim.animal] <- sample(1:num.states, size = 1)
        
        if(is.cyclic(sub.obs[1, col.state, ]))
            bk <- TRUE
    }
    
    if(num.switches == 0)
    {
        new.switches <- sub.switches
        
        #if there are switches but none are sampled - reject
        if(sub.obs[nrow(sub.obs), col.state, sim.animal] != sub.obs[1, col.state, sim.animal])
        {
            bk <- TRUE
        }
    }
    else
    {
        #initialize simulated switches
        tmp.switches <- switches[1:num.switches, , , drop = FALSE]
        tmp.switches[, col.time, ] <- runif(num.switches, time.beg, time.end)
        tmp.switches[, col.type, ] <- sim.animal
        
        #combine tmp.switches with sub.switches
        new.switches <- array.bind(sub.switches[, , , drop = F], tmp.switches[, , , drop = F])
    }
    
    #all simulated switches are in the first part, data are in second part
    sim.data <- array.bind(new.switches[, , , drop = F], sub.obs)
    len.sim <- nrow(sim.data)
    
    #order data in time
    ord <- order(sim.data[, col.time, 1])
    sim.data <- sim.data[ord, ,]
    
    ind.obs <- (1:len.sim)[sim.data[, col.type, 1] < 1]
    ind.switch <- (1:len.sim)[-ind.obs]
    
    t <- 2
    
    while (t <= len.sim & !bk) {
        #if type is not an obs or partial
        if(sim.data[t, col.type, 1] > 0)
        {
            if(sim.data[t, col.type, 1] == sim.animal)
            {
                #fill in states of non sim-animals
                sim.data[t, col.state, -sim.animal] <- sim.data[t - 1, col.state, -sim.animal]

                #sample new state of sim.animal
                prev.state <- sim.data[t - 1, col.state, sim.animal]
                jump.now <- runif(1) < (sum(lambda[prev.state, ]) / max.rate)
                
                if(jump.now)
                {
                    new.state <- sample(1:num.states, size = 1, prob = lambda[prev.state, ])
                    sim.data[t, col.state, sim.animal] <- new.state
                    sim.data[t, col.jump, sim.animal] <- 1
                    
                    if(is.cyclic(sim.data[t, col.state, ]))
                        bk <- TRUE
                }
                else
                {
                    sim.data[t, col.state, sim.animal] <- sim.data[t - 1, col.state, sim.animal]
                }
            }
            else
            {
                sim.data[t, col.state, sim.animal] <- sim.data[t - 1, col.state, sim.animal]
                
                #if new states are cyclic, stop
                if(any(sim.data[t, col.jump, -sim.animal] == 1))
                {
                    if(is.cyclic(sim.data[t, col.state, ]))
                        bk <- TRUE
                }
            }
        }
        else
        {
            if(t == len.sim & !last.obs)
            {
                if(sim.data[len.sim - 1, col.state, sim.animal] != sim.data[len.sim, col.state, sim.animal])
                    bk <- TRUE
            }
            else
            {
                sim.data[t, col.state, ] <- sim.data[t - 1, col.state, ]
            }
        }
        
        t <- t + 1
    }

    return(list(sim.data = sim.data, ind.obs = ind.obs, ind.switch = ind.switch, bk = bk))
}