#' Calculate the MH ratio for the new trajectory
#'
#' @param sim.animal integer of the animal that is being sampled
#' @param sim.data array of the newly simulated trajectory. Contains the 
#' observations and switching times for all animals
#' @param ind.obs vector of indices of observations in \code{sim.data}
#' @param sub.obs array of observations for the trajectory for all animals
#' @param a.switches array of all switches for all animals
#' @param first.obs boolean of whether the trajectory contains the first obs
#' @param last.obs boolean of whether the trajectory contains the last obs
#' @return a list of \code{hr} - hastings ratio; \code{sim.data}- the potentially
#' updated sim.data
move.like.slim <- function(sim.animal, sim.data, ind.obs, sub.obs, a.switches, first.obs, last.obs)
{
    old.like.move <- 0
    new.like.move <- 0
    
    #################################
    ## Piece together the old traj ##
    #################################
    #indices of "accepted switches" happening between time.beg and time.end
    which.accepted <- which(a.switches[, col.time, 1] < time.end & a.switches[, col.time, 1] > time.beg)
    old.sim.data <- array.bind(a.switches[c(0, which.accepted), , , drop = FALSE], sub.obs)
    
    #order data in time
    ord <- order(old.sim.data[, col.time, 1])
    old.sim.data <- old.sim.data[ord, , ]
    
    old.len.sim <- nrow(old.sim.data)
    old.ind.obs <- (1:old.len.sim)[old.sim.data[, col.type, 1] < 1]
    
    #######################################################################
    ## sample missing locations in first obs (if partial) in the new sim ##
    #######################################################################
    if(sim.data[1, col.type, sim.animal] == type.partial & first.obs)
    {
        sim.data[1, col.x, sim.animal] <- rmvnorm_cut(1, partial.var) + sim.data[1, col.x, sim.animal]
        sim.data[1, col.y, sim.animal] <- rmvnorm_cut(1, partial.var) + sim.data[1, col.y, sim.animal]
    }
    
    ################
    ## Likelihood ##
    ################
    sub.length <- length(ind.obs) #same for both old and new
    
    for(i in 1:(sub.length - 1))
    {
        start.obs <- ind.obs[i]
        end.obs <- ind.obs[i + 1]
        
        old.start.obs <- old.ind.obs[i]
        old.end.obs <- old.ind.obs[i + 1]
        
        #if both old and new sim has exact same states in step, no need to calc likelihood
        any.info <- !same.states(sim.data[start.obs:end.obs, col.state, sim.animal],
                                 old.sim.data[old.start.obs:old.end.obs, col.state, sim.animal])
        #but we want partial obs in middle of sim and the last obs (if partial)
        any.partial <- ((last.obs | (i + 1) < sub.length) & sim.data[end.obs, col.type, sim.animal] == type.partial) |
            ((first.obs | i > 1) & sim.data[start.obs, col.type, sim.animal] == type.partial)
        
        if(any.info | any.partial)
        {
            all.points <- (start.obs:end.obs)       
            switch.points <- all.points[c(-1, -length(all.points))]
            #only need to use switches with a jump as there is not info to
            #be gained from non-switches
            info.points <- start.obs
            for(j in switch.points)
            {
                if(any(sim.data[j, col.jump, ] == 1))
                    info.points <- c(info.points, j)
            }
            info.points <- c(info.points, end.obs)
            
            old.all.points <- (old.start.obs:old.end.obs)
            old.switch.points <- old.all.points[c(-1, -length(old.all.points))]
            #only need to use switches with a jump as there is not info to
            #be gained from non-switches
            old.info.points <- old.start.obs
            for(j in old.switch.points)
            {
                if(any(old.sim.data[j, col.jump, ] == 1))
                    old.info.points <- c(old.info.points, j)
            }
            old.info.points <- c(old.info.points, old.end.obs)
            
            #only need to calculate movement distributions for associated animals
            all.states <- rbind(sim.data[info.points, col.state, ],
                                old.sim.data[old.info.points, col.state, ])
            mates <- find.mates(sim.animal, all.states)
            
            ##################
            ## new sim calc ##
            ##################
            #init path with obs data (covar = 0 as it is known)
            exp.x <- sim.data[start.obs, col.x, mates]
            exp.y <- sim.data[start.obs, col.y, mates]
            covar <- matrix(0, nrow = length(mates), ncol = length(mates))
            
            for(j in 2:length(info.points))
            {
                states <- sim.data[info.points[j - 1], col.state, mates]
                deltat <- sim.data[info.points[j], col.time, 1] - sim.data[info.points[j - 1], col.time, 1]
                
                next.step <- step.dist(states, deltat, para.move, exp.x, exp.y, covar, mates)
                
                exp.x <- next.step$exp.x
                exp.y <- next.step$exp.y
                covar <- next.step$covar
            }
            
            ##################
            ## old sim calc ##
            ##################
            #init path with obs data (covar = 0 as it is known)
            old.exp.x <- old.sim.data[old.start.obs, col.x, mates]
            old.exp.y <- old.sim.data[old.start.obs, col.y, mates]
            old.covar <- matrix(0, nrow = length(mates), ncol = length(mates))
            
            for(j in 2:length(old.info.points))
            {
                states <- old.sim.data[old.info.points[j - 1], col.state, mates]
                deltat <- old.sim.data[old.info.points[j], col.time, 1] - old.sim.data[old.info.points[j - 1], col.time, 1]
                
                next.step <- step.dist(states, deltat, para.move, old.exp.x, old.exp.y, old.covar, mates)
                
                old.exp.x <- next.step$exp.x
                old.exp.y <- next.step$exp.y
                old.covar <- next.step$covar
            }
            
            #####################################
            ## likelihood and partial obs sims ##
            #####################################
            partial <- sim.data[end.obs, col.type, sim.animal] == type.partial
            
            #if i is a partial obs for sim.animal; i is in the middle of the
            #trajectory or at the end and last.obs == T; sim.animal is an isolated animal
            if((last.obs | (i + 1) < sub.length) & partial & length(mates) == 1)
            {
                sim.data[end.obs, col.x, sim.animal] <- rmvnorm_cut(1, covar) + exp.x
                sim.data[end.obs, col.y, sim.animal] <- rmvnorm_cut(1, covar) + exp.y
            }
            else
            {
                #if i is a partial obs for sim.animal; i is in the middle of the
                #trajectory or at the end and last.obs == T
                if((last.obs | (i + 1) < sub.length) & partial)
                {
                    miss <- 1 #first animal will always be sim animal in exp x etc
                    obs <- -1
                    #just want the likelihood for the animals associated with sim.animal
                    #not sim.animal itself as it's not observed
                    mates.like <- mates[obs]
                    
                    mu.miss.x <- exp.x[miss] +
                        covar[miss, obs] %*% solve(covar[obs, obs]) %*%
                        (sim.data[end.obs, col.x, mates.like] - exp.x[obs])
                    mu.miss.y <- exp.y[miss] +
                        covar[miss, obs] %*% solve(covar[obs, obs]) %*%
                        (sim.data[end.obs, col.y, mates.like] - exp.y[obs])
                    covar.miss <- covar[miss, miss] - 
                        covar[miss, obs] %*% solve(covar[obs, obs]) %*% covar[obs, miss]
                    
                    sim.data[end.obs, col.x, sim.animal] <- rmvnorm_cut(1, covar.miss) + mu.miss.x
                    sim.data[end.obs, col.y, sim.animal] <- rmvnorm_cut(1, covar.miss) + mu.miss.y
                }
                else
                {
                    #sim.animal is observed and we want likelihood for all animals associated with sim.animal
                    obs <- 1:length(mates)
                    mates.like <- mates
                }
                
                new.like.move <- new.like.move + dmvnorm_cut_log(sim.data[end.obs, col.x, mates.like],
                                                                 mean = exp.x[obs],
                                                                 sigma = covar[obs, obs])
                
                new.like.move <- new.like.move + dmvnorm_cut_log(sim.data[end.obs, col.y, mates.like],
                                                                 mean = exp.y[obs],
                                                                 sigma = covar[obs, obs])
                
                old.like.move <- old.like.move + dmvnorm_cut_log(old.sim.data[old.end.obs, col.x, mates.like],
                                                                 mean = old.exp.x[obs],
                                                                 sigma = old.covar[obs, obs])
                
                old.like.move <- old.like.move + dmvnorm_cut_log(old.sim.data[old.end.obs, col.y, mates.like],
                                                                 mean = old.exp.y[obs],
                                                                 sigma = old.covar[obs, obs])
            }
        }
    }
    
    #hastings ratio
    hr <- exp(new.like.move - old.like.move)
    
    return(list(hr = hr, sim.data = sim.data))
}