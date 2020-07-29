for(j in 1:iter)
{
    for(k in 1:traj.iter)
    {
        ###########################
        ####prep sim trajectory####
        ###########################
        #length of considered interval of observations
        sim.len <- sample(traj.min:traj.max, size = 1)
        #indices of first and last selected obs
        sim.point1 <- sample(1:(num.obs - sim.len + 1), size = 1)
        sim.point2 <- sim.point1 + sim.len - 1
        
        #which animal to sim
        sim.animal <- sample(1:total.animals, size = 1)
        
        first.obs <- FALSE
        last.obs <- FALSE
        
        if(sim.point1 == 1)
            first.obs <- TRUE
        
        if(sim.point2 == num.obs)
            last.obs <- TRUE

        #selected observations
        sub.obs <- animals[sim.point1:sim.point2, , ]
        #times and states of beginning and end of sim
        time.beg <- sub.obs[1, col.time, 1]
        time.end <- sub.obs[sim.len, col.time, 1]

        #sub.switches contains the switches in the interval that don't relate to the sim animal
        which.switches.time <- (1:nrow(a.switches))[a.switches[, col.time, 1] < time.end]
        which.switches.time <- which.switches.time[a.switches[which.switches.time, col.time, 1] > time.beg]
        which.switches <- which.switches.time[a.switches[which.switches.time, col.type, 1] != sim.animal]
        sub.switches <- a.switches[c(0, which.switches), , , drop = FALSE]

        ######################
        ####sim trajectory####
        ######################
        sim <- sim.move(sim.animal, sub.obs, sub.switches, first.obs, last.obs)
        sim.data <- sim$sim.data
        ind.obs <- sim$ind.obs
        ind.switches <- sim$ind.switch
        bk <- sim$bk

        if(!bk)
        {
            #compute the likelihood of the trajectory
            like <- move.like.slim(sim.animal, sim.data, ind.obs, sub.obs, a.switches, first.obs, last.obs)
            sim.data <- like$sim.data

            if(runif(1) < like$hr)
            {
                #cat(file = file.updates, sim.point1, " ", sim.point2, " ", sim.len, " ", length(ind.switches), " ", j, "\n", append = TRUE)
                acc.rates[2, 1] <- acc.rates[2, 1] + 1

                #update data states and partial obs
                animals[sim.point1:sim.point2, , ] <- sim.data[ind.obs, , ]
                
                #drop old switches
                if(length(which.switches.time) > 0)
                    a.switches <- a.switches[-which.switches.time, , , drop = FALSE]
                #combine new switches with existing outside the interval
                new.switches <- sim.data[c(0, ind.switches), , , drop = FALSE]
                if(nrow(new.switches) > 0)
                {
                    a.switches <- array.bind(a.switches, new.switches)
                }
            }
        }
    }
    
    num.accepted <- nrow(a.switches)
    
    ################################
    ## Update movement parameters ##
    ################################
    #accepted switches and observations
    all.data <- array.bind(a.switches, animals)
    #order data in time
    ord <- order(all.data[, col.time, 1])
    all.data <- all.data[ord, , ]
    
    num.all.obs <- num.obs + num.accepted
    ind.obs.all <- (1:num.all.obs)[all.data[, col.type, 1] < 1]
    ind.switch.all <- (1:num.all.obs)[-ind.obs.all]
    
    ################################
    ## Update movement parameters ##
    ################################
    #update lambda
    if(num.accepted > 0)
        lambda <- update.rate(all.data, num.all.obs)
    
    #movement param proposals
    new.alpha <- alpha + rnorm(1, 0, prop.alpha)
    new.sigma <- sigma + rnorm(1, 0, prop.sigma)
    
    new.rho <- c(NA, NA)
    new.rho[1] <- rho[1] + rnorm(1, 0, prop.slow)
    new.rho[2] <- rho[2] + rnorm(1, 0, prop.fast)
    
    if(new.alpha > 0 & new.sigma > 0 & all(new.rho > 0) & new.sigma < new.rho[2] & new.rho[2] > new.rho[1])
    {
        para.prop <- set.params(para.prop, new.alpha, new.sigma, new.rho)

        ret <- update.params(para.prop, para.move, all.data, ind.obs.all)

        if(ret$params.update)
        {
            alpha <- new.alpha
            sigma <- new.sigma
            rho <- new.rho

            para.move <- ret$params
            acc.rates[1, 1] <- acc.rates[1, 1] + 1
        }
    }
    
    if(j%%thin.move == 0)
    {
        cat(file = file.alpha, alpha, "\n", append = TRUE)
        cat(file = file.sigma, sigma, "\n", append = TRUE)
        cat(file = file.rho, rho, "\n", append = TRUE)
        # cat(file = file.rates, c(t(lambda)), "\n", append = TRUE)
    }
    
    if(j%%thin.behav == 0)
    {
        #move params percent
        acc.rates[1, 2] <- (acc.rates[1, 1] / thin.behav) * 100
        #traj updates percent
        acc.rates[2, 2] <- (acc.rates[2, 1] / thin.behav) * 100 / traj.iter
        
        obs <- all.data[, col.type, 1]
        obs[obs == type.obs | obs == type.partial] <- "o"
        obs <- which(obs == "o")
        
        cat(file = file.states1, all.data[, col.state, 1], "\n", append = TRUE)
        cat(file = file.states2, all.data[, col.state, 2], "\n", append = TRUE)
        cat(file = file.states3, all.data[, col.state, 3], "\n", append = TRUE)
        cat(file = file.states4, all.data[, col.state, 4], "\n", append = TRUE)
        cat(file = file.states5, all.data[, col.state, 5], "\n", append = TRUE)
        cat(file = file.obs, obs, "\n", append = TRUE)
        cat(file = file.times, all.data[, col.time, 1], "\n", append = TRUE)
        cat(file = file.acc.rates, round(acc.rates[, 2], 2), "\n", append = TRUE)
        cat(file = file.rates, c(t(lambda)), "\n", append = TRUE)

        acc.rates <- matrix(0, nrow = 2, ncol = 2)
    }
}
