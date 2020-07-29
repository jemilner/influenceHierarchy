#' Calculate the MH ratio for the new movement parameter proposals
#'
#' @param new.params vector of the new sample of movement params
#' @param old.params vector of the current sampple of movement params
#' @param all.data array of augmented data for all animals
#' @param ind.obs vector of indices of the observations in \code{all.data}
#' @return a list of \code{params} - the accepted parameters;
#' \code{params.update} - whether the params were udpated or not
update.params <- function(new.params, old.params, all.data, ind.obs)
{
    old.like <- 0
    new.like <- 0
    
    for(i in 1:(length(ind.obs) - 1))
    {
        start.obs <- ind.obs[i]
        end.obs <- ind.obs[i + 1]
        
        all.points <- (start.obs:end.obs)
        switch.points <- all.points[c(-1, -length(all.points))]

        #only need to use switches with a jump as there is not info to
        #be gained from non-switches
        info.points <- start.obs
        for(j in switch.points)
        {
            if(any(all.data[j, col.jump, ] == 1))
                info.points <- c(info.points, j)
        }
        info.points <- c(info.points, end.obs)
        
        #init path with obs data (covar = 0 as it is known)
        old.exp.x <- new.exp.x <- all.data[start.obs, col.x, ]
        old.exp.y <- new.exp.y <- all.data[start.obs, col.y, ]
        old.covar <- new.covar <- mat.a0
        
        for(j in 2:length(info.points))
        {
            states <- all.data[info.points[j - 1], col.state, ]
            deltat <- all.data[info.points[j], col.time, 1] - all.data[info.points[j - 1], col.time, 1]
            
            #calc movement step distribution with old parameter sample
            old.next.step <- step.dist(states, deltat, old.params, old.exp.x, old.exp.y, old.covar)
            
            old.exp.x <- old.next.step$exp.x
            old.exp.y <- old.next.step$exp.y
            old.covar <- old.next.step$covar
            
            #calc movement step distribution with new parameter sample
            new.next.step <- step.dist(states, deltat, new.params, new.exp.x, new.exp.y, new.covar)
            
            new.exp.x <- new.next.step$exp.x
            new.exp.y <- new.next.step$exp.y
            new.covar <- new.next.step$covar
        }
        
        old.like <- old.like + dmvnorm_cut_log(all.data[end.obs, col.x, ],
                                       mean = old.exp.x,
                                       sigma = old.covar)
        
        old.like <- old.like + dmvnorm_cut_log(all.data[end.obs, col.y, ],
                                       mean = old.exp.y,
                                       sigma = old.covar)
        
        new.like <- new.like + dmvnorm_cut_log(all.data[end.obs, col.x, ],
                                       mean = new.exp.x,
                                       sigma = new.covar)
        
        new.like <- new.like + dmvnorm_cut_log(all.data[end.obs, col.y, ],
                                       mean = new.exp.y,
                                       sigma = new.covar)
    }
    
    hastings <- exp(new.like - old.like)
    
    params.update <- FALSE
    if(runif(1) < hastings)
    {
        params <- new.params
        params.update <- TRUE
    } else {
        params <- old.params
    }
    
    return(list(params = params, params.update = params.update))
}