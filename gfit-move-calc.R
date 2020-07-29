#' Calculate the distribution of a movement step
#'
#' @param params vector of movement params
#' @param states vector of states for all \code{mates}
#' @param deltat time span of the movement step
#' @param current.x vector of current x values for all animals in \code{mates}
#' @param current.y vector of current y values for all animals in \code{mates}
#' @param mates vector of animals that depend on each other in this step
#' @return A list of \code{mu.x} - expected x values for all \code{mates};
#' \code{mu.y} - expected y values for all \code{mates};
#' \code{covar} - covariance matrix of \code{mates};
#' \code{expFt} - exponential matrix of the attraction matrix of \code{mates};
#' \code{texpFt} - transpose of \code{expFt}
move.fn <- function(params, states, deltat, current.x, current.y, mates)
{
    code <- move.code[states]
    off <- move.offset[states]
    
    num.mates <- length(mates)
    leaders <- (1:num.mates)[code == 'B']
    
    #simple cases first
    if(num.mates == 1)
    {
        rho <- params[off + 3]
        
        centre.x <- mu.x <- current.x
        centre.y <- mu.y <- current.y
        
        expFt <- texpFt <- diag.o
        covar <- matrix((rho ^ 2) * deltat, 1, 1)
    }
    else if(num.mates == length(leaders))
    {
        rho <- params[off + 3]
        
        centre.x <- mu.x <- current.x
        centre.y <- mu.y <- current.y
        
        covar <- expFt <- texpFt <- diag.a[1:num.mates, 1:num.mates]
        covar <- covar * (rho ^ 2) * deltat
    }
    else if(num.mates == 2)
    {
        followers <- (1:num.mates)[-leaders]
        
        # only ok if onlt 1 alpha and sigma
        alpha <- params[1]
        sigma <- params[2]
        rho <- params[off[leaders] + 3]
        
        eat <- exp(-alpha * deltat)
        
        #set centre for bm animals
        centre.x <- current.x
        centre.y <- current.y
        centre.x[followers] <- centre.x[leaders]
        centre.y[followers] <- centre.y[leaders]
        
        #each animals dominant
        dom.each <- 1:num.mates
        dom.each[followers] <- leaders
        
        #init cov and expFt matrices
        covar <- expFt <- texpFt <- diag.a[1:2, 1:2]
        covar[leaders, leaders] <- (rho ^ 2) * deltat
        
        #common expressions
        rho.norm <- (rho ^ 2) * (1 / (2 * alpha))
        sigma.norm <- (sigma ^ 2) * (1 / (2 * alpha))
        
        expFt[followers, followers] <- texpFt[followers, followers] <- eat
        expFt[followers, leaders] <- texpFt[leaders, followers] <- 1 - eat
        
        covar[leaders, followers] <- covar[followers, leaders] <- covar[leaders, leaders] -
            (rho ^ 2) * (1 / alpha) * expFt[followers, leaders]
        covar[followers, followers] <- covar[followers, leaders] -
            rho.norm * (expFt[followers, leaders] ^ 2) +
            sigma.norm * (1 - expFt[followers, followers] * expFt[followers, followers])
        
        mu.x <- expFt %*% (current.x - centre.x) + centre.x
        mu.y <- expFt %*% (current.y - centre.y) + centre.y
    }
    else
    {
        followers <- (1:num.mates)[-leaders]
        
        #if only 1 alpha and sigma
        alpha <- params[1]
        sigma <- params[2]
        rho <- params[off[leaders] + 3]
        
        #common expressions
        rho.norm <- (rho ^ 2) * (1 / (2 * alpha))
        sigma.norm <- (sigma ^ 2) * (1 / (2 * alpha))
        eat <- exp(-alpha * deltat)
        
        #set centre for bm animals
        centre.x <- current.x
        centre.y <- current.y
        
        #each animals dominant
        dom.each <- match(states, mates)
        dom.each[leaders] <- leaders
        
        #init cov and expFt matrices
        covar <- expFt <- texpFt <- diag.a[1:num.mates, 1:num.mates]
        
        for(i in 1:length(leaders))
        {
            #the index of the leader
            monarch <- leaders[i]
            upper <- NULL
            covar[monarch, monarch] <- (rho[i] ^ 2) * deltat
            
            #which indexes are sub to the leader
            next.level <- followers[states[followers] == mates[monarch]]
            
            while(length(next.level) > 0)
            {
                tmp.level <- NULL
                
                for(j in next.level)
                {
                    ##########################
                    ## centre of attraction ##
                    ##########################
                    centre.x[j] <- centre.x[monarch]
                    centre.y[j] <- centre.y[monarch]
                    j.dom <- dom.each[j]
                    
                    ################
                    ## expFt calc ##
                    ################
                    expFt[j, j] <- texpFt[j, j] <- eat
                    
                    prev <- j
                    tmp.dom <- j.dom
                    count <- 1
                    
                    #all expFt values for i except with it's leader
                    while(tmp.dom != monarch)
                    {
                        expFt[j, tmp.dom] <- texpFt[tmp.dom, j] <- ((alpha * deltat) / count) * expFt[j, prev]
                        
                        prev <- tmp.dom
                        tmp.dom <- dom.each[tmp.dom]
                        count <- count + 1
                    }
                    
                    #expFt value for i and it's leader
                    expFt[j, monarch] <- texpFt[monarch, j] <- expFt[j.dom, monarch] - expFt[j, prev]
                    
                    #####################
                    ## covariance calc ##
                    #####################
                    rho.norm.vec <- rho.norm[i] * expFt[j , monarch]
                    sigma.norm.vec <- sigma.norm * expFt[j , -monarch]
                    
                    #cov of animal j with it's leader
                    covar[monarch, j] <- covar[j, monarch] <- covar[j.dom, monarch] - (rho.norm.vec * 2)
                    
                    #cov of animal j with all others above it
                    for(k in upper)
                    {
                        covar[j, k] <- covar[k, j] <- (covar[j, dom.each[k]] + covar[j.dom, k]) / 2 -
                            rho.norm.vec * expFt[k, monarch] - sum(sigma.norm.vec * expFt[k, -monarch])
                    }
                    
                    #variance of animal j
                    covar[j, j] <- covar[j, j.dom] - rho.norm.vec * expFt[j, monarch] +
                        sigma.norm - sum(sigma.norm.vec * expFt[j, -monarch])
                    
                    #remove j from list of followers that need dealing with
                    followers <- followers[followers != j]
                    upper <- c(upper, j)
                    tmp.level <- c(tmp.level, followers[states[followers] == mates[j]])
                }

                next.level <- tmp.level
            }
        }
        
        mu.x <- expFt %*% (current.x - centre.x) + centre.x
        mu.y <- expFt %*% (current.y - centre.y) + centre.y
    }
    
    
    list(mu.x = mu.x, mu.y = mu.y, covar = covar, expFt = expFt, texpFt = texpFt)
}