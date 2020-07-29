#' Combine two arrays through the third dimension
#'
#' @param first a 3D array
#' @param second a 3D array
#' @return A combined array of \code{first} and \code{second}
array.bind <- function(first, second)
{
    len1 <- nrow(first)
    len2 <- nrow(second)
    len <- len1 + len2
    
    if(len1 == 0)
    {
        both <- second
    }
    else if(len2 == 0)
    {
        both <- first
    }
    else
    {
        both <- array(0, c(len, num.cols, total.animals))
        
        both[1:len1, , ] <- first
        both[(len1 + 1):len, , ] <- second
    }
    
    return(both)
}

#' Calculate the distribution of a movement step
#' Regardless of a stand alone calculation or part of a recursion
#'
#' @param states vector of the states of the animals over the movement step
#' @param deltat the time span of the movement step
#' @param prev.x vector of the prev locations in the x axis, observed or expected, for all animals
#' @param prev.y vector of the prev locations in the y axis, observed or expected, for all animals
#' @param prev.covar matrix of the prev covar, 0 if the previous data was observed
#' @param mates vector of which animals the distribution needs to be calculated for
#' @return list of \code{exp.x} - vector of the expected x values;
#' \code{exp.y} - vector of the expected y values;
#' \code{covar} - covariance matrix
step.dist <- function(states, deltat, params, prev.x, prev.y, prev.covar, mates = 1:total.animals)
{
    move.info <- move.fn(params, states, deltat, prev.x, prev.y, mates)

    #new expected value and covariance
    exp.x <- move.info$mu.x
    exp.y <- move.info$mu.y
    covar <- move.info$covar + move.info$expFt %*% prev.covar %*% move.info$texpFt
    
    return(list(exp.x = exp.x, exp.y = exp.y, covar = covar))
}

#' Return the modes of a vector
#'
#' @param x a vector
#' @return the modes of \code{x}
my.mode <- function(x)
{
    values <- unique(x)
    tab <- tabulate(match(x, values))
    values[tab == max(tab)]
}

#' Set the movement params vector
#'
#' @param params the movement params vector
#' @param alpha
#' @param sigma
#' @param rho
#' @return the movement params vector
set.params <- function(params, alpha, sigma, rho)
{
    params[move.offset[move.code == 'B'] + 1] <- alpha
    params[move.offset[move.code == 'B'] + 2] <- sigma
    params[move.offset[move.code == 'B'] + 3] <- rho
    
    params
}

#' Returns a boolean on whether two behaviour trajectories contain the same
#' singular state
#'
#' @param x a vector of states
#' @param y a vector of states
#' @return a boolean on whether there is only a singular value between both params
same.states <- function(x, y)
{
    total <- c(x, y)
    unique.total <- unique(total)

    #true if only 1 state across old and new sim
    return(length(unique.total) == 1)
}

#' Returns a boolean on whether the states of a group produces a cyclic structure
#'
#' @param states vector of the states of a group
#' @return a boolean on whether \code{states} produces a cyclic structure
is.cyclic <- function(states)
{
    cyclic <- FALSE
    
    code <- move.code[states]
    leaders <- (1:total.animals)[code == 'B']
    
    if(length(leaders) == 0)
    {
        cyclic <- TRUE
    }
    else
    {
        states[leaders] <- 0 #0 indicates BM here
        dom <- states
        followers <- (1:total.animals)[dom > 0]
        
        while(length(followers) > 0 & cyclic == FALSE)
        {
            if(any(dom == 1:total.animals))
            {
                cyclic <- TRUE
            }
            else
            {
                dom[followers] <- states[dom[followers]]
                followers <- (1:total.animals)[dom > 0]
            }
        }
    }
    
    return(cyclic)
}

#' Takes the states of the group over a period of time and returns a vector of
#' all the animals that sim.animal is associated with. That is, animals that were
#' at some point in the same subgroup as sim.animal, and their own associates.
#'
#' @param sim.animal integer of the animal in focus
#' @param states a matrix of the states of a group over a period of time
#' @return a vector of animals sim.animal is associated with in \code{states}
find.mates <- function(sim.animal, states)
{
    if(all(states[, sim.animal] >= bm.states[1]) &
       all(states[, -sim.animal] != sim.animal))
    {
        return(sim.animal)
    }
    else
    {
        mates <- sim.animal
        
        new.sub <- T
        new.dom <- T
        while(new.sub | new.dom)
        {
            new.sub <- F
            new.dom <- F
            
            #foes are animals not (yet) linked to sim.animal
            foes <- (1:total.animals)[-mates]    
            for(i in foes)
            {
                #if i is subordinate to anyone in mates
                if(any(states[, i] %in% mates))
                {
                    mates <- c(mates, i)
                    new.sub <- T
                }
            }
            
            #foes are animals not (yet) linked to sim.animal
            foes <- (1:total.animals)[-mates]
            for(i in foes)
            {
                #if i is the dominant to anyone in mates
                if(any(states[, mates] == i))
                {
                    mates <- c(mates, i)
                    new.dom <- T
                }
            }
        }
        
        return(mates)
    }
}

#' Determines in two animals are in the same subgroup
#'
#' @param a1 integer of an animal
#' @param a2 integer of an animal
#' @param states the states of the group
#' @return a boolean on whether \code{a1} and \code{a2} are in the same subgroup
same.group <- function(a1, a2, states)
{
    same <- FALSE
    
    leaders <- (1:total.animals)[states >= bm.states[1]]
    
    if(states[a1] >= bm.states[1] & states[a2] >= bm.states[1])
    {
        same <- FALSE
    }
    else
    {
        leader.each <- 1:total.animals
        followers <- (1:total.animals)[-leaders]
        
        remain <- followers
        
        while(length(remain) > 0)
        {
            leader.each[remain] <- states[leader.each[remain]]
            
            remain <- (1:total.animals)[leader.each %in% followers]
        }
        
        if(leader.each[a1] == leader.each[a2])
            same <- TRUE
    }
    
    return(same)
}