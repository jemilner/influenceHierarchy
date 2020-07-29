#' Samples heterogeneous transition rates from analytical posterior
#'
#' @param all.data an array of the augmented data for all animals
#' @param num.data number of augmented data
#' @return a matrix of transition rates
update.rate <- function (all.data, num.data)
{
    #locate all switching times for each animal(whether actual switches or not)
    poi <- lapply(1:total.animals, function(x) {(1:num.data)[all.data[, col.type, 1] == x]})
    #find all the states before the switches and after
    before <- unlist(lapply(1:total.animals, function(x) {all.data[poi[[x]] - 1, col.state, x]}))
    after <- unlist(lapply(1:total.animals, function(x) {all.data[poi[[x]], col.state, x]}))
    
    #tabultate how many switches there are from one state to another
    #including switches that don't switch state
    #counts[i, j] = number of switches from state i to state j
    counts <- state.tab(before, after, num.states)
    
    alpha0 <- rep(shape1, num.states)
    lambda <- mat.s0
    
    for (i in 1:num.states)
    {
        #combine tabulated switches with priors
        alpha <- alpha0
        alpha[i] <- shape2
        alpha <- alpha + counts[i, ]
        
        #sample from the posterior, scale by max.rate
        lambda[i, -i] <- max.rate * rdirichlet(n = 1, alpha)[-i]
    }
    
    return(lambda)
}
