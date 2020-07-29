#' Tabulates how many switches go from one state to another
#'
#' @param a1 vector of the states of all animals before their switches
#' @param a2 vector the states of all animals after their switches
#' @param num.states number of states
#' @return a matrix of transitions
state.tab <- function(a1, a2, num.states)
{
    dims <- c(num.states, num.states)
    pd <- prod(dims)
    
    #gives each transition a unique value
    bin <- a1 + num.states*(a2 - 1)
    
    return(array(tabulate(bin, pd), dims))
}
