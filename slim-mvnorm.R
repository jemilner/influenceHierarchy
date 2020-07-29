# radically cut-down version of rmvnorn & dmvnorm from mvtnorm package
# Only uses Cholesky method and hence only upper triangular part 
# of cov matrix - much quicker to calculate than whole matrix 

rmvnorm_cut <- function (ndim,sigma)
    # ndim: number of values to generate - typically 1 per dimension 
    # sigma: version of cov matrix (e.g. for Brownian bridge) but only 
    # diagonal and upper triangle need to be correct
{
    R <- chol(sigma, pivot = TRUE)
    R <- R[, order(attr(R, "pivot"))]
    return(matrix(rnorm(ndim * ncol(sigma)), nrow = ndim, byrow = TRUE) %*% R)
}

# dmvnorm_cut_log <- function (x, mean, sigma) 
#     # x: a vector
# {
#     #p <- length(x)
#     x <- matrix(x, ncol = length(x))
#     p <- ncol(x)
#     dec <- chol(sigma)#, pivot = TRUE)
#     #cat("sigma\n")
#     #print(sigma)
#     #print(chol(sigma))
#     #print(cov2cor(sigma))
#     #cat("dec\n")
#     #print(dec)
#     tmp <- backsolve(dec, t(x) - mean, k=p, transpose = TRUE)
#     rss <- sum(tmp^2)
#     logretval <- -sum(log(diag(dec))) - 0.5 * p * log(2 * pi) - 0.5 * rss
#     return(logretval)
# }

dmvnorm_cut_log <- function (x, mean, sigma) 
{
    p <- length(x)
    dec <- chol(sigma)#, pivot = TRUE)
    tmp <- backsolve(dec, x - mean, k = p, transpose = TRUE)
    rss <- sum(tmp^2)
    logretval <- -sum(log(diag(dec))) - 0.5 * p * log(2 * pi) - 0.5 * rss
    return(logretval)
}