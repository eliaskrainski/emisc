#' Computes the Gaussian approximation
#' given the likelihood function, design matrix and
#' the prior precision matrix.
#' @param f a function that computes the likelihood
#' @param A the design matrix
#' @param Q the prior precision matrix
#' @param k the number of iterations
#' @param h the change around the avaluation point
#' @param x vector with the initial value
#' @param ... additional arguments to f
#' @export
#' @examples
#' ## Consider confirmed cases of COVID19 in NS
#' data(NScases)
#' with(NScases, plot(day, cases, pch=8))
#' ### likelihood
#' llf <- function(b)
#'     dpois(NScases$cases, exp(B%*%b), log=TRUE)
#' n <- length(NScases$cases)
#' c(sum(llf(rep(0,n))), sum(llf(rep(1,n))))
#' ## precision structure matrix of a RW2 model
#' R <- crossprod(diff(Diagonal(n), differences=2))
#' tau <- 5 ## a reasonable value for the precision
#' ### perform GA
#' ga <- GaussApprox(llf, Diagonal(n), R*tau)
#' lines(NScases$day, ga$x)
#'
#'\dontrun{
#' ## Consider Tokyo data from INLA package
#' data('Tokyo', package='INLA')
#' if (any(ls()=='Tokyo')) {
#'   llb <- function(x)
#'      dbinom(Tokyo$y, Tokyo$n, plogis(x), log=TRUE)
#'   n <- nrow(Tokyo)
#' ### cyclic RW1 structure matrix
#'   R1c <- sparseMatrix(
#'     i=c(1:n, 2:n,n,1,1:(n-1)),
#'     j=c(1:n, 1:(n-1),1,n,2:n),
#'     x=rep(c(2,-1), c(n, 2*n)))
#'   ga1 <- GaussApprox(llb, Diagonal(n), R1c*40)
#' ### cyclic RW2 structure matrix
#'   R2c <- crossprod(R1c)
#'   ga2 <- GaussApprox(llb, Diagonal(n), R2c*15000)
#'### visualize
#'   with(Tokyo, plot(time, y/n, pch=19))
#'   lines(Tokyo$time, plogis(ga1$mu))
#'   lines(Tokyo$time, plogis(ga2$mu), col=2)
#' }
#' }
GaussApprox <- function(f, A=NULL, Q, k=5,
                        h=.Machine$double.eps^0.2,
                        x=NULL, ...) {
    dsmall <- .Machine$double.eps^0.5
    if (is.null(x))
        x <- double(ncol(A))
    for (j in 1:k) {
        if (!is.null(A))
            x <- drop(A %*% x)
        fx <- f(x)
        fa  <- f(x - h)
        fb  <- f(x + h)
        cc <- (2*fx -fa -fb)/(h^2)
        cc[cc<dsmall] <- dsmall
        bb <- (fb - fa)/(2*h) + x*cc
        Qe <- Diagonal(length(cc), cc)
        AcA <- crossprod(A, Qe)%*%A
        Qn <- Q + AcA
        L <- chol(Qn)
        x <- drop(solve(Qn, crossprod(A, bb)))
    }
    return(list(x=x, bb=bb, cc=cc, Q=Qn))
}
