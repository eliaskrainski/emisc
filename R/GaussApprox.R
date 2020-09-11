#' Computes the Gaussian approximation
#' given the likelihood function, design matrix and
#' the prior precision matrix.
#' @param f a function that computes the likelihood
#' @param n the number of data for the likelihood
#' @param Q the joint precision matrix for the problem
#' @param k the number of iterations
#' @param h the change around the avaluation point
#' @param x initial value vector
#' @export
#' @examples
#' ## Consider confirmed cases of COVID19 in NS
#' n <- 70
#' x <- sin(2*pi*1:n/n)
#' y <- rpois(n, exp(1+x))
#' par(mfrow=c(1,1), mar=c(2,3,0.5,0.5), mgp=c(2,0.5,0))
#' plot(y, pch=19)
#' ## define the likeihood function
#'  llp <- function(e)
#'   dpois(y, exp(e), log=TRUE)
#' ## define the precision matrix structure
#' Q2 <- crossprod(diff(Diagonal(n), differences=2))
#' ga <- GaussApprox(llp, n, Q2*5000)
#' ga$V <- solve(ga$Q)
#' ga$sd <- sqrt(diag(ga$V))
#'
#' ## visualize the result
#' plot(y, pch=19)
#' lines(exp(1+x), col=4, lwd=3)
#' lines(exp(ga$x), col=2)
#' lines(exp(ga$x - 1.96*ga$sd), col=2, lty=2)
#' lines(exp(ga$x + 1.96*ga$sd), col=2, lty=2)
#'
#' ## do it with basis functions
#' x0 <- seq(0, n, 10)
#' B <- bs2f(1:n, x0)
#'
#' ## precision (as in the INLA review paper)
#' Qb <- crossprod(diff(Diagonal(length(x0)-1)))
#' tau.e <- exp(12)
#' Qx <- rbind(
#'   cbind(Diagonal(n, tau.e), tau.e*B),
#'   cbind(tau.e*t(B), 2*Qb + tau.e*crossprod(B)))
#'
#' ga.b <- GaussApprox(llp, n, Qx)
#' ga.b$sd <- sqrt(diag(solve(ga.b$Q)))
#'
#' lines(exp(ga.b$x[1:n]), col=4)
#' lines(exp(ga.b$x[1:n] - 1.96*ga.b$sd[1:n]), col=4, lty=2)
#' lines(exp(ga.b$x[1:n] + 1.96*ga.b$sd[1:n]), col=4, lty=2)
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
GaussApprox <- function(f, n, Q, k=5,
                        h=.Machine$double.eps^0.2,
                        x=NULL, ...) {
    dsmall <- .Machine$double.eps^0.5
    m <- ncol(Q)
    if (is.null(x))
        x <- double(m)
    for (j in 1:k) {
        fx <- f(x[1:n])
        fa  <- f(x[1:n] - h)
        fb  <- f(x[1:n] + h)
        cc <- (2*fx -fa -fb)/(h^2)
        cc[cc<dsmall] <- dsmall
        bb <- (fb - fa)/(2*h) + x[1:n]*cc
        Qn <- Q + Diagonal(m, c(cc, rep(0, m-n)))
        L <- chol(Qn)
        x <- drop(solve(Qn, c(bb, rep(0, m-n))))
    }
    return(list(x=x, bb=bb, cc=cc, Q=Qn))
}
