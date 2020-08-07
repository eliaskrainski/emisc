#' @title GaussApprox
#' \code{GaussApprox} computes the Gaussian approximation
#' given the likelihood function, design matrix and
#' the prior precision matrix.
#' @param f the function that computes the likelihood
#' @param A the design matrix
#' @param Q the prior precision matrix
#' @param k the number of iterations
#' @param h the change around the avaluation point
#' @param ... additional arguments to f
#' @examples
#' data(NSCases)
#' llf <- function(x)
#'     dgamma(cases, exp(x), log=TRUE)
#' n <- length(NSCases$cases)
#' library(Matrix)
#' Rt <- crossprod(diff(Diagonal(n)))
#' tau <- 1500
#' ga <- GaussApprox(llf, Diagonal(n), Rt*tau)
#' with(NSCases, plot(day, cases, pch=8))
#' lines(NSCases$day, ga$mu)
GaussApprox <- function(f, A, Q, k=5,
                        h=.Machine$double.eps^0.2,
                        ...) {
    b <- double(ncol(A))
    for (j in 1:k) {
        fx <- f(b)
        fa  <- f(b - h)
        fb  <- f(b + h)
        cc <- (2*fx -fa -fb)/(h^2)
        cc[cc<0] <- 0
        xx <- (fb - fa)/(2*h) + b*cc
        Qn <- Q + Diagonal(x=cc)
        L <- chol(Qn)
        new <- drop(solve(L, solve(t(L), xx)))
        b <- new
    }
    return(list(mu=b, bb=xx, cc=cc, Q=Qn, L=L,
                sldL=sum(log(diag(L)))))
}
