#' Computes the second order basis function
#' for a set of values at a given set of knots.
#' @param x set of values
#' @param x0 set of knots
#' @export
#' @examples
#' ## set of knots
#' x0 <- 0:5
#' ## set of points to evaluate
#' x1 <- seq(0, 5, 0.1)
#' B <- bs2f(x1, x0)
#' plot(x1, B[,1], type='n', xlab='x', ylab='f(x)')
#' for (j in 1:ncol(B)) lines(x1, B[,j])
bs2f <- function(x, x0) {
  i <- findInterval(
    x, x0, all.inside=TRUE)
  k <- length(x0)
  x0 <- c(x0[1]-(x0[2]-x0[1]), x0,
          x0[k]+(x0[k]-x0[k-1]))
  d <- diff(x0)
  d2 <- diff(x0, 2)
  ul <- (x-x0[i+1])/d[i+1]
  ur <- 1-ul
  xl <- ul^2*d[i+1]/d2[i+1]
  xr <- ur^2*d[i+1]/d2[i]
  a <- sparseMatrix(
    i=rep(1:length(i), 3),
    j=c(i+2L, i+1L, i),
    x=c(xl, 1-(xl+xr), xr))
  a[, 2] <- a[,1]+a[,2]
  a[, ncol(a)-1] <-
    a[, ncol(a)-1] + a[, ncol(a)]
  return(a[,2:(ncol(a)-1)])
}
