#' Functions for the ar2 model
#'
#' Compute the autocorrelation function
#' given the partial correlation coefficients.
#'
#' @param p1 first lag correlation.
#' @param p2 partial correlation of second lag.
#' @param k max lag to compute it.
#' @return vector of length k or matrix of k columns
#' with autocorrelation up to lag k for each
#' parameter pairs, each line if matrix.
#' @describeIn ar2acf marginal correlation.
#' @examples
#' ar2acf(-1, 0.5, 10)
#' ar2acf(c(-1,-1.8), c(0.5, 0.5))
#' @export
#' @md
ar2acf <- function(p1, p2, k) {
  r <- cbind(p1/(1-p2), (p1^2+p2-p2^2)/(1-p2))
  if (k>2) {
    r <- cbind(r, matrix(NA, nrow(r), k-2))
    for (j in 3:k)
      r[,j] <- p1*r[,j-1]+p2*r[,j-2]
  }
  return(drop(r))
}

#' @title ar2q
#' @describeIn precision matrix.
#' @param n is the size of the precision matrix
#' @param a is a length three vector of parameters
#'
#' @examples
#' ar2q(10, c(1, -1, 0.5))
#' @export
#' @md
ar2q <- function(n, a) {
  if (n<1) return(NULL)
  if (n==1) return(a[1]^2)
  if (n==2) return(matrix(
    c(a[1]^2, sum(a[1:2]))[c(1,2,2,1)], 2))
  if (n==3) return(matrix(
    c(a[1]^2, a[1]*a[2], a[1]*a[3],
      a[1]*a[2], sum(a[1:2]^2), a[1]*a[2],
      a[1]*a[3], a[1]*a[2], a[1]^2), 3))
  if (n>3)
    return(Matrix::sparseMatrix(
      i=c(1:n, 1:(n-1), 1:(n-2)),
      j=c(1:n, 2:n, 3:n),
      x=c(a[1]^2, a[1]^2+a[2]^2,
          rep(sum(a^2), max(0,n-4)),
          a[1]^2+a[2]^2, a[1]^2,
          a[1]*a[2],
          rep(a[2]*(a[1]+a[3]), n-3),
          a[1]*a[2],
          rep(a[1]*a[3], n-2)),
      symmetric = TRUE))
}
