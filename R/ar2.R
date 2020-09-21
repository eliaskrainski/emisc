#' Functions for the ar2 model
#'
#' Compute the autocorrelation function
#' given the partial correlation coefficients.
#'
#' @param a1 first difference parameter.
#' @param a2 second difference parameter or
#' partial correlation of order two.
#' @param k max lag to compute the acf.
#' @return vector of length k or matrix of k columns
#' with autocorrelation up to lag k for each
#' parameter pairs, each line if matrix.
#' @describeIn ar2acf marginal correlation.
#' @examples
#' ar2acf(1, -0.5, 10)
#' ar2acf(c(1,1.8), c(-0.5, -0.5), 10)
#' @export
#' @md
ar2acf <- function(a1, a2, k) {
  r <- cbind(-a1/(1+a2), (a1^2-a2+a2^2)/(1+a2))
  if (k>2) {
    r <- cbind(r, matrix(NA, nrow(r), k-2))
    for (j in 3:k)
      r[,j] <- -a1*r[,j-1]-a2*r[,j-2]
  }
  return(drop(r))
}

#' @title ar2q
#' @describeIn precision matrix.
#' @param n is the size of the precision matrix
#' @param a length three vector with the parameters
#' a_0, a_1 and a_2 in
#'  a_0 y_t + a_1 y_{t-1} + a_2 y_{t-2} = w_t
#' stochastic difference equation.
#'
#' @examples
#' ar2q(10, c(1, 1, -0.5))
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
