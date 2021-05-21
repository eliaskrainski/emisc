#' Function for Real to Real transform
#' of Real to Positive functions
#' and its inverse transformation
#'
#' Compute the the desired transformation
#' for the abs of the input and multiply
#' the output by the sign of the input.
#'
#' @param x numeric
#' @param transf character name of the function
#' that implements the transformation to be applied.
#' @param inverse logical indicating if the inverse
#' of the transformation is to be performed
#' @param base numeric to be considered as the
#'  base of the log transformation
#' @return numeric with the transformed input
#' @examples
#' transfR(-5:5, 'sqrt')
#' transfR(-5:5, 'log')
#' plot(function(x) transfR(x, 'log'), -10, 10)
#' plot(function(x) transfR(x, 'log', base=10),
#'      -10, 10, col=4, add=TRUE)
#' x <- -30:30/30
#' plot(x, transfR(transfR(x, 'log'),
#'      'log', inverse=TRUE))
#' @export
#' @md
transfR <- function(x, transf, inverse=FALSE, base=2) {
  s <- sign(x)
  x <- abs(x)
  if (transf=='sqrt') {
    if (inverse) {
      return(s*(x^2))
    } else {
      return(s*sqrt(x))
    }
  }
  if (transf=='log') {
    if (inverse) {
      x <- base^x
      ii <- x<2
      x[ii] <- (x[ii] -1)*2
      return(x*s)
    } else {
      ii <- x<2
      x[ii] <- 1 + x[ii]*0.5
      return(log(x, base)*s)
    }
  }
}
