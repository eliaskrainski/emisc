#' Functions for Real to Real transform
#' of Real to Positive functions
#'
#' Compute the the desired transformation
#' for the abs of the input and multiply
#' the output by the sign of the input.
#'
#' @param x numeric
#' @param transf character name of the function
#' that implements the transformation to be applied.
#' @return numeric with the transformed input
#' @examples
#' transfR(-5:5, 'sqrt')
#' transfR(-5:5, 'log10')
#' plot(function(x) transfR(x, 'log'), -10, 10)
#' plot(function(x) transfR(x, 'log2'), -10, 10,
#'      col=2, add=TRUE)
#' plot(function(x) transfR(x, 'log10'), -10, 10,
#'      col=4, add=TRUE)
#' @export
#' @md
transfR <- function(x, transf) {
  s <- sign(x)
  x <- abs(x)
  if (substr(transf,1,3)=='log') {
    i1 <- which(x<2)
    x[i1] <- 1.0 + 0.5*x[i1]
  }
  return(s * do.call(transf, list(x)))
}
