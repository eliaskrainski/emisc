#' A more complete summary
#' @description a more complete summary that
#' applyes for both numeric and factor vectors.
#' @param x numeric vector or factor
#' @param p percentiles for computing quantiles
#' when \code{x} is a numeric vector.
#' @return a named vector with the length of \code{x},
#' number of NA`s and mean, sd an quantiles if \code{x}
#' is numeric or frequencies if \code{x} is factor.
#' @export
#' @examples
#' mysummary(iris[,1])
#' mysummary(iris[,2])
mysummary <- function(x, p=c(0, 0.025, 0.25,
                             0.5, 0.75, 0.975, 1)) {
    n <- length(x)
    x <- x[complete.cases(x)]
    r <- c(n=n, na=n-length(x))
    if (is.numeric(x)) {
        q <- quantile(x, p)
        return(c(r, m=mean(x), sd=sd(x), q))
    } else {
        return(c(r, table(x)))
    }
}
