#' @title rgb colors from a numeric vector
#' @description create rgb colors from a numeric vector
#' @param x numeric vector
#' @param breaks non decreasing numeri vector
#' @param u logical. Uses the rank of \code{x} as \code{x}.
#' @param reverse logic. If FALSE colors are blue for lower
#' values of \code{x} and red for hight values of \code{x}.
#' If TRUE, colors are red for lower values of \code{x} and
#' blue for blue for hight values of \code{x}.
#' @param ... additional parameters passed to \code{\link{rgb}}.
#' @section Warning:
#' 'transparent' is returned for NA
#' @seealso \code{\link{rgb}}
#' @return 'transparent' if NA or the output of the
#' \code{\link{rgb}} function
#' @export
#' @examples
#' plot(0:5, pch=19, col=x2rgb(0:5))
#' plot(0:5, pch=19, col=x2rgb(0:5, c(0,3,6)))
x2rgb <- function(x, breaks=NULL,
                  u=FALSE, reverse=FALSE, ...) {
### define rgb colors from a numeric vector
    r <- rep('transparent', length(x))
    i.ok <- which(complete.cases(x))
    x <- x[i.ok]
    if (u) {
        x <- rank(x)/length(x)
    } else {
        if (!is.null(breaks)) {
            x <- findInterval(x, breaks)
        }
        rx <- range(x)
        dx <- diff(rx)
        if (dx<sqrt(.Machine$double.eps))
            dx <- sqrt(.Machine$double.eps)
        x <- (x-rx[1])/dx
    }
    if (reverse) {
        r[i.ok] <- rgb(1-x, 1-2*abs(0.5-x), x, ...)
    } else {
        r[i.ok] <- rgb(x, 1-2*abs(x-0.5), 1-x, ...)
    }
    return(r)
}

