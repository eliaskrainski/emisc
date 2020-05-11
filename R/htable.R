#' \code{htable} provides a frequency table using
#' the \code{hist} function internally.
#' @param x numeric vector
#' @param ... additional arguments passed to \code{hist}
#' @return a table of frequencies.
#' @examples
#' htable(rnorm(100))
htable <- function(x, ...) {
    h <- hist(x, ...)
    nb <- length(h$breaks)
    l <- paste(h$breaks[-nb], h$breaks[-1], sep='-')
    return(data.frame(
        mid=h$mids, freq=h$counts,
        rel=100*h$counts/sum(h$counts),
        row.names=l))
}

