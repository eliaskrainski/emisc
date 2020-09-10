#' Frequency table for a numeric vector
#' @description Uses the \code{hist} function
#' to build a frequency table.
#' @param x numeric vector.
#' @param sep character to be used to build the
#' rownames from the breaks.
#' @param digits integer used to format the breaks.
#' @param ... additional arguments passed to \code{hist}
#' @return a table of frequencies.
#' @export
#' @examples
#' htable(rnorm(400))
htable <- function(x, sep=' to ', digits=2, ...) {
    h <- hist(x, ...)
    nb <- length(h$breaks)
    fbk <- format(h$breaks, digits=digits)
    l <- paste(fbk[-nb], fbk[-1], sep=sep)
    return(data.frame(
        mid=h$mids, freq=h$counts,
        rel=100*h$counts/sum(h$counts),
        row.names=l))
}

