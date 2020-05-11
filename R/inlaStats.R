#' @title inlaStats
#' \code{inlaStats} is to extract dic, waic and cpo
#' summaries from an output of an object returned
#' by the 'inla' function of the 'INLA' package.
#' @param r the result of 'inla'
#' @param idx integer vector to specify which elements
#' of the local dic, waic and cpo is to be considered.
inlaStats <- function(r, idx) {
    if(missing(idx))
        idx <- 1:nrow(r$summary.fitted.values)
    c(dic=sum(r$dic$local.dic[idx]),
      waic=sum(r$waic$local.waic[idx]),
      cpo=-sum(log(r$cpo$cpo[idx])))
}
