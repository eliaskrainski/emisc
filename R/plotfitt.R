#' Visualize fit output against observed values
#' @param o vector of observed values
#' @param sfitt an output of \code{lm}, \code{glm} or
#' a \code{data.frame} containing columns with
#' 'mean', 'sd', lower limit IC, median and upper limit IC.
#' @param asp the aspect ratio
#' @param xlab the x axis title
#' @param ylab the y axis title
#' @param xlim the limits for the x axis scale
#' @param ylim the limits for the y axis scale
#' @param ... additional arguments passed to
#' \code{plot} and \code{arrows}
#' @export
#' @examples
#' d <- list(x=seq(pi, 3*pi, 0.1))
#' d$y <- rpois(length(d$x), exp(1 + d$x))
#' r <- glm(y~x, poisson, d)
#' plotfitt(d$y, r)
plotfitt <- function(o, sfitt, asp=1,
                     xlab='Fitted', ylab='Observed',
                     xlim, ylim, ...) {
  if (any(class(sfitt)%in%c('lm', 'glm'))) {
    o <- sfitt$model[,1]
    p <- predict(
      sfitt, se.fit=TRUE, type='response')
    sfitt <- data.frame(mean=p$fit, se=p$se.fit)
    sfitt$li <- sfitt$mean - 1.96*sfitt$se
    sfitt$median <- sfitt$mean
    sfitt$ls <- sfitt$mean + 1.96*sfitt$se
  }
  oo <- order(o)
  if (missing(xlim))
    xlim <- range(o, sfitt[, c(3,5)])
  if (missing(ylim))
    ylim <- range(o, sfitt[, c(3,5)])
  plot(sfitt[, 1], o, asp=asp,
       xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, ...)
  arrows(sfitt[oo, 3], o[oo],
         sfitt[oo, 5], o[oo],
         lty=3, length=0.03, angle=90, code=3, ...)
  abline(0:1)
  r <- cor(sfitt[, 1], o)
  legend('topleft', '',
         title=paste('R2:', format(r^2, digits=4)), bty='n')
}
