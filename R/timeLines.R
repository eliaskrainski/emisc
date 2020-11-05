#' Function to plot multiple time series.
#'
#' @param w numeric matrix, each line a time series data
#' @param x vector with the x-axis values
#' @param xlim numeric length two vector with x-axis limits
#' @param ylim numeric length two vector with y-axis limits
#' @param xlab label for the x-axis
#' @param ylab label for the y-axis
#' @param axis1args arguments list for x-axis
#' @param axis2args arguments list for y-axis
#' @param col colors for each line drawn
#' @param x.expand numeric to be used as fraction of
#' diff(xlim) as extension of it
#' @param y.labs labels to be added like a legend
#' @param y.ylabs numeric vector with y-values to
#' place y.labs
#' @param yh numeric used as h in cutree(hclust(...))
#' to arrange the y.labs towards a more
#' equally spaced sequence over y-limits.
#' NOTE: This is only used if y.ylabs is not NULL.
#' @param ... further arguments passed further
#' @return plot
#' @examples
#' timeLines(t(Seatbelts[,3:4]))
#' @export
#' @md
timeLines <- function(
  w,
  x=1:ncol(w),
  xlim=range(x),
  ylim=range(w, na.rm=TRUE),
  xlab='', ylab='',
  axis1args=list(at=pretty(x)),
  axis2args=list(at=pretty(ylim), las=1),
  col=NULL,
  x.expand=0.35,
  y.labs=rownames(w),
  y.ylabs=NULL,
  yh=0.7*diff(range(
    w[, ncol(w)], na.rm=TRUE))/nrow(w),
  ...)  {
  ny <- nrow(w)
  nt <- ncol(w)
  yy <- w[, nt]
  oy <- order(yy)
  if (is.null(col))
    col <- rgb(
      oy/ny,
      1-2*abs(1:ny/(ny+1)-0.5),
      rev(oy)/ny)
  if(!is.null(y.labs)) {
    x.ylabs <- rep(
      max(xlim) +
        diff(xlim)*x.expand/20, ny)
    xlim[2] <- xlim[2] +
      diff(xlim)*x.expand
    if (is.null(y.ylabs)) {
      dw <- dist(yy)
      hc <- hclust(dw, method='single')
      gg <- cutree(hc, h=yh)
      g.y <- tapply(yy, gg, mean)
      gg <- factor(
        gg, names(g.y)[order(g.y)])
      g.y <- tapply(yy, gg, mean)
      g.n <- table(gg)
      n.g <- length(g.y)
      yy0 <- c(
        mean(c(ylim[1], min(w[, nt]))),
        (g.y[1:(n.g-1)] + g.y[2:n.g])/2,
        mean(c(max(w[, nt]), ylim[2])))
      y.ylabs <- unlist(
        lapply(1:length(g.y), function(j) {
        if (g.n[j]==1)
          return(g.y[j])
        if (j==1) {
          a <- yy0[1]
        } else {
          a <- seq(
            max(yy[gg==names(g.y)[j-1]]),
            min(yy[gg==names(g.y)[j]]),
            length=4)[3]
        }
        if (g.n[j+1]==1) {
          b <- yy0[j+1]
        } else {
          b <- seq(
            max(yy[gg==names(g.y)[j]]),
            min(yy[gg==names(g.y)[j+1]]),
            length=4)[2]
        }
        return(seq(a, b, length=g.n[j]))
      }))
    }
  }
  plot(x,
       w[1, ],
       xlim=xlim,
       ylim=ylim,
       type='l',
       col=col[1],
       axes=FALSE,
       xlab=xlab,
       ylab=ylab,
       ...)
  for (j in 2:nrow(w))
    lines(x, w[j,], col=col[j], ...)
  do.call('axis', c(list(side=1), axis1args))
  do.call('axis', c(list(side=2), axis2args))
  if (!is.null(y.labs))
    text(x.ylabs, y.ylabs, y.labs[oy],
         col=col[oy], adj=c(0, 0.5))
  return(invisible())
}
