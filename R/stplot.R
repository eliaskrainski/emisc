#' Time series plot for each location
#' @description given a matrix where each line is a
#' time series of each spatial unit, the time series
#' are shown over the map
#' @param x a matrix with each column containing data
#' at each time for a set of geographical units
#' @param sp an object as defined in the 'sp' package
#' @param d a vector to specify a \code{GridTopology}
#' object to do some subset over the \code{sp} object.
#' @param col a vector of colors for each time series.
#' If not provided, the collors will be defined considering
#' the relative growth of each time series and series with
#' similar growth will have similar colors. By default
#' red color is for the series with most positive growth
#' and blue for those with least growth.
#' @param ce vector of coordinates to center the time
#' the time series. If not provided, the \code{coordinates}
#' output applied to \code{sp}.
#' @param tsub if NULL, default, the colors consider
#' the average of the time series. If integer,
#' specifies the number of temporal subdivisions
#' to evaluate the growth of each time series,
#' used to define the colors.
#' @param ex is a relative expansion factor for
#' the x-axis extent
#' @param ey is a relative expansion factor for
#' the y-axis extent
#' @param add logical indicating if the plot is to
#' be added to a current plot
#' @param legend logical indicatint if the legend
#' is bo be drawn
#' @param leg.args legend arguments, used if legend=TRUE
#' @param verbose logical to specify if some steps
#' are being reported during the process work
#' @param ... additional arguments passed to
#' \code{line} used to plot each time series.
#' @export
#' @examples
#'\dontrun{
#' ## Consider data from the surveillance package
#' data(measlesWeserEms, package='surveillance')
#' if (any(ls()=='measlesWeserEms'))
#'   stplot(t(measlesWeserEms@observed),
#'          measlesWeserEms@map)
#' }
stplot <- function(x, sp, d, col, ce,
                   tsub=NULL, ex=1, ey=1,
                   add=FALSE, legend=TRUE,
                   leg.args=list(x='bottomleft', lty=1),
                   verbose=FALSE, ...) {
  n <- nrow(x)
  p <- ncol(x)
  b <- bbox(sp)
  r <- apply(b, 1, diff)
  if (legend) {
      if (!any(names(leg.args)=='legend')) {
          rrxx <- apply(x, 1, range, na.rm=TRUE)
          leg.args$legend <- paste0(
              rownames(x), ' ',
              format(rrxx[1,], digits = 1),
              ' : ',
              format(rrxx[2,], digits = 1)
          )[rev(order(rrxx[2,]-rrxx[1,]))]
      }
  }
  if (missing(col)) {
    if (length(tsub)>0) {
      mx <- mean(x, na.rm=TRUE)
      jsub <- sort(((0:(ncol(x)-1))%%tsub) + 1)
      usub <- unique(jsub)
      ccod <- (mx<rowMeans(x[, jsub==usub[1]],
                           na.rm=TRUE))+0
      if (length(usub)>1)
        for (j in usub[-1])
          ccod <- ccod +
        (mx<rowMeans(x[, jsub==usub[j]],
                     na.rm=TRUE))*(10^(j-1))
      ox <- pmatch(ccod, unique(sort(ccod)),
                   duplicates.ok=TRUE)
      if (verbose) print(table(ox))
      ox.ok <- !is.na(ox)
      nc <- max(ox, na.rm=TRUE)
      col <- rep('transparent', n)
      col[ox.ok] <- rgb(
        ox[ox.ok]/nc,
        1-2*abs(ox[ox.ok]-(0.5+nc/2))/nc,
        1-ox[ox.ok]/nc)
    } else {
      col <- x2rgb(rowMeans(x, na.rm=TRUE), u=TRUE)
      if (verbose) cat('Defined colors as function of mean.\n')
    }
  }
  if (legend)
      if (!any(names(leg.args)=='col'))
          leg.args$col <- col
  yl <- range(x, na.rm=TRUE)
  ry <- diff(yl)
  if (missing(d)) {
      if (class(sp)%in%
          c('SpatialPolygons',
            'SpatialPolygonsDataFrame')) {
          if (missing(ce))
              ce <- coordinates(sp)
          if (!add)
              sp::plot(sp, border=gray(.5))
          ss <- sqrt(mean(sapply(
              1:length(sp@polygons), function(i)
                  sp@polygons[[i]]@area)))*0.5
          ssxy <- rowMeans(sapply(1:n, function(i)
              apply(bbox(sp[i,]), 1, diff)))
          xi0 <- seq(-0.5, 0.5, length=p)
          for (i in 1:n) {
              b <- bbox(sp[i,])
              xi <- xi0*ex*ssxy[1] + mean(b[1,])
              xi <- xi-mean(xi) + ce[i, 1]
              yi <- ssxy[2]*((x[i,]-yl[1])/ry - 0.5)*ey +
                  mean(b[2,])
              yi <- yi-mean(yi)+ce[i, 2]
              lines(xi, yi, col=col[i], ...)
              segments(xi[1], mean(yi), xi[length(xi)], mean(yi), lty=3)
          }
      } else {
          if (class(sp)%in%
              c('SpatialPoints',
                'SpatialPointsDataFrame')) {
              ssxy <- apply(bbox(sp), 1, diff)/sqrt(n)
              if (!add)
                  sp::plot(sp, cex=0.01, ...)
              xi0 <- seq(-0.5, 0.5, length=p)
              for (i in 1:n) {
                  xi <- xi0 * ex * ssxy[1] + sp@coords[i, 1]
                  yi <- ssxy[2] * ((x[i,] - yl[1])/ry - 0.5)*ey +
                      sp@coords[i, 2]
                  lines(xi, yi, col=col[i], ...)
                  segments(xi[1], mean(yi), xi[length(xi)], mean(yi), lty=3)
              }
          }
      }
  } else {
    g <- GridTopology(b[,1]+r/(2*d), r/d, d)
    G <- SpatialGrid(g, sp@proj4string)
    isg <- over(SpatialPoints(
      coordinates(sp), sp@proj4string), G)
    G <- as.SpatialPolygons.GridTopology(g)
    par(mfrow=rev(d), mar=c(0,0,0,0),
        xaxs='i', yaxs='i')
    for (i in 1:prod(d)) {
      jj <- which(i==isg)
      gi <- G[i, ]
      if (!add)
          sp::plot(gi, asp=0, border=gray(.7), lty=2, ...)
      sp::plot(sp, add=TRUE, border=gray(.3))
      xyll <- bbox(gi)[c(1,3,2,4)]
      xx <- seq(xyll[1], xyll[2], length=p)
      if (length(jj)>0) {
        for (j in jj) {
          yy <- xyll[3] + (xyll[4]-xyll[3])*
            (x[j,]-yl[1])/diff(yl)
          lines(xx, yy,
                col=col[which(j==jj)])
          segments(xx[1], mean(yy), xx[length(xx)], mean(yy), lty=3)
        }
      }
    }
  }
  if (legend)
      do.call('legend', leg.args)
  invisible()
}
