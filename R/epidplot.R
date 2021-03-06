#' Visualize an epdemic outbreak
#' @description visualize the cummulated time
#' series and its derivative (news).
#' Add smoothed lines of the news.
#' Additional plots shows the 1st and 2nd
#' derivatives of the smoothed, interpreted
#' as the absolute increase and velocity.
#' @param x a date or numeric vector
#' @param y the cummulated number of cases
#' @param log10 length one or two logical to
#' indicate if accumulated and new series are
#' shown in the log10 transformation.
#' @param which integer vector or list of integers
#' to specify which plots will be produced in each
#' figure. See detais.
#' @param w numeric vector to compute the expected
#' number of new cases given past cases. See details.
#' @param leg.args a named list with arguments
#' passed to the legend of the first plot.
#' @param lxlab a list with the xlab for each plot.
#' @param lylab a list with the xlab for each plot.
#' @param ask logical if ‘TRUE’, the user is _ask_ed
#' before each plot, see \code{par(ask=.)}.
#' @param ... additional arguments passed to \code{hist}
#' @export
#' @details
#' This function produces up to 6 plots.
#' 1st: accumulated, provided y, with smoothed line
#' 2nd: new cases with smoothed line
#' 3nd: derivative, change on the trend, of the
#' smoothed line for new cases.
#' 4th: derivative divided by the smoothed new cases line
#' 5th: second derivative of the smoothed line for the
#' new cases series, representing the velocity of
#' the absolute change.
#' 6th: effective reproduction number.
#' The first and second plots can be joined into one.
#' The smoothed line for the second plot is computed
#' using the \code{mgcv::gam} function with
#' \code{mgcv::s()} on \code{x} considering
#' \code{family=poisson()}.
#' The line for the 1st plot is the cumsum of the
#' line being fitted for the 2nd plot.
#' The last plot consider
#'     R_t = s_t/E_t
#'   where E_t = sum_j w_j s_{t-j},
#' for j = 1, ..., k=length(w).
#' In this plot it will be shown the last (n-k)
#' values for R_t, if type!='n'.
#' @return up to six plots
#' @examples
#' ## COVID19 cases in Nova Scotia, Canada
#' ns.cases <- c(5, 7, 12, 14, 15, 21, 28, 41,
#'  51, 68, 73, 90, 110, 122, 127, 147, 173,
#'  193, 207, 236, 262, 293, 310, 310, 342,
#'  407, 428, 445, 474, 517, 549, 579, 606,
#'  649, 675, 721, 737, 772, 827, 850, 865,
#'  873, 900, 915, 935, 947, 959, 963, 971,
#'  985, 991, 998, 1007, 1008, 1011, 1018,
#'  1019, 1020, 1024, 1026, 1034, 1037, 1040,
#'  1043, 1044, 1045, 1046, 1048, 1049, 1050,
#'  1051, 1052, 1053, 1055, 1055, 1056, 1056,
#'  1057, 1057, 1058, 1058, 1058, 1058, 1059,
#'  1059, 1060, 1061, 1061, 1061, 1061)
#' date <- as.Date('2020-03-15') + 1:length(ns.cases)
#' pw <- pgamma(0:14, shape=(5/3)^2, scale=3^2/5)
#' w <- diff(pw)/sum(diff(pw))
#' par(mfrow=c(2,2), mar=c(2,3,1,1),
#'     mgp=c(2, 0.5, 0), las=1)
#' epidplot(date, ns.cases, logbase=10, w=w)
epidplot <-
  function(x, y, log10 = TRUE,
           which = list(1:2, c(3,5), 4, 6),
           w=exp(-sqrt(1:14))*0.686,
           leg.args=NULL,
           lxlab=NULL,
           lylab=NULL,
           ask = prod(par("mfcol")) < length(which) &&
             dev.interactive(),
           ...)
    {
    if (is.null(y)) {
      y <- x
      x <- 1:length(y)
    }
    xl <- list(x=pretty(x))
    if (class(x)%in%c('Date', 'POSIXct', 'POSIXt')) {
      xl$l <- format(xl$x, '%b %d')
    } else {
      xl$l <- xl$x
    }
    dolog <- rep(log10, 2)[1:2]
    logbase <- 10 ##logbase>(1+sqrt(.Machine$double.eps))
    if (any(dolog)) {
      y0 <- 0.3
      yl0 <- c(1, 3)
      k <- ceiling(log(logbase + 1 + max(y), logbase))
      kk <- rep(0:k, each=length(yl0))
      yl <- list(y=c(log(y0, logbase),
                     rep(log(yl0, logbase), k+1) + kk))
      yl$l <- c(0, paste0(yl0*(logbase^(kk%%3)),
                          c('', 'K', 'M', 'B')[(kk%/%3)+1]))
    } else {
      yl <- list(y=pretty(y), l=pretty(y))
    }
    dy <- diff(c(0, y))
    dyr <- dy/y
    ff <- mgcv::gam(y ~ s(i), poisson(),
                    data=list(y=dy,
                              i=as.numeric(x))
                    )$fitted
    cff <- cumsum(ff)
    df1 <- c(NA, diff(ff))
    df1r <- df1/ff
    df2 <- c(NA, diff(df1))
    nwin <- prod(par('mfcol'))
    if (ask) {
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))
    }
    for (wp in 1:length(which)) {
      show <- 1:6 %in% which[[wp]]
      if (show[1]) {
        if (is.null(leg.args[[wp]])) {
          leg.args[[wp]]$x <- 'topleft'
          leg.args[[wp]]$legend <- c('accumulated')
          leg.args[[wp]]$pch <- c(19)
          leg.args[[wp]]$col <- c(1)
          leg.args[[wp]]$lty <- c(1)
          leg.args[[wp]]$lwd <- c(1)
        }
        if (dolog[1]) {
          yplot <- log(ifelse(y==0, y0, y), logbase)
          ylm <- log(c(y0, max(y, na.rm=TRUE)), logbase)
          cffplot <- log(cff, logbase)
        } else {
          yplot <- y
          ylm <- c(0, max(y, na.rm=TRUE))
          cffplot <- cff
        }
        plot(x, yplot, ylim=ylm, axes=FALSE,
             pch=leg.args[[wp]]$pch[1],
             col=leg.args[[wp]]$col[1],
             xlab=ifelse(is.null(lxlab[[wp]]),
                         '', lxlab[[wp]]),
             ylab=ifelse(is.null(lylab[[wp]]),
                         '', lylab[[wp]]),
             ...)
        axis(1, xl$x, xl$l)
        axis(2, yl$y, yl$l)
        lines(x, cffplot)
        if (show[2]) {
          leg.args[[wp]]$legend <- c(
            leg.args[[wp]]$legend, 'new')
          leg.args[[wp]]$pch <- c(
            leg.args[[wp]]$pch, 8)
          leg.args[[wp]]$col <- c(
            leg.args[[wp]]$col, 2)
          leg.args[[wp]]$lty <- c(
            leg.args[[wp]]$lty, 1)
          leg.args[[wp]]$lwd <- c(
            leg.args[[wp]]$lwd, 1)
          if (dolog[2]) {
            dyplot <- log(ifelse(dy<y0, y0, dy), logbase)
            ylm <- c(y0, max(dy, na.rm=TRUE))
            ffplot <- log(ff, logbase)
          } else {
            dyplot <- dy
            ylm <- c(y0, max(dyplot, na.rm=TRUE))
            ffplot <- ff
          }
          points(x, dyplot,
                 pch=leg.args[[wp]]$pch[2],
                 col=leg.args[[wp]]$col[2],
                 ...)
          lines(x, ffplot,
                lty=leg.args[[wp]]$lty[2],
                lwd=leg.args[[wp]]$lwd[2],
                col=leg.args[[wp]]$col[2])
        }
        if (show[3]) {
          leg.args[[wp]]$legend <- c(
            leg.args[[wp]]$legend, 'change')
          leg.args[[wp]]$pch <- c(
            leg.args[[wp]]$pch, NA)
          leg.args[[wp]]$col <- c(
            leg.args[[wp]]$col, 3)
          leg.args[[wp]]$lty <- c(
            leg.args[[wp]]$lty, 1)
          leg.args[[wp]]$lwd <- c(
            leg.args[[wp]]$lwd, 1)
          lines(x, df1,
                leg.args[[wp]]$lty[3],
                leg.args[[wp]]$lwd[3],
                leg.args[[wp]]$col[3],
                xlab=lxlab[[wp]],
                 ylab=lylab[[wp]], ...)
          abline(h=0, lty=2, col=gray(0.5, 0.5))
        }
        if (any(show[c(2,3)]))
          do.call('legend', leg.args[[wp]])
      }
      if ((!show[1]) & show[2]) {
        if (dolog[2]) {
          dyplot <- log(ifelse(dy<y0, y0, dy), logbase)
          ylm <- log(c(y0, max(dy, na.rm=TRUE)), logbase)
          ffplot <- log(ff, logbase)
        } else {
          dyplot <- dy
          ylm <- c(0, max(dyplot, na.rm=TRUE))
          ffplot <- ff
          yl <- list(y=pretty(dyplot), l=pretty(dyplot))
        }
        plot(x, dyplot, ylim=ylm,
             pch=19, axes=FALSE,
             xlab=lxlab[[wp]],
             ylab=lylab[[wp]], ...)
        axis(1, xl$x, xl$l)
        axis(2, yl$y, yl$l)
        lines(x, ffplot, ...)
      }
      if (show[3]) {
        plot(x, df1, type='l',
             ylim=range(0, df1, na.rm=TRUE),
             xlab=lxlab[[wp]],
             ylab=lylab[[wp]], ...)
        abline(h=0, lty=2, col=gray(0.5, 0.5))
      }
      if (show[4]) {
        plot(x, df1r, type='l',
             ylim=range(0, df1r, na.rm=TRUE),
             xlab=lxlab[[wp]],
             ylab=lylab[[wp]], ...)
        abline(h=0, lty=2, col=gray(0.5, 0.5))
      }
      if (show[5]) {
        plot(x, df2, type='l',
             ylim=range(0, df2, na.rm=TRUE),
             xlab=lxlab[[wp]],
             ylab=lylab[[wp]], ...)
        abline(h=0, lty=2, col=gray(0.5, 0.5))
      }
      if (show[6]) {
        if (missing(w)) {
          warning("Missing 'w', defined 14 exponential weights!")
          w <- exp(-sqrt(1:14))
          w <- w/sum(w)
          cat('Its four decimal digits are:\n')
          print(round(1e4*w))
        }
        n <- length(ff)
        k <- length(w)
        if ((n-k)<k)
          stop(paste("'length(y)<2*length(w)':",
                     "Too few data to fit R_t"))
        ee <- rep(0.01, n+k)
        for (i in 1:n)
          ee[i+1:k] <- ee[i+1:k] + dy[i] * w
        lee <- log(ifelse(ee<0.1, 0.1, ee))
        rt.hat <- c(NA, dy[2:n]/exp(lee[2:n]))
        dtemp <- list(y=dy[2:n], i=2:n)
        m2rt <- mgcv::gam(
          y ~ s(i), poisson(),
          data=dtemp,
          offset=lee[2:n])
        prd <- predict(m2rt, se.fit=TRUE)
        plot(x, c(rep(NA, k), rt.hat[(k+1):n]), pch=8,
             xlab=lxlab[[wp]], ylab=lylab[[wp]], ...)
        m <- c(NA, exp(prd$fit))
        lines(x[(k+1):n], m[(k+1):n])
        lo <- exp(prd$fit-1.96*prd$se.fit)
        up <- exp(prd$fit+1.96*prd$se.fit)
        polygon(c(x[(k+1):n], x[n:(k+1)], x[k+1]),
                c(lo[k:(n-1)], up[(n-1):k], lo[k]),
                col=gray(0.7,0.5),
                border=gray(0.7, 0.5))
        abline(h=1)
        axis(2)
        axis(1, pmatch(xl$x, x), xl$l)
      }
  }
   invisible()
}
