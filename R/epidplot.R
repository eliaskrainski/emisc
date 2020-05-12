#' Visualize an epdemic outbreak
#' @description visualize the cummulated time
#' series and its derivative (news).
#' Add smoothed lines of the news.
#' Additional plots shows the 1st and 2nd
#' derivatives of the smoothed, interpreted
#' as the absolute increase and velocity.
#' @param x a date or numeric vector
#' @param y the cummulated number of cases
#' @param which integer vector to pecify which
#' plots will be produced. See details.
#' @param w numeric vector to compute the expected
#' number of new cases given past cases. See details.
#' @param leg.args a named list with arguments
#' passed to the legend of the first plot.
#' @param lxlab a list with the xlab for each plot.
#' @param lylab a list with the xlab for each plot.
#' @param ask logical if ‘TRUE’, the user is _ask_ed
#' before each plot, see \code{par(ask=.)}.
#' @param ... additional arguments passed to \code{hist}
#' @details
#' The first plot shows the series of accumulated cases,
#' the number of new cases along with a line which is the
#' smoothed version of the new cases counts.
#' The smoothed is computed using the \code{mgcv:::gam}
#' function with \code{mgcv:::s()} on \code{x} and
#' \code{family=poisson()}.
#' The second plot is the, numerical, derivative
#' of the fitted in the previous plot.
#' The third plot is the, numerical, derivative
#' of the fist derivative in the previous plot.
#' The fourth plot is a smoothed version of R_t.
#' This is computed as R_t = s_t/E_t, where
#' s_t is smoothed version of the (new) series,
#' E_t = sum_j w_j s_{t-j}.
#' @return up to four plots.
#' @examples
#' ns.cases <- c(5, 7, 12, 14, 15, 21, 28, 41,
#'  51, 68, 73, 90, 110, 122, 127, 147, 173,
#'  193, 207, 236, 262, 293, 310, 310, 342,
#'  407, 428, 445, 474, 517, 549, 579, 606,
#'  649, 675, 721, 737, 772, 827, 850, 865,
#'  873, 900, 915, 935, 947, 959, 963, 971,
#'  985, 991, 998, 1007, 1008, 1011, 1018, 1019)
#' date <- as.Date('2020-03-15') + 1:length(ns.cases)
#' pw <- pgamma(0:14, shape=(5/3)^2, scale=3^2/5)
#' w <- diff(pw)/sum(diff(pw))
#' par(mfrow=c(2,2), mar=c(2,3,1,1),
#'     mgp=c(2, 0.5, 0), las=1)
#' epidplot(date, ns.cases, w=w)
epidplot <-
  function(x, y, which=1:4, w,
           leg.args=list(
             x='topleft',
             legend=c('accumulated', 'new',
                      'smoothed new'),
             col=c(1, 2, 2),
             pch=c(19, 8, NA),
             lty=c(NA, NA, 1),
             bty='n'),
           lxlab=list('', '', '', ''),
           lylab=list('Cases', 'Absolute growth',
                      'Growth velocity',
                      'Reproduction number'),
           ask = prod(par("mfcol")) < length(which) &&
             dev.interactive(),
           ...)
  {
    if (is.null(y)) {
      y <- x
      x <- 1:length(y)
      xl <- list(x=pretty(x))
      xl$l <- xl$x
    } else {
      xl <- list(x=pretty(x))
      xl$l <- format(xl$x, '%b %d')
    }
    y <- cummax(y)
    logbase <- 10; y0 <- 0.3
    yl0 <- c(1, 3)
    k <- ceiling(log(logbase + 1 + max(y), logbase))
    kk <- rep(0:k, each=length(yl0))
    yl <- list(y=c(log(y0, logbase),
                   rep(log(yl0, logbase), k+1) + kk))
    yl$l <- c(0, paste0(yl0*(logbase^(kk%%3)),
                        c('', 'K', 'M', 'B')[(kk%/%3)+1]))
    dy <- diff(c(0, y))
    ff <- mgcv::gam(n ~ s(d), poisson(),
                    data=list(n=dy, d=as.numeric(x))
    )$fitted
    df1 <- c(NA, diff(ff))
    df2 <- c(NA, diff(df1))
    logbase <- 10
    nwin <- prod(par('mfcol'))==1
    if (ask) {
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))
    }
    show <- 1:4 %in% which
    if (show[1]) {
      plot(x,
           log(ifelse(y==0, y0, y), logbase),
           ylim=log(c(y0, max(y)), logbase),
           pch=leg.args$pch[1],
           col=leg.args$col[1], axes=FALSE,
           xlab=lxlab[[1]],
           ylab=lylab[[1]], ...)
      axis(1, xl$x, xl$l)
      axis(2, yl$y, yl$l)
      points(x, log(ifelse(dy==0, y0, dy), logbase),
             pch=leg.args$pch[2],
             col=leg.args$col[2], ...)
      lines(x, log(ff, logbase),
            lwd=leg.args$lwd[3],
            lty=leg.args$lty[3],
            col=leg.args$col[3])
      do.call('legend', leg.args)
    }
    if (show[2]) {
      plot(x, df1, type='l',
           ylim=range(0, df1, na.rm=TRUE),
           xlab=lxlab[[2]],
           ylab=lylab[[2]], ...)
      abline(h=0, lty=2, col=gray(0.5, 0.5))
    }
    if (show[3]) {
      plot(x, df2, type='l',
           ylim=range(0, df2, na.rm=TRUE),
           xlab=lxlab[[3]],
           ylab=lylab[[3]], ...)
      abline(h=0, lty=2, col=gray(0.5, 0.5))
    }
    if (show[4]) {
      if (missing(w))
        stop("To compute R_t 'w' must be provided!")
      n <- length(ff)
      k <- length(w)
      Rs <- Rt <- es <- ee <- rep(NA, n)
      for (i in (k+1):n) {
        ee[i] <- sum(dy[i-1:k] * w)
        Rt[i] <- dy[i]/ee[i]
        es[i] <- sum(ff[i-1:k] * w)
        Rs[i] <- ff[i]/ee[i]
      }
      plot(x, Rt, pch=8,
           ylim=range(0:1, Rt, Rs, na.rm=TRUE),
           xlab=lxlab[[4]],
           ylab=lylab[[4]], ...)
      lines(x, Rs, ...)
      abline(h=1, lty=2, col=gray(0.3, 0.5), lwd=2)
    }
    invisible()
  }
