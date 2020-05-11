#' @title epidplot
#' @description `epidplot` visualize data
#' from an epidemic outbreak
#' @param x a date or numeric vector
#' @param y the cummulated number of cases
#' @param which to specify which plots will be produced.
#' See details.
#' @param leg.args a named list with arguments passed to
#' the legend of the first plot.
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
#' @return up to three plots.
epidplot <-
  function(x, y, which=1:3,
           leg.args=list(
             x='topleft',
             legend=c('accumulated', 'new', 'smoothed new'),
             col=c(1, 2, 2),
             pch=c(19, 8, NA),
             lty=c(NA, NA, 1),
             bty='n'),
           lxlab=list('', '', ''),
           lylab=list('Cases', 'Absolute growth',
                      'Growth velocity'),
           ask = prod(par("mfcol")) < length(which) &&
             dev.interactive(),
           ...)
  {
    if (is.null(y)) {
      y <- x
      x <- 1:length(y)
      xl <- list(x=pretty(x, 15))
      xl$l <- xl$x
    } else {
      xl <- list(x=pretty(x, 15))
      xl$l <- format(xl$x, '%d%b')
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
    show <- 1:3 %in% which
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
    if (show[1]) {
      plot(x, log(ifelse(y==0, y0, y), logbase),
           pch=leg.args$pch[1],
           col=leg.args$col[1], axes=FALSE,
           xlab=lxlab[[1]],
           ylab=lylab[[1]], ...)
      axis(1, xl$x, xl$l, las=2)
      axis(2, yl$y, yl$l, las=1)
      points(x, log(ifelse(dy==0, y0, dy), logbase),
             pch=leg.args$pch[2],
             col=leg.args$col[2])
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
           xlab=lxlab[[2]],
           ylab=lylab[[3]], ...)
      abline(h=0, lty=2, col=gray(0.5, 0.5))
    }
    invisible()
  }
