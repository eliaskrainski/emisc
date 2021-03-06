% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotfitt.R
\name{plotfitt}
\alias{plotfitt}
\title{Visualize fit output against observed values}
\usage{
plotfitt(
  o,
  sfitt,
  asp = 1,
  xlab = "Fitted",
  ylab = "Observed",
  xlim,
  ylim,
  ...
)
}
\arguments{
\item{o}{vector of observed values}

\item{sfitt}{an output of \code{lm}, \code{glm} or
a \code{data.frame} containing columns with
'mean', 'sd', lower limit IC, median and upper limit IC.}

\item{asp}{the aspect ratio}

\item{xlab}{the x axis title}

\item{ylab}{the y axis title}

\item{xlim}{the limits for the x axis scale}

\item{ylim}{the limits for the y axis scale}

\item{...}{additional arguments passed to
\code{plot} and \code{arrows}}
}
\description{
Visualize fit output against observed values
}
\examples{
d <- list(x=seq(pi, 3*pi, 0.1))
d$y <- rpois(length(d$x), exp(1 + d$x))
r <- glm(y~x, poisson, d)
plotfitt(d$y, r)
}
