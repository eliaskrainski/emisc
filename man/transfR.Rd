% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transfR.R
\name{transfR}
\alias{transfR}
\title{Function for Real to Real transform
of Real to Positive functions
and its inverse transformation}
\usage{
transfR(x, transf, inverse = FALSE, base = 2)
}
\arguments{
\item{x}{numeric}

\item{transf}{character name of the function
that implements the transformation to be applied.}

\item{inverse}{logical indicating if the inverse
of the transformation is to be performed}

\item{base}{numeric to be considered as the
base of the log transformation}
}
\value{
numeric with the transformed input
}
\description{
Compute the the desired transformation
for the abs of the input and multiply
the output by the sign of the input.
}
\examples{
transfR(-5:5, 'sqrt')
transfR(-5:5, 'log')
plot(function(x) transfR(x, 'log'), -10, 10)
plot(function(x) transfR(x, 'log', base=10),
     -10, 10, col=4, add=TRUE)
x <- -30:30/30
plot(x, transfR(transfR(x, 'log'),
     'log', inverse=TRUE))
}
