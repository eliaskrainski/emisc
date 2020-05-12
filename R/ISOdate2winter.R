#' Convert date to winter relative pertain
#' @description Takes any date in "POSIXct"
#' fomat and convert it to a relative number
#' to the mid of the winter period.
#' @param Date POSIXct format date
#' @param type if it is computed using a
#' cosine function or linear.
#' @return numeric in (0,1)
#' @details consider the winter season in the
#' northern hemisphere and compute a pertain
#' function defined such that it will be near 1
#' in the center of the winter period and near
#' zero in the center of the summer period.
#' Considering meteorological seasons as:
#'  Winter 1 Dec to 28(29) Feb;
#'  Spring: 1 Mar to 31 May;
#'  Summer: 1 Jun to 31 Aug,
#'  Autumn: 1 Sep to 30 Nov.
#' This definition of the seasons gives 0
#' around July 16 and 1 around January 15.
#' @examples
#' date <- ISOdate(2019, 1:12, 15)
#' ISOdate2winter(date)
#' ISOdate2winter(dae, 'linear')
ISOdate2winter <- function(Date,
                           type=c('cos', 'linear'))
  {
  type <- match.arg(type)
  if (type=='cos') {
    y0 <- round(mean(range(as.integer(substr(Date, 1, 4)))))
    d0 <- ISOdate(y0, 1, 15, 0, 0)
    dt <- as.numeric(difftime(Date, d0, units='days'))/365.25
    return(0.5*(1+cos(2*pi*dt)))
  }
  if (type=='linear') {
    y0 <- as.integer(substr(Date, 1, 4))
    d0 <- abs(difftime(Date-ISOdate(y0-1, 7, 16, 0, 0),
                       units='days'))/183
    d1 <- abs(difftime(Date-ISOdate(y0, 7, 16, 0, 0),
                       units='days'))/183
    return(pmin(as.numeric(d0), as.numeric(d1)))
  }
}
