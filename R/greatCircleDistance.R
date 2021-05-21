#' Great circle distance from longlat coordinates
#' @description given the coordinates in the longitude
#' and latitude it computes the distance between these
#' points. The default distance computed is in kilometers.
#' @param x1 longitude of the first location set
#' @param y1 latitude of the second location set
#' @param x2 longitude of the first location set
#' @param y2 latitude of the second location set
#' @param R the radius. Default is 6371 gives the
#' the resulting distance in kilometers.
#' @export
#' @examples
#'
#' ## Consider some locations
#' locs <- rbind(
#'   Paris=c(2.3295489, 48.8588377),
#'   Berlin=c(13.4072169, 52.5067614),
#'   Curitiba=c(-49.2582225, -25.4420471),
#'   Tokyo=c(139.6708224, 35.6681625))
#' locs
#'
#' apply(locs, 1, function(xy)
#'   gcDist(xy[1], xy[2], locs[,1], locs[,2]))
#'
gcDist <- function(x1, y1, x2, y2, R=6371) {
  x1 <- x1*pi/180
  x2 <- x2*pi/180
  y1 <- y1*pi/180
  y2 <- y2*pi/180
  cdl <- cos(abs(x1-x2))
  R*acos(sin(y1)*sin(y2) + cos(y1)*cos(y2)*cdl)
}
