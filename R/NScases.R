#' Confirmed cases of COVID19 in
#' Nova Scotia, Canada.
#'
#' @docType data
#'
#' @usage data(NScases)
#'
#' @format An object of class \code{"list"}
#' containing 'day' (in Date format) and
#' 'cases', the integer vector of daily counts.
#'
#' @keywords datasets
#'
#' @source \href{https://novascotia.ca/coronavirus/}{Oficial website}
#'
#' @examples
#' data(NScases)
#' with(NScases, plot(day, cases, pch=8))
#'
"NScases"
