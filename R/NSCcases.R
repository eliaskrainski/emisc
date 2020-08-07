#' Confirmed cases of COVID19 in
#' Nova Scotia, Canada.
#'
#' @docType data
#'
#' @usage data(NSCcases)
#'
#' @format An object of class \code{"list"}
#' containing 'day' (in Date format) and
#' 'cases', the integer vector of daily counts.
#'
#' @keywords datasets
#'
#' @source \href{https://novascotia.ca/coronavirus/}
#'
#' @examples
#' data(NSCcases)
#' with(NSCcases, plot(day, cases, pch=8))
#'
"NSCcases"
