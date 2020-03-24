#' Values of auxiliary variables for 64 domains
#'
#' Values of auxiliary variables for units within 64 domains of data set \code{\link{erosion}}.
#'
#' @format A data frame with 16,580 rows and 10 variables.
#' \itemize{
#'   \item{ctylab}: {county name}
#'   \item{cty}: {county code}
#'   \item{mukey}: {soil map unit key}
#'   \item{crop}: {crop category from the CDL data}
#'   \item{logR}: {\code{log} USLE rainfall factor}
#'   \item{logK}: {\code{log} USLE erosion erodibility factor}
#'   \item{logS}: {\code{log} USLE soil slope factor}
#'   \item{crop2}: {indicator of soybean}
#'   \item{crop3}: {indicator of spring wheat}
#'   \item{cnt}: {number of crop pixels in soil map unit segment within county, from the CDL data}
#' }
#' @source \url{https://lyux.shinyapps.io/viscover/}
"Xaux"

#' Simulated erosion data
#'
#' Simulated data that mimics (the main features of) the real soil erosion data and
#' other related variables for 64 South Dakota counties.
#'
#' @format A data frame with 646 rows and 10 variables.
#' \itemize{
#'   \item{ctylab}: {county name}
#'   \item{cty}: {county code}
#'   \item{mukey}: {soil map unit key}
#'   \item{crop}: {crop category from the CDL data}
#'   \item{logR}: {\code{log} USLE rainfall factor}
#'   \item{logK}: {\code{log} USLE erosion erodibility factor}
#'   \item{logS}: {\code{log} USLE soil slope factor}
#'   \item{crop2}: {indicator of soybean}
#'   \item{crop3}: {indicator of spring wheat}
#'   \item{RUSLE2}: {soil erosion loss in tons per year}
#' }
#' @seealso \code{\link{Xaux}}
"erosion"
