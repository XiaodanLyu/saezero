#' Values of auxiliary variables for 64 domains
#'
#' Values of auxiliary variables for full population units within 64 domains of data set \code{\link{erosion}}.
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
#' @source \href{https://lyux.shinyapps.io/viscover/}{viscover}: a live Shiny tool
#'   featured in the RStudion Shiny gallery,
#'   which demonstrates how we integrate the USDA-NRCS Soil data and
#'   the USDA-NASS Cropland data in order to
#'   produce this dataset of auxiliary information for the cropland population in South Dakota.
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

#' South Dakota map data
#'
#' South Dakota county map coordinates in longitude and latitude with FIPS.
#'
"map_sd"
