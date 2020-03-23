#' MLE model parameter estimator
#'
#' Fits by ML method the unit level model of Lyu, Berg and Hofmann to a transformation of
#' the specified response variable. The specified link function is used in the binary part.
#'
#' @param f_pos an object of class \code{\link[stats]{formula}}:
#'   a symbolic description of the model to be fitted to the positive part.
#' @param f_zero an object of class \code{\link[stats]{formula}}:
#'   a symbolic description of the model to be fitted to the binary part.
#' @param area vector with area codes.
#' @param data optional data frame containing the variables named in \code{f_pos}, \code{f_zero} and \code{area}.
#' @param link a specification for the link function used to model the binary part.
#'   The accepted link functions are \code{logit}, \code{probit}, \code{cauchit}, \code{log} and \code{cloglog}.
#'
#' @return The function returns a list with the following objects:
#' \itemize{
#'   \item \code{fixed}: list with the estimated values of the fixed regression coefficient
#'   in the positive part (\code{p1}) and in the binary part (\code{p0}).
#'   \item \code{random}: data frame with the predicted random effects in the positive part (\code{p1})
#'   and the binary part (\code{p0}).
#'   \item \code{errorvar}: estimated model error variance in the positive part.
#'   \item \code{refvar1}: estimated random effects variance in the positive part.
#'   \item \code{refvar0}: estimated random effects variance in the binary part.
#'   \item \code{refcor}: estimated correlation coefficient of the random effects between the two parts.
#'   \item \code{loglike}: log-likelihood.
#'   \item \code{residuals}: data frame with the raw residuals from the model fit in the positive part
#'   (\code{p1}, cases with zero response are \code{NA}.) and the binary part (\code{p0}).
#' }
#'@export
mleLBH <- function(f_pos, f_zero, area, data, link = "logit"){

}
