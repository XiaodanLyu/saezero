#' Simulate unit responses
#'
#' Simulate responses given the model parameter and covariates
#' under the unit level model of Lyu, Berg and Hofmann.
#'
#' @param fit a list of model parameter
#'   estimates containing at least fixed effects coefficients and variance components
#'   (named as the return value of \code{\link{mleLBH}}).
#' @param Xs covariates matrix or data frame containing the variables named in
#'   \code{f_pos}, \code{f_zero} and \code{f_area}.
#' @param f_pos an object of class \code{\link[stats]{formula}}:
#'   a symbolic description of the fixed effect model to be fitted to the positive part.
#' @param f_zero an object of class \code{\link[stats]{formula}}:
#'   a symbolic description of the fixed effect model to be fitted to the binary part.
#'   Default value is to using the same formula as the positive part (\code{f_pos}).
#' @param f_area an object of class \code{\link[stats]{formula}}:
#'   a symbolic description of the area code to be fitted to both the positive part and the negative part.
#' @param link a specification for the link function used to model the binary part.
#'   The accepted link functions are \code{logit}, \code{probit}, \code{cauchit}, \code{log} and \code{cloglog}.
#'   Default value is "logit".
#'
#' @return response vector
#'
#' @examples
#'   erosion_2p <- as.2pdata(f_pos = RUSLE2~logR+logK+logS,
#'                           f_zero = ~logR+logS+crop2+crop3,
#'                           f_area = ~cty, data = erosion)
#'   fit <- mleLBH(erosion_2p)
#'   responses <- simLBH(fit, Xaux, f_pos = ~logR+logK+logS,
#'                       f_zero = ~logR+logS+crop2+crop3, f_area = ~cty)
#' @export
simLBH <- function(fit, Xs, f_pos, f_zero = f_pos, f_area, link = "logit"){
  ## model parameters
  beta <- fit$fixed$p1; alpha <- fit$fixed$p0
  sig2le <- fit$errorvar; sig2lu <- fit$refvar1
  sig2lb <- fit$refvar0; rho <- fit$refcor

  ## model matrices
  Xs <- as.data.frame(Xs)
  area <- Xs[, all.vars(f_area)]
  if(length(f_pos)==3) f_pos <- f_pos[-2]
  Xs1 <- model.matrix(f_pos, data = Xs)
  if(length(f_zero)==3) f_zero <- f_zero[-2]
  Xs0 <- model.matrix(f_zero, data = Xs)

  Nis <- tapply(area, area, length)
  N <- sum(Nis)
  D <- length(Nis)
  GN <- Gmat(Nis)

  us <- rnorm(D)*sqrt(sig2lu)
  luN <- GN%*%us
  leN <- rnorm(N)*sqrt(sig2le)
  lyN <- Xs1%*%beta + luN + leN

  mubu <- rho*sqrt(sig2lb/sig2lu)*us
  sbu <- sqrt((1-rho^2)*sig2lb)
  bs <- mubu + rnorm(D)*sbu

  lbN <- GN%*%bs
  mup <- Xs0%*%alpha + lbN
  p <- make.link(link)$linkinv(mup)
  oind <- rbinom(N, size = 1, prob = p)

  return(exp(lyN)*oind)
}
