#' EB estimator of area means and associated MSE estimator
#'
#' Obtains numerical approximations of EB estimators of area means
#' under the unit level model of Lyu, Berg and Hofmann
#' when the values of auxiliary variables for population units and
#' the model parameter estimates are available.
#'
#' @param Xaux matrix or data frame containing covariates, the area code and
#'   the variables named in \code{f_q} for population units.
#' @param f_q an object of class \code{\link[stats]{formula}}:
#'   a symbolic description of the number of population units with the same covariates.
#'   Default value is \code{~1}.
#' @param data_2p a two-part data object returned by \code{\link{as.2pdata}}.
#' @param fit a list of model parameter
#'   estimates containing at least fixed effects coefficients and variance components
#'   (named as the return value of \code{\link{mleLBH}}).
#' @param fullpop a boolean variable indicating whether \code{Xaux} contains covariates information
#'   for the full population (\code{TRUE}) or just the out-of-sample units (\code{FALSE}).
#'   Default value is \code{FALSE}. The details of this indicator are given under Details.
#'
#' @details
#'   When \code{Xaux} contains only the covariates of the out-of-sample units (\code{fullpop = FALSE}),
#'   observed response is used for the sampled units when calculating the EB estimator.
#'   When \code{Xaux} contains the covariates of the full population (\code{fullpop = TRUE}),
#'   unit-level EB prediction is used for the sampled units. This is
#'   reasonable when the sampling fraction is extremely small in each area
#'   (e.g., \href{https://doi.org/10.1080/01621459.1988.10478561}{Battese, Harter and Fuller (1988)}).
#'
#' @return A data frame with the number of rows equal to the number of unique areas in \code{Xaux}:
#' \itemize{
#'  \item \code{area}: area codes
#'  \item \code{eb}: EB estimator of area means
#'  \item \code{mse}: the One-step MSE estimator
#' }
#'
#' @seealso \code{\link{as.2pdata}}, \code{\link{mleLBH}}
#' @examples
#'   erosion_2p <- as.2pdata(f_pos = RUSLE2~logR+logK+logS,
#'                           f_zero = ~logR+logS+crop2+crop3,
#'                           f_area = ~cty, data = erosion)
#'   fit <- mleLBH(erosion_2p)
#'   predictions <- ebLBH(Xaux, f_q = ~cnt, erosion_2p, fit, fullpop = TRUE)
#' @export
#'
ebLBH <- function(Xaux, f_q = ~1, data_2p, fit, fullpop = FALSE){

  ## argument parsing and checking ####
  Xaux <- data.frame(Xaux)
  if (length(all.vars(f_q))) {
    q <- Xaux[, all.vars(f_q)]
  } else {
    q <- rep(1, nrow(Xaux))
  }
  fs <- attributes(data_2p)
  area_oos <- Xaux[, all.vars(fs$f_area)]
  Xs1_oos <- model.matrix(fs$f_pos[-2], data = Xaux)
  Xs0_oos <- model.matrix(fs$f_zero, data = Xaux)
  link <- attributes(fit)$link
  if(class(link)!="link-glm") link <- make.link(link)

  ## only suitable for log-transformed positive responses
  if(attributes(data_2p)$lambda!=0){
    stop("This method only applies to log-transformed positive responses.")
  }

  attach(data_2p, warn.conflicts = FALSE)
  attach(fit, warn.conflicts = FALSE)

  nis <- tapply(deltas, area, length)
  mis <- tapply(area_oos, area_oos, length)
  Gs <- Gmat(nis)
  Gns <- Gmat(mis)

  if (fullpop){
    fis <- 0
  } else{
    fis <- nis
  }

  ## parameter estimates
  beta <- fixed$p1; alpha <- fixed$p0
  sig2le <- errorvar; sig2lu <- refvar1
  sig2lb <- refvar0; rho <- refcor
  ## conditional distributions ####
  D <- length(nis)
  nti <- tapply(deltas, area, sum)
  rt <- rep(0, length(area))
  rt[deltas==1] <- lys - Xs1%*%beta

  gamma <- (1-rho^2)*sig2lu/((1-rho^2)*sig2lu+sig2le/nti)
  rbar <- drop(t(Gs)%*%rt)/nti
  rbar[nti==0] <- 0
  mu_u <- gamma*rbar
  mu_u[nti==0] <- 0
  eta <- exp((1-gamma)*rho*sqrt(sig2lu)/sqrt(sig2lb))
  sig2_u <- gamma*sig2le/nti
  sig2_u[nti==0] <- (1-rho^2)*sig2lu
  b_m <- rho*sqrt(sig2lu*sig2lb)*rbar/(sig2lu+sig2le/nti)
  b_v <- (1-rho^2)/(1-(1-gamma)*rho^2)*sig2lb

  ## Gauss-Hermite method
  ghwts <- do.call("rbind", lapply(1:D, function(i){
    wts <- statmod::gauss.quad.prob(20, dist = "normal", mu = b_m[i], sigma = sqrt(b_v[i]))
    as.data.frame(wts)
  }))
  bmat <- matrix(ghwts$nodes, nrow = D, ncol = 20, byrow = T)
  wmat <- matrix(ghwts$weights, nrow = D, ncol = 20, byrow = T)
  thres.gen <- function(b){
    mu_p <- Xs0%*%alpha + Gs%*%b
    p <- link$linkinv(mu_p)
    pq <- as.vector(ifelse(deltas, p, 1-p))
    thres <- apply(Gs*pq, 2, function(x) prod(x[x>0]))
    thres
  }
  den.eb <- rowSums(apply(bmat, 2, thres.gen)*wmat)

  ## EB estimator ####
  num.emp <- function(b){
    thres <- thres.gen(b)
    mu_p <- Xs0_oos%*%alpha + Gns%*%b
    pb <- link$linkinv(mu_p)
    num <- Gns%*%(thres*(eta^b))*pb
    drop(num)
  }
  num.eb <- rowSums(apply(bmat, 2, num.emp)*(Gns%*%wmat))
  eb.peta <- num.eb/(Gns%*%den.eb)
  cij <- exp(Xs1_oos%*%beta + Gns%*%(mu_u+sig2_u/2) + sig2le/2)
  nbaris <- tapply(q, area_oos, sum)
  ybar.r <- drop(t(Gns) %*% (eb.peta*cij*q))/nbaris
  ybar.s <- drop(t(Gs[deltas==1,]) %*% exp(lys))/nis
  eb <- drop(fis*ybar.s+nbaris*ybar.r)/(fis+nbaris)

  ## MSE estimator ####
  num.y2eb <- function(b){
    mu_p <- Xs0_oos%*%alpha + Gns%*%b
    pb <- link$linkinv(mu_p)
    V <- pb*cij*q
    Wij <- function(i){
      id <- which(area_oos==unique(area_oos)[i])
      W <- V[id,] %*% t(V[id,]) * exp(sig2_u[i]) * eta[i]^(2*b[i])
      diag(W) <- diag(W)*(1 + (exp(sig2le)/pb[id]-1)/q[id])
      diag(W)[q[id]==0] <- 0
      sum(W)
    }
    Wi <- sapply(1:D, Wij)
    thres.gen(b)*Wi
  }
  num.vareb <- rowSums(apply(bmat, 2, num.y2eb)*wmat)
  ybar2.r <- (1/nbaris)^2*num.vareb/den.eb
  a <- nbaris/(nbaris+fis)
  mse <- a^2*(ybar2.r-(ybar.r)^2)

  detach(data_2p)
  detach(fit)
  return(data.frame(area = unique(area_oos), eb = eb, mse = mse))
}
