#' EB estimators of area means with non-sample values of auxiliary variables
#'
#' Fits by ML method the unit level model of Lyu, Berg and Hofmann and
#' obtains numerical approximations of EB estimators of the specified area means,
#' when the values of auxiliary variables for out-of-sample units are available.
#'
#' @param Xnonsample matrix or data frame containing covariates, the area code and
#'   the variables named in \code{f_q} for the out-of-sample units.
#' @param f_q an object of class \code{\link[stats]{formula}}:
#'   a symbolic description of the number of out-of-sample units with the same covariates.
#'   Default value is \code{~1}.
#' @param data_2p a two-part data object returned by \code{\link{as.2pdata}}.
#' @param fit a list from model fitting with components as returned by \code{\link{mleLBH}}.
#'
#' @return A data frame with the number of rows equal to the number of unique areas in \code{Xnonsample}:
#' \describe{
#'  \item{\code{area}}{area codes}
#'  \item{\code{eb}}{EB estimator of area means}
#' }
#'
#' @examples
#'   erosion_2p <- as.2pdata(f_pos = RUSLE2~logR+logK+logS,
#'                           f_zero = ~logR+logS+crop2+crop3,
#'                           f_area = ~cty, data = erosion)
#'   fit <- mleLBH(erosion_2p)
#'   predictions <- ebLBH(Xaux, f_q = ~cnt, erosion_2p, fit)
#' @export
#'
ebLBH <- function(Xnonsample, f_q = ~1, data_2p, fit){

  ## argument parsing and checking ####
  Xnonsample <- data.frame(Xnonsample)
  if (length(all.vars(f_q))) {
    q <- Xnonsample[, all.vars(f_q)]
  } else {
    q <- rep(1, nrow(Xnonsample))
  }
  fs <- attributes(data_2p)
  area_oos <- Xnonsample[, all.vars(fs$f_area)]
  Xs1_oos <- model.matrix(fs$f_pos[-2], data = Xnonsample)
  if(length(fs$f_zero)==3) fs$f_zero <- fs$f_zero[-2]
  Xs0_oos <- model.matrix(fs$f_zero, data = Xnonsample)
  link <- attributes(fit)$link

  attach(data_2p, warn.conflicts = FALSE)
  attach(fit, warn.conflicts = FALSE)

  nis <- tapply(deltas, area, length)
  mis <- tapply(area_oos, area_oos, length)
  Gs <- Gmat(nis)
  Gns <- Gmat(mis)

  ## parameter estimates
  beta <- fixed$p1; alpha <- fixed$p0
  sig2le <- errorvar; sig2lu <- refvar1
  sig2lb <- refvar0; rho <- refcor
  ## conditional distribution center and scale
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

  cij <- exp(Xs1_oos%*%beta + Gns%*%(mu_u+sig2_u/2) + sig2le/2)

  ## Gauss-Hermite method
  ghwts <- do.call("rbind", lapply(1:D, function(i){
    wts <- statmod::gauss.quad.prob(20, dist = "normal", mu = b_m[i], sigma = sqrt(b_v[i]))
    as.data.frame(wts)
  }))
  bmat <- matrix(ghwts$nodes, nrow = D, ncol = 20, byrow = T)
  wmat <- matrix(ghwts$weights, nrow = D, ncol = 20, byrow = T)
  thres.gen <- function(b){
    mu_p <- Xs0%*%alpha + Gs%*%b
    p <- make.link(link)$linkinv(mu_p)
    pq <- as.vector(ifelse(deltas, p, 1-p))
    thres <- apply(Gs*pq, 2, function(x) prod(x[x>0]))
    thres
  }
  num.emp <- function(b){
    thres <- thres.gen(b)
    mu_p <- Xs0_oos%*%alpha + Gns%*%b
    pb <- make.link(link)$linkinv(mu_p)
    num <- Gns%*%(thres*(eta^b))*pb
    drop(num)
  }
  den.eb <- rowSums(apply(bmat, 2, thres.gen)*wmat)
  num.eb <- rowSums(apply(bmat, 2, num.emp)*(Gns%*%wmat))
  eb.peta <- num.eb/(Gns%*%den.eb)

  ybar.r <- drop(t(Gns) %*% (eb.peta*cij*q))/sum(q)
  ybar.s <- drop(t(Gs[deltas==1,]) %*% expo(lys))/nis
  eb <- drop(nis*ybar.s+sum(q)*ybar.r)/(nis+sum(q))

  return(data.frame(area = unique(area_oos), eb = eb))
}
