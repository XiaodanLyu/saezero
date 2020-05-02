#' MLE model parameter estimator
#'
#' Fits by ML method the unit level model of Lyu, Berg and Hofmann.
#' The specified link function is used in the binary part.
#'
#' @param data_2p a two-part data object returned by \code{\link{as.2pdata}}.
#' @param link a specification for the link function used to model the binary part.
#'   The accepted link functions are \code{logit}, \code{probit}, \code{cauchit}, \code{log} and \code{cloglog}.
#'   Default value is "logit".
#' @return The function returns a list with the following objects:
#' \itemize{
#'   \item \code{fixed}: list with the estimated values of the fixed regression coefficient
#'     in the positive part (\code{p1}) and in the binary part (\code{p0}).
#'   \item \code{random}: data frame with the predicted random effects in the positive part (\code{p1})
#'     and the binary part (\code{p0}).
#'   \item \code{errorvar}: estimated model error variance in the positive part.
#'   \item \code{refvar1}: estimated random effects variance in the positive part.
#'   \item \code{refvar0}: estimated random effects variance in the binary part.
#'   \item \code{refcor}: estimated correlation coefficient of the random effects between the two parts.
#'   \item \code{loglik}: log-likelihood accomodating a general transformation family and \code{lambda}
#'     in the positive part as specified in \code{\link{as.2pdata}}.
#'   \item \code{residuals}: the marginal (\code{mar}) and conditional (\code{con}) residuals
#'     from the model fit in the positive part (\code{p1}, cases with nonpositive response are NA.).
#'   \item \code{vcov}: list of the estimated variance-covariance matrices of the linear coefficients
#'     in the positive part (\code{p1}) and in the binary part (\code{p0}).
#'   \item \code{cirefcor}: function that has one argument -- confidence level (\code{alpha}) and
#'     returns \eqn{1-\alpha} confidence interval of the correlation coefficient.
#'   \item \code{fit0}: model parameter estimator under independence assumption.
#' }
#'
#' @seealso \code{\link{as.2pdata}}
#' @examples
#'   erosion_2p <- as.2pdata(f_pos = RUSLE2~logR+logK+logS,
#'                           f_zero = ~logR+logS+crop2+crop3,
#'                           f_area = ~cty, data = erosion)
#'   fit <- mleLBH(erosion_2p)
#' @import stats
#' @export
mleLBH <- function(data_2p, link = "logit"){

  result <- list(fixed = NA, random = NA, errorvar = NA, refvar1 = NA,
                 refvar0 = NA, refcor = NA, loglik = NA, residuals = NA,
                 vcov = NA, cirefcor = NA, fit0 = NA)
  attr(result, "link") <- link

  ## argument parsing and checking ####
  if (!(link %in% c("logit", "probit", "cauchit", "log", "cloglog")))
    stop("Argument link must be \"logit\", \"probit\", \"cauchit\", \"log\" or \"cloglog\".")

  attach(data_2p, warn.conflicts = FALSE)

  ## model fit ####
  # 1. fit independence model to get appropriate initial values
  area_pos <- area[which(deltas==1)]
  fit_pos <- lme4::lmer(lys~Xs1-1+(1|area_pos))
  b_ind <- matrix(lme4::fixef(fit_pos), ncol=1)
  sig2u_ind <- unname(unlist(lme4::VarCorr(fit_pos)))
  sige_ind <- sigma(fit_pos)
  uihat_ind <- unlist(lme4::ranef(fit_pos))

  suppressMessages(suppressWarnings(
    fit_zero <- lme4::glmer(deltas~Xs0-1+(1|area), family = binomial(link), nAGQ = 10)
  ))

  a_ind <- matrix(lme4::fixef(fit_zero), ncol = 1)
  sig2b_ind <- unlist(lme4::VarCorr(fit_zero))
  bihat_ind <- unlist(lme4::ranef(fit_zero))

  nti <- tapply(deltas, area, sum)
  suppressWarnings(
    rho_ind <- cor(uihat_ind, bihat_ind[nti>0])
  )
  if(is.na(rho_ind)) rho_ind <- 0
  # 2. maximize log-likelihood over parameter space
  theta0 <- c(b_ind, sqrt(sig2u_ind), sige_ind,
              a_ind, sqrt(max(0.0098, sig2b_ind)), tan(pi/2*rho_ind))

  result$fit0 <- list(fixed = list(p1 = b_ind, p0 = a_ind),
                      errorvar = sige_ind^2, refvar1 = sig2u_ind,
                      refvar0 = sig2b_ind, refcor = 0)

  res.theta <- optim(theta0, fn = loglik.neg, gr = loglik.neg.grad,
                     method = "BFGS", hessian = TRUE,
                     lys, Xs1, deltas, Xs0, area, link)

  if(res.theta$convergence!=0){
    warning("convergence error code: ", res.theta$convergence, ", please see ?optim")
  }

  thetahat <- res.theta$par
  p1 <- ncol(Xs1)
  p2 <- ncol(Xs0)
  result$fixed <- list(p1 = thetahat[1:p1], p0 = thetahat[(p1+3):(p1+p2+2)])
  result$errorvar <- thetahat[p1+2]^2
  result$refvar1 <- thetahat[p1+1]^2
  result$refvar0 <- thetahat[p1+p2+3]^2
  result$refcor <- atan(thetahat[p1+p2+4])*2/pi
  ## loglik term of lambda
  type <- attributes(data_2p)$transform
  lambda <- attributes(data_2p)$lambda
  ys <- sae::bxcx(lys, lambda, InverseQ = TRUE, type = type)
  lllambda <- (lambda-1)*sum(log(ys)) + ifelse(type=="power", sum(nti)*log(lambda), 0)
  ## full log-likelihood
  result$loglik <- -1*res.theta$value + lllambda
  vcovmat <- solve(res.theta$hessian)
  result$vcov <- list(p1 = vcovmat[1:p1, 1:p1],
                      p0 = vcovmat[(p1+3):(p1+p2+2), (p1+3):(p1+p2+2)])
  result$cirefcor <- function(alpha = 0.05){
    ci_t <- thetahat[p1+p2+4]+c(-1,1)*qnorm(1-alpha/2)*sqrt(vcovmat[p1+p2+4, p1+p2+4])
    atan(ci_t)*2/pi
  }

  ## EB predictor of random effects
  result$random <- ranefLBH(thetahat, lys, Xs1, deltas, Xs0, area, link)
  ## fitted values
  lyshat.mar <- Xs1 %*% (result$fixed$p1)
  lyshat.con <- lyshat.mar + Gmat(nti) %*% (result$random$p1[nti>0])
  result$residuals <- data.frame(mar = rep(NA, length(area)),
                                 con = rep(NA, length(area)))
  result$residuals[deltas==1, "mar"] <- lys - lyshat.mar
  result$residuals[deltas==1, "con"] <- lys - lyshat.con

  detach(data_2p)
  return(result)
}

ranefLBH <- function(theta, lys, Xs1, deltas, Xs0, area, link = "logit"){

  p1 <- ncol(Xs1)
  p2 <- ncol(Xs0)
  beta <- theta[1:p1]
  sig2lu <- theta[p1+1]^2
  sig2le <- theta[p1+2]^2
  alpha <- theta[(p1+3):(p1+p2+2)]
  sig2lb <- theta[p1+p2+3]^2
  rho <- atan(theta[p1+p2+4])*2/pi

  nis <- tapply(deltas, area, length)
  D <- length(nis)
  Gs <- Gmat(nis)
  nti <- tapply(deltas, area, sum)
  rt <- rep(0, length(area))
  rt[deltas==1] <- lys - Xs1%*%beta

  gamma <- (1-rho^2)*sig2lu/((1-rho^2)*sig2lu+sig2le/nti)
  rbar <- drop(t(Gs)%*%rt)/nti
  rbar[nti==0] <- 0
  b_m <- rho*sqrt(sig2lu*sig2lb)*rbar/(sig2lu+sig2le/nti)
  b_v <- (1-rho^2)/(1-(1-gamma)*rho^2)*sig2lb

  ## Gauss-Hermite method
  ghwts <- do.call("rbind", lapply(1:D, function(i){
    wts <- statmod::gauss.quad.prob(20, dist = "normal", mu = b_m[i], sigma = sqrt(b_v[i]))
    as.data.frame(wts)
  }))
  bmat <- matrix(ghwts$nodes, nrow = D, ncol = 20, byrow = T)
  wmat <- matrix(ghwts$weights, nrow = D, ncol = 20, byrow = T)
  numint <- function(g){
    arr <- lapply(1:20, function(r){
      mu_p <- Xs0%*%alpha + Gs%*%bmat[,r]
      p <- make.link(link)$linkinv(mu_p)
      pq <- as.vector(ifelse(deltas, p, 1-p))
      thres <- apply(Gs*pq, 2, function(x) prod(x[x>0]))
      return(as.matrix(thres*g(bmat[,r])*wmat[,r]))
    })
    return(apply(simplify2array(arr), 1:2, sum))
  }

  bhat <- numint(function(b) b)/numint(function(b) 1)

  mu_u <- gamma*rbar
  mu_u[nti==0] <- 0
  uhat <- mu_u + (1-gamma)*rho*sqrt(sig2lu/sig2lb)*bhat[,1]

  return(data.frame(p1 = uhat, p0 = bhat))
}

loglik.neg <- function(theta, lys, Xs1, deltas, Xs0, area, link = "logit"){

  p1 <- ncol(Xs1)
  p2 <- ncol(Xs0)
  beta <- theta[1:p1]
  sig2lu <- theta[p1+1]^2
  sig2le <- theta[p1+2]^2
  alpha <- theta[(p1+3):(p1+p2+2)]
  sig2lb <- theta[p1+p2+3]^2
  rho <- atan(theta[p1+p2+4])*2/pi

  nis <- tapply(deltas, area, length)
  D <- length(nis)
  Gs <- Gmat(nis)
  nti <- tapply(deltas, area, sum)
  rt <- rep(0, length(area))
  rt[deltas==1] <- lys - Xs1%*%beta

  gamma <- (1-rho^2)*sig2lu/((1-rho^2)*sig2lu+sig2le/nti)
  l1 <- log(1-gamma) + log(1-rho^2) - log(1-rho^2*(1-gamma)) - nti*log(sig2le)
  rbar <- drop(t(Gs)%*%rt)/nti
  rbar[nti==0] <- 0
  l2 <- gamma*rbar^2/(sig2le/nti)
  b_m <- rho*sqrt(sig2lu*sig2lb)*rbar/(sig2lu+sig2le/nti)
  b_v <- (1-rho^2)/(1-(1-gamma)*rho^2)*sig2lb
  l3 <- b_m^2/b_v
  l4 <- t(Gs)%*%(rt^2)/sig2le

  ## Gauss-Hermite method
  ghwts <- do.call("rbind", lapply(1:D, function(i){
    wts <- statmod::gauss.quad.prob(20, dist = "normal", mu = b_m[i], sigma = sqrt(b_v[i]))
    as.data.frame(wts)
  }))
  bmat <- matrix(ghwts$nodes, nrow = D, ncol = 20, byrow = T)
  wmat <- matrix(ghwts$weights, nrow = D, ncol = 20, byrow = T)
  numint <- function(g){
    arr <- lapply(1:20, function(r){
      mu_p <- Xs0%*%alpha + Gs%*%bmat[,r]
      p <- make.link(link)$linkinv(mu_p)
      pq <- as.vector(ifelse(deltas, p, 1-p))
      thres <- apply(Gs*pq, 2, function(x) prod(x[x>0]))
      return(as.matrix(thres*g(bmat[,r])*wmat[,r]))
    })
    return(apply(simplify2array(arr), 1:2, sum))
  }

  l5 <- log(numint(function(b) 1))
  loglik <- sum(1/2*(l1+l2+l3-l4[,1])+l5[,1])

  return(-1*loglik)
}

loglik.neg.grad <- function(theta, lys, Xs1, deltas, Xs0, area, link = "logit"){

  p1 <- ncol(Xs1)
  p2 <- ncol(Xs0)
  beta <- theta[1:p1]
  sig2lu <- theta[p1+1]^2
  sig2le <- theta[p1+2]^2
  alpha <- theta[(p1+3):(p1+p2+2)]
  sig2lb <- theta[p1+p2+3]^2
  t <- theta[p1+p2+4]
  rho <- atan(t)*2/pi

  nis <- tapply(deltas, area, length)
  D <- length(nis)
  Gs <- Gmat(nis)
  nti <- tapply(deltas, area, sum)
  rt <- rep(0, length(area))
  rt[deltas==1] <- lys - Xs1%*%beta
  rbar <- drop(t(Gs)%*%rt)/nti
  rbar[nti==0] <- 0
  xt <- matrix(0, length(area), p1)
  xt[deltas==1,] <- Xs1
  xbar <- (t(Gs)%*%xt)/as.vector(nti)
  xbar[nti==0,] <- 0

  gamma <- (1-rho^2)*sig2lu/((1-rho^2)*sig2lu+sig2le/nti)
  b_m <- rho*sqrt(sig2lu*sig2lb)*rbar/(sig2lu+sig2le/nti)
  b_v <- (1-rho^2)/(1-(1-gamma)*rho^2)*sig2lb

  ## Gauss-Hermite method
  ghwts <- do.call("rbind", lapply(1:D, function(i){
    wts <- statmod::gauss.quad.prob(20, dist = "normal", mu = b_m[i], sigma = sqrt(b_v[i]))
    as.data.frame(wts)
  }))
  bmat <- matrix(ghwts$nodes, nrow = D, ncol = 20, byrow = T)
  wmat <- matrix(ghwts$weights, nrow = D, ncol = 20, byrow = T)
  numint <- function(g){
    arr <- lapply(1:20, function(r){
      mu_p <- Xs0%*%alpha + Gs%*%bmat[,r]
      p <- make.link(link)$linkinv(mu_p)
      pq <- as.vector(ifelse(deltas, p, 1-p))
      thres <- apply(Gs*pq, 2, function(x) prod(x[x>0]))
      return(as.matrix(thres*g(bmat[,r])*wmat[,r]))
    })
    return(apply(simplify2array(arr), 1:2, sum))
  }

  den <- numint(function(b) 1)

  ## partial beta
  mb <- b_m/rbar
  mb[nti==0] <- 0
  tmp <- -gamma*rbar/(sig2le/nti) - b_m/b_v*mb + drop(numint(function(b) -(b-b_m)/sqrt(b_v))/den)*mb/sqrt(b_v)
  pbeta <- drop(t(xt)%*%rt/sig2le) + colSums(as.vector(tmp)*xbar)

  ## partial alpha
  pppa <- function(b){
    mu_p <- Xs0%*%alpha + Gs %*% b
    p <- make.link(link)$linkinv(mu_p)
    qp <- as.vector(ifelse(deltas, 1-p, -p))
    t(Gs)%*%(Xs0*qp)
  }
  palpha <- colSums(numint(pppa)/drop(den))

  pgpu <- 2*gamma*(1-gamma)/sqrt(sig2lu)
  pgpe <- -2*gamma*(1-gamma)/sqrt(sig2le)
  pvpu <- -rho^2/(1-(1-gamma)*rho^2)*b_v*pgpu
  pvpe <- -rho^2/(1-(1-gamma)*rho^2)*b_v*pgpe
  pvpb <- 2*b_v/sqrt(sig2lb)
  pmpu <- (sig2le/nti-sig2lu)/(sig2le/nti+sig2lu)*b_m/sqrt(sig2lu)
  pmpu[nti==0] <- 0
  pmpe <- -2*b_m*(sqrt(sig2le)/nti)/(sig2le/nti+sig2lu)
  pmpe[nti==0] <- 0
  pmpb <- b_m/sqrt(sig2lb)

  ## partial sigu
  psigu <- sum(1/2*(-pgpu/(1-gamma)+rbar^2*pgpu/(sig2le/nti)+
                      2*pmpu*b_m/b_v-(b_m/b_v)^2*pvpu))+
    colSums(numint(function(b) (b-b_m)/b_v*pmpu+(b-b_m)^2/(2*b_v^2)*pvpu)/den)
  ## partial sige
  psige <- sum(1/2*(-pgpe/(1-gamma)-2*nti/sqrt(sig2le)+
                      rbar^2/(sig2le/nti)*pgpe-2*gamma*rbar^2/(sig2le^(3/2)/nti)+
                      2*b_m/b_v*pmpe-(b_m/b_v)^2*pvpe+
                      2/(sig2le^(3/2))*drop(t(Gs)%*%(rt^2))))+
    colSums(numint(function(b) (b-b_m)/b_v*pmpe+(b-b_m)^2/(2*b_v^2)*pvpe)/den)
  ## partial sigb
  psigb <- sum(1/2*(2*b_m/b_v*pmpb-(b_m/b_v)^2*pvpb-2/sqrt(sig2lb)))+
    colSums(numint(function(b) (b-b_m)/b_v*pmpb+(b-b_m)^2/(2*b_v^2)*pvpb)/den)

  ## partial rho
  pgpr <- -2*rho/(1-rho^2)*gamma*(1-gamma)
  pvpr <- -2*rho*gamma*sig2lb
  pmpr <- b_m/rho
  prpt <- (2/pi)/(1+t^2)
  pr <- sum(1/2*(-pgpr/(1-gamma)+rbar^2*pgpr/(sig2le/nti)+
                   2*pmpr*b_m/b_v-(b_m/b_v)^2*pvpr))+
    colSums(numint(function(b) (b-b_m)/b_v*pmpr+(b-b_m)^2/(2*b_v^2)*pvpr)/den)
  pt <- pr*prpt

  return(-c(pbeta, psigu, psige, palpha, psigb, pt))
}
