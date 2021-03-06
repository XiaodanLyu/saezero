% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mleLBH.r
\name{mleLBH}
\alias{mleLBH}
\title{MLE model parameter estimator}
\usage{
mleLBH(data_2p, link = "logit")
}
\arguments{
\item{data_2p}{a two-part data object returned by \code{\link{as.2pdata}}.}

\item{link}{a specification for the link function used to model the binary part.
The accepted link functions are \code{logit}, \code{probit}, \code{cauchit}, \code{log}, \code{cloglog}
and an object of the class \code{link-glm}, see the documentation of \code{\link[stats]{make.link}}.
Default value is "logit".}
}
\value{
The function returns a list with the following objects:
\itemize{
  \item \code{fixed}: list with the estimated values of the fixed regression coefficient
    in the positive part (\code{p1}) and in the binary part (\code{p0}).
  \item \code{random}: data frame with the predicted random effects in the positive part (\code{p1})
    and the binary part (\code{p0}).
  \item \code{errorvar}: estimated model error variance in the positive part.
  \item \code{refvar1}: estimated random effects variance in the positive part.
  \item \code{refvar0}: estimated random effects variance in the binary part.
  \item \code{refcor}: estimated correlation coefficient of the random effects between the two parts.
  \item \code{loglik}: log-likelihood accomodating a general transformation family and \code{lambda}
    in the positive part as specified in \code{\link{as.2pdata}}.
  \item \code{residuals}: the marginal (\code{mar}) and conditional (\code{con}) residuals
    from the model fit in the positive part (\code{p1}, cases with nonpositive response are NA.).
  \item \code{fit0}: model parameter estimator under independence assumption.
}
}
\description{
Fits by ML method the unit level model of Lyu, Berg and Hofmann.
The specified link function is used in the binary part.
}
\examples{
  erosion_2p <- as.2pdata(f_pos = RUSLE2~logR+logK+logS,
                          f_zero = ~logR+logS+crop2+crop3,
                          f_area = ~cty, data = erosion)
  fit <- mleLBH(erosion_2p)
}
\seealso{
\code{\link{as.2pdata}}
}
