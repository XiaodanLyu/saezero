#' Converts a data frame to a list made for fitting LBH model
#'
#' The default method transforms the data into a format that fits the framework of the unit level model of Lyu, Berg and Hofmann.
#'
#' @param f_pos an object of class \code{\link[stats]{formula}}:
#'   a symbolic description of the fixed effect model to be fitted to the positive part.
#' @param f_zero an object of class \code{\link[stats]{formula}}:
#'   a symbolic description of the fixed effect model to be fitted to the binary part.
#'   Default value is to using the same formula as the positive part (\code{f_pos}).
#' @param f_area an object of class \code{\link[stats]{formula}}:
#'   a symbolic description of the area code to be fitted to both the positive part and the negative part.
#' @param data sample data containing the variables named in \code{f_pos}, \code{f_zero} and \code{f_area}.
#'   Any sampling units containing missing entries in the model data frame are removed.
#' @param transform type of transformation for the positive responses to be chosen between the
#'   "BoxCox" and "power" families. Default value is "Boxcox".
#' @param lambda value for the parameter of the family of transformations specified in \code{transform}.
#'   Default value is 0, which gives the log transformation for the two possible families.
#'
#' @return An object of the class "2pdata" which is a list with the following components:
#' \itemize{
#'   \item \code{lys}: transformed response vector for the positive part
#'   \item \code{Xs1}: model matrix for the positive part
#'   \item \code{deltas}: binary response vector for the binary part
#'   \item \code{Xs0}: model matrix for the binary part
#'   \item \code{area}: vector with the area code
#' }
#'
#' @details The response variable in the formula \code{f_zero} is ignored.
#'   \code{I(y>0)} will be used for the binary part where \code{y} is the
#'   response variable in the formula \code{f_pos}. \cr
#'   Although this function
#'   allows a general transformation family for the positive
#'   part, the Empirical Bayes predictor of area means developed in the paper of
#'   Lyu, Berg and Hofmann (\code{\link{ebLBH}}) assumes a log transformation
#'   is used. A general transformation family is
#'   provided in this function so as to facilitate a profile likelihood analysis
#'   together with the function \code{\link{mleLBH}} to check the assumption of a log transformation
#'   for the positive part, i.e., whether lambda is significantly different from 0.
#'   See \href{https://github.com/XiaodanLyu/zero-sae-bj-2020/blob/master/case_study/case_study.pdf}{here}
#'   for an example of profiling lambda with respect to the BoxCox transformation family.
#' @seealso \code{\link[sae]{bxcx}}
#'
#' @examples
#'   erosion_2p <- as.2pdata(f_pos = RUSLE2~logR+logK+logS,
#'                           f_zero = ~logR+logS+crop2+crop3,
#'                           f_area = ~cty, data = erosion)
#' @export
as.2pdata <- function(f_pos, f_zero = f_pos, f_area, data, transform = "BoxCox", lambda = 0){

  ## argument parsing and checking

  if (!missing(data)){
    data <- as.data.frame(data)
    fposdata <- model.frame(f_pos, na.action = na.pass, data)
    fzerodata <- model.frame(f_zero, na.action = na.pass, data)
    area <- data[, all.vars(f_area)]
  } else{
    fposdata <- model.frame(f_pos, na.action = na.pass)
    fzerodata <- model.frame(f_zero, na.action = na.pass)
    area <- model.frame(f_area, na.action = na.pass)
  }

  if (!is.vector(area)) area <- as.vector(area)

  if (nrow(fposdata)!=length(area))
    stop("Arguments f_pos [nrows=",nrow(fposdata),"] and area [nrows=",length(area),"] must be the same length.\n")
  if (nrow(fzerodata)!=length(area))
    stop("Arguments f_zero [nrows=",nrow(fzerodata),"] and area [nrows=",length(area),"] must be the same length.\n")

  # Delete rows with NA values
  rowNA1 <- which(apply(fposdata, 1, function(x) any(is.na(x))))
  rowNA0 <- which(apply(fzerodata, 1, function(x) any(is.na(x))))
  omitted <- unique(c(rowNA1, rowNA0, which(is.na(area))))
  if (length(omitted)) {
    area <- area[-omitted]
    fposdata <- fposdata[-omitted,]
    fzerodata <- fzerodata[-omitted,]
  }

  if(any(fposdata[,1]<0)) stop("Responses must be nonnegative.")

  deltas <- ifelse(fposdata[,1]>0, 1, 0)
  Xs0  <- model.matrix(f_zero, fzerodata)

  fposdata <- fposdata[which(deltas==1),]
  # lys <- log(fposdata[,1])
  lys <- sae::bxcx(fposdata[,1], lambda = lambda, type = transform)
  Xs1  <- model.matrix(f_pos, fposdata)

  result <- list(lys = lys, Xs1 = Xs1, deltas = deltas, Xs0 = Xs0, area = area)
  attr(result, "f_pos") <- f_pos
  if(length(f_zero)==3) f_zero <- f_zero[-2]
  attr(result, "f_zero") <- f_zero
  attr(result, "f_area") <- f_area
  attr(result, "transform") <- transform
  attr(result, "lambda") <- lambda

  return(result)
}
