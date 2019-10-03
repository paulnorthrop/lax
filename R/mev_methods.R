# Methods for class "mev_gev", returned by mev::fit.gev()

#' @export
nobs.mev_gev <- function(object, ...) {
  return(length(object$xdat))
}

#' @export
coef.mev_gev <- function(object, ...) {
  return(object$estimate)
}

#' @export
vcov.mev_gev <- function(object, ...) {
  return(object$vcov)
}

#' @export
logLik.mev_gev <- function(object, ...) {
  val <- -object$nllh
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(coef(object))
  class(val) <- "logLik"
  return(val)
}
