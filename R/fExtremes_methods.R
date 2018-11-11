# nobs, coef, vcov and logLik methods for class "fGEVFIT", produced by
# fExtremes::gevFit()

#' @export
nobs.fGEVFIT <- function(object, ...) {
  return(object@fit$n)
}

#' @export
coef.fGEVFIT <- function(object, ...) {
  return(object@fit$par.ests)
}

#' @export
vcov.fGEVFIT <- function(object, ...) {
  return(object@fit$varcov)
}

#' @export
logLik.fGEVFIT <- function(object, ...) {
  val <- -object@fit$nllh.final
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(coef(object))
  class(val) <- "logLik"
  return(val)
}

