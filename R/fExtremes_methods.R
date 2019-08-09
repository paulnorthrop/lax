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

# nobs, coef, vcov and logLik methods for class "fGPDFIT", produced by
# fExtremes::gpdFit()

#' @export
nobs.fGPDFIT <- function(object, ...) {
  return(sum(object@fit$data > object@fit$threshold))
}

#' @export
coef.fGPDFIT <- function(object, ...) {
  return(object@fit$par.ests)
}

#' @export
vcov.fGPDFIT <- function(object, ...) {
  vc <- object@fit$varcov
  par_names <- names(coef(object))
  dimnames(vc) <- list(par_names, par_names)
  return(vc)
}

#' @export
logLik.fGPDFIT <- function(object, ...) {
  # Note: value is the *negated* loglikelihood at the MLE
  val <- -object@fit$fit$value
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(coef(object))
  class(val) <- "logLik"
  return(val)
}

