# nobs, coef, vcov and logLik methods for class "gev", produced by evir::gev()

#' @export
nobs.gev <- function(object, ...) {
  # If this is an "evd" object then use the "evd" method
  if (inherits(object, "evd")) {
    class(object) <- "evd"
    return(nobs(object, ...))
  }
  return(object$n)
}

#' @export
coef.gev <- function(object, ...) {
  # If this is an "evd" object then use the "evd" method
  if (inherits(object, "evd")) {
    class(object) <- "evd"
    return(coef(object, ...))
  }
  return(object$par.ests)
}

#' @export
vcov.gev <- function(object, ...) {
  # If this is an "evd" object then use the "evd" method
  if (inherits(object, "evd")) {
    class(object) <- "evd"
    return(vcov(object, ...))
  }
  vc <- object$varcov
  par_names <- names(coef(object))
  dimnames(vc) <- list(par_names, par_names)
  return(vc)
}

#' @export
logLik.gev <- function(object, ...) {
  # If this is an "evd" object then use the "evd" method
  if (inherits(object, "evd")) {
    class(object) <- "evd"
    return(logLik(object, ...))
  }
  val <- -object$nllh.final
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(coef(object))
  class(val) <- "logLik"
  return(val)
}

# nobs, coef, vcov and logLik methods for class "gpd", produced by evir::gpd()

#' @export
nobs.gpd <- function(object, ...) {
  return(object$n.exceed)
}

#' @export
coef.gpd <- function(object, ...) {
  return(object$par.ests)
}

#' @export
vcov.gpd <- function(object, ...) {
  vc <- object$varcov
  par_names <- names(coef(object))
  dimnames(vc) <- list(par_names, par_names)
  return(vc)
}

#' @export
logLik.gpd <- function(object, ...) {
  val <- -object$nllh.final
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(coef(object))
  class(val) <- "logLik"
  return(val)
}
