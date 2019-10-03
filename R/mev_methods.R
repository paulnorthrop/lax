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

# Methods for class "mev_gpd"

#' @export
nobs.mev_gpd <- function(object, ...) {
  return(object$nat)
}

#' @export
coef.mev_gpd <- function(object, ...) {
  return(object$estimate)
}

#' @export
vcov.mev_gpd <- function(object, ...) {
  vc <- object$vcov
  dimnames(vc) <- list(names(coef(object)), names(coef(object)))
  return(vc)
}

#' @export
logLik.mev_gpd <- function(object, ...) {
  val <- -object$nllh
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(coef(object))
  class(val) <- "logLik"
  return(val)
}

# Methods for class "mev_pp", returned by mev::fit.pp()

#' @export
nobs.mev_pp <- function(object, ...) {
  return(object$nat / object$pat)
}

#' @export
coef.mev_pp <- function(object, ...) {
  return(object$estimate)
}

#' @export
vcov.mev_pp <- function(object, ...) {
  return(object$vc)
}

#' @export
logLik.mev_pp <- function(object, ...) {
  val <- -object$nllh
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(coef(object))
  class(val) <- "logLik"
  return(val)
}

