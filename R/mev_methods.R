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

# Methods for class "mev_egp", returned by mev::fit.egp()

#' @export
nobs.mev_egp <- function(object, ...) {
  return(object$nat)
}

#' @export
coef.mev_egp <- function(object, ...) {
  return(object$estimate)
}

#' @export
vcov.mev_egp <- function(object, ...) {
  return(object$vcov)
}

#' @export
logLik.mev_egp <- function(object, ...) {
  val <- -object$deviance / 2
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(coef(object))
  class(val) <- "logLik"
  return(val)
}

# Methods for class "mev_rlarg", returned by mev::fit.rlarg()

#' @export
nobs.mev_rlarg <- function(object, ...) {
  return(nrow(object$xdat))
}

#' @export
coef.mev_rlarg <- function(object, ...) {
  return(object$estimate)
}

#' @export
vcov.mev_rlarg <- function(object, ...) {
  vc <- object$vcov
  dimnames(vc) <- list(names(coef(object)), names(coef(object)))
  return(vc)
}

#' @export
logLik.mev_rlarg <- function(object, ...) {
  val <- -object$nllh
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(coef(object))
  class(val) <- "logLik"
  return(val)
}
