# Methods for class "gev.fit", returned by ismev::gev.fit()
# ismev_gev_names() is in ismev_gev_fit.R

#' @export
nobs.gev.fit <- function(object, ...) {
  return(length(object$data))
}

#' @export
coef.gev.fit <- function(object, ...) {
  val <- object$mle
  names(val) <- ismev_gev_names(object)
  return(val)
}

#' @export
vcov.gev.fit <- function(object, ...) {
  vc <- object$cov
  par_names <- ismev_gev_names(object)
  dimnames(vc) <- list(par_names, par_names)
  return(vc)
}

#' @export
logLik.gev.fit <- function(object, ...) {
  val <- -object$nllh
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(coef(object))
  class(val) <- "logLik"
  return(val)
}

# Methods for class "gpd.fit"
# ismev_gpd_names() is in ismev_gpd_fit.R

#' @export
nobs.gpd.fit <- function(object, ...) {
  return(object$nexc)
}

#' @export
coef.gpd.fit <- function(object, ...) {
  val <- object$mle
  names(val) <- ismev_gpd_names(object)
  return(val)
}

#' @export
vcov.gpd.fit <- function(object, ...) {
  vc <- object$cov
  dimnames(vc) <- list(ismev_gpd_names(object), ismev_gpd_names(object))
  return(vc)
}

#' @export
logLik.gpd.fit <- function(object, ...) {
  val <- -object$nllh
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(coef(object))
  class(val) <- "logLik"
  return(val)
}

# Methods for class "pp.fit"
# ismev_gev_names() is in ismev_gev_fit.R

#' @export
nobs.pp.fit <- function(object, ...) {
  return(nrow(object$vals))
}

#' @export
coef.pp.fit <- function(object, ...) {
  val <- object$mle
  names(val) <- ismev_gev_names(object)
  return(val)
}

#' @export
vcov.pp.fit <- function(object, ...) {
  vc <- object$cov
  dimnames(vc) <- list(ismev_gev_names(object), ismev_gev_names(object))
  return(vc)
}

#' @export
logLik.pp.fit <- function(object, ...) {
  val <- -object$nllh
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(coef(object))
  class(val) <- "logLik"
  return(val)
}

# Methods for class "pp.fit"
# ismev_gev_names() is in ismev_gev_fit.R

#' @export
nobs.rlarg.fit <- function(object, ...) {
  return(nrow(object$vals))
}

#' @export
coef.rlarg.fit <- function(object, ...) {
  val <- object$mle
  names(val) <- ismev_gev_names(object)
  return(val)
}

#' @export
vcov.rlarg.fit <- function(object, ...) {
  vc <- object$cov
  dimnames(vc) <- list(ismev_gev_names(object), ismev_gev_names(object))
  return(vc)
}

#' @export
logLik.rlarg.fit <- function(object, ...) {
  val <- -object$nllh
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(coef(object))
  class(val) <- "logLik"
  return(val)
}
