# nobs, coef, vcov and logLik methods for class "gev", produced by evir::gev()

#' @export
nobs.gev <- function(object, ...) {
  # If this is an "evd" object then use the "evd" method
  if (inherits(object, "evd")) {
    # c("gev", "uvevd", "evd") becomes c("evd", "uvevd", "gev")
    class(object) <- rev(class(object))
    return(nobs(object, ...))
  }
  return(object$n)
}

#' @export
coef.gev <- function(object, ...) {
  # If this is an "evd" object then use the "evd" method
  if (inherits(object, "evd")) {
    # c("gev", "uvevd", "evd") becomes c("evd", "uvevd", "gev")
    class(object) <- rev(class(object))
    return(coef(object, ...))
  }
  return(object$par.ests)
}

#' @export
vcov.gev <- function(object, ...) {
  # If this is an "evd" object then use the "evd" method
  if (inherits(object, "evd")) {
    # c("gev", "uvevd", "evd") becomes c("evd", "uvevd", "gev")
    class(object) <- rev(class(object))
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
    # c("gev", "uvevd", "evd") becomes c("evd", "uvevd", "gev")
    class(object) <- rev(class(object))
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

# nobs, coef, vcov and logLik methods for class "potd", produced by evir::pot()

#' @export
nobs.potd <- function(object, ...) {
  return(object$n)
}

#' @export
coef.potd <- function(object, complete = FALSE, ...) {
  if (complete) {
    val <- object$par.ests
  } else {
    val <- object$par.ests[c("xi", "sigma", "mu")]
  }
  return(val)
}

#' @export
vcov.potd <- function(object, complete = FALSE, ...) {
  vc <- object$varcov
  par_names <- names(coef(object))
  if (complete) {
    vc <- rbind(cbind(vc, NA), NA)
    par_names <- c(names(coef(object)), "beta")
  }
  dimnames(vc) <- list(par_names, par_names)
  return(vc)
}

#' @export
logLik.potd <- function(object, ...) {
  val <- -object$nllh.final
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(coef(object))
  class(val) <- "logLik"
  return(val)
}
