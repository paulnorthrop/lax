# ================================ evir::gev ================================ #

# Methods for class evir_gev

#' @export
logLikVec.evir_gev <- function(object, pars = NULL, ...) {
  if (!missing(...)) {
    warning("extra arguments discarded")
  }
  # If the parameter estimates have not been provided in pars then extract
  # them from the fitted object
  if (is.null(pars)) {
    pars <- coef(object)
  }
  n_pars <- length(pars)
  response_data <- object$data
  mu <- pars["mu"]
  sigma <- pars["sigma"]
  xi <- pars["xi"]
  # Calculate the loglikelihood contributions
  if (sigma <= 0) {
    val <- -Inf
  } else {
    val <- revdbayes::dgev(object$data, loc = mu, scale = sigma, shape = xi,
                           log = TRUE)
  }
  # Return the usual attributes for a "logLik" object
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(pars)
  class(val) <- "logLikVec"
  return(val)
}

#' @export
nobs.evir_gev <- function(object, ...) {
  return(object$n)
}

#' @export
coef.evir_gev <- function(object, ...) {
  return(object$par.ests)
}

#' @export
vcov.evir_gev <- function(object, complete = FALSE, ...) {
  vc <- object$varcov
  par_names <- names(coef(object))
  dimnames(vc) <- list(par_names, par_names)
  return(vc)
}

#' @export
logLik.evir_gev <- function(object, ...) {
  return(logLik(logLikVec(object)))
}

# See evir_methods.R for nobs, coef, vcov, logLik methods for class "gev"
