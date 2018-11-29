# ============================ fExtremes::gevFit ============================ #

# Methods for class fExtremes_gev

#' @export
logLikVec.fExtremes_gev <- function(object, pars = NULL, ...) {
  if (!missing(...)) {
    warning("extra arguments discarded")
  }
  # If the parameter estimates have not been provided in pars then extract
  # them from the fitted object
  if (is.null(pars)) {
    pars <- coef(object)
  }
  n_pars <- length(pars)
  response_data <- object@fit$data
  mu <- pars["mu"]
  sigma <- pars["beta"]
  xi <- pars["xi"]
  # Calculate the loglikelihood contributions
  if (sigma <= 0) {
    val <- -Inf
  } else {
    val <- revdbayes::dgev(response_data, loc = mu, scale = sigma, shape = xi,
                           log = TRUE)
  }
  # Return the usual attributes for a "logLik" object
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(pars)
  class(val) <- "logLikVec"
  return(val)
}

#' @export
nobs.fExtremes_gev <- function(object, ...) {
  return(object@fit$n)
}

#' @export
coef.fExtremes_gev <- function(object, ...) {
  return(object@fit$par.ests)
}

#' @export
vcov.fExtremes_gev <- function(object, ...) {
  return(object@fit$varcov)
}

#' @export
logLik.fExtremes_gev <- function(object, ...) {
  return(logLik(logLikVec(object)))
}

# See fExtremes_methods.R for nobs, coef, vcov, logLik methods for
# class "fGEVFIT"
