# ============================ fExtremes::gpdFit ============================ #

# Methods for class fExtremes_gpd

#' @export
logLikVec.fExtremes_gpd <- function(object, pars = NULL, ...) {
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
  sigma <- pars["beta"]
  xi <- pars["xi"]
  # Calculate the loglikelihood contributions
  if (sigma <= 0) {
    val <- -Inf
  } else {
    val <- revdbayes::dgp(response_data, loc = object@fit$threshold,
                          scale = sigma, shape = xi, log = TRUE)
  }
  # Return the usual attributes for a "logLik" object
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(pars)
  class(val) <- "logLikVec"
  return(val)
}

#' @export
nobs.fExtremes_gpd <- function(object, ...) {
  return(sum(object@fit$data > object@fit$threshold))
}

#' @export
coef.fExtremes_gpd <- function(object, ...) {
  return(object@fit$par.ests)
}

#' @export
vcov.fExtremes_gpd <- function(object, ...) {
  vc <- object@fit$varcov
  par_names <- names(coef(object))
  dimnames(vc) <- list(par_names, par_names)
  return(vc)
}

#' @export
logLik.fExtremes_gpd <- function(object, ...) {
  return(logLik(logLikVec(object)))
}

# See fExtremes_methods.R for nobs, coef, vcov, logLik methods for
# class "fGPDFIT"
