# ================================= evir::gpd =============================== #

# Methods for class evir_gpd

#' @export
logLikVec.evir_gpd <- function(object, pars = NULL, ...) {
  if (!missing(...)) {
    warning("extra arguments discarded")
  }
  # If the parameter estimates have not been provided in pars then extract
  # them from the fitted object
  if (is.null(pars)) {
    pars <- coef(object)
  }
  n_pars <- length(pars)
  # Threshold exceedances (values that lie above the threshold)
  response_data <- object$data
  # If trans = FALSE then there are no covariates and object$data contains
  # the response data
  sigma <- pars["beta"]
  xi <- pars["xi"]
  # Calculate the loglikelihood contributions
  if (any(sigma <= 0)) {
    val <- -Inf
  } else {
    val <- revdbayes::dgp(response_data, loc = object$threshold, scale = sigma,
                          shape = xi, log = TRUE)
  }
  # Return the usual attributes for a "logLik" object
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- n_pars
  class(val) <- "logLikVec"
  return(val)
}

#' @export
nobs.evir_gpd <- function(object, ...) {
  return(object$n.exceed)
}

#' @export
coef.evir_gpd <- function(object, ...) {
  return(object$par.ests)
}

#' @export
vcov.evir_gpd <- function(object, ...) {
  vc <- object$varcov
  par_names <- names(coef(object))
  dimnames(vc) <- list(par_names, par_names)
  return(vc)
}

#' @export
logLik.evir_gpd <- function(object, ...) {
  return(logLik(logLikVec(object)))
}

# See evir_methods.R for nobs, coef, vcov, logLik methods for class "gpd"
