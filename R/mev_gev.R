# =============================== mev::mev_gev ============================= #

# Methods for class laxmev_gev

#' @export
logLikVec.laxmev_gev <- function(object, pars = NULL, ...) {
  if (!missing(...)) {
    warning("extra arguments discarded")
  }
  # If the parameter estimates have not been provided in pars then extract
  # them from the fitted object
  if (is.null(pars)) {
    pars <- coef(object)
  }
  n_pars <- length(pars)
  response_data <- object$xdat
  mu <- pars[1]
  sigma <- pars[2]
  xi <- pars[3]
  # Calculate the loglikelihood contributions
  if (sigma <= 0) {
    val <- -Inf
  } else {
    val <- revdbayes::dgev(response_data, loc = mu, scale = sigma,
                           shape = xi, log = TRUE)
  }
  # Return the usual attributes for a "logLik" object
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- n_pars
  class(val) <- "logLikVec"
  return(val)
}

#' @export
nobs.laxmev_gev <- function(object, ...) {
  return(length(object$xdat))
}

#' @export
coef.laxmev_gev <- function(object, ...) {
  return(object$estimate)
}

#' @export
vcov.laxmev_gev <- function(object, ...) {
  return(object$vcov)
}

#' @export
logLik.laxmev_gev <- function(object, ...) {
  val <- -object$nllh
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(coef(object))
  class(val) <- "logLik"
  return(val)
}
