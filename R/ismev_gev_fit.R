# ============================== ismev::gev.fit ============================= #

# Methods for class ismev_gev

#' @export
logLikVec.ismev_gev <- function(object, pars = NULL, ...) {
  if (!missing(...)) {
    warning("extra arguments discarded")
  }
  # If the parameter estimates have not been provided in pars then extract
  # them from the fitted object
  if (is.null(pars)) {
    pars <- coef(object)
  }
  n_pars <- length(pars)
  #
  if (object$trans & is.null(object$xdat)) {
    stop("Covariate data are needed.  Refit the model using oolax::oogev.fit")
  }
  if (!object$trans) {
    response_data <- object$data
    # If trans = FALSE then there are no covariates and object$data contains
    # the response data
    mu <- pars[1]
    sigma <- pars[2]
    xi <- pars[3]
  } else {
    # If trans = TRUE then there are covariates
    response_data <- object$xdat
    # The numbers of parameters for mu, sigma, xi
    reg_pars <- sapply(object$model, length)
    npmu <- reg_pars[1] + 1
    npsc <- reg_pars[2] + 1
    npsh <- reg_pars[3] + 1
    # object$mumat, object$sigmat, object$shmat contain design matrices
    # Values of mu, sigma, xi for each observation
    mu <- object$mulink(object$mumat %*% (pars[1:npmu]))
    sigma <- object$siglink(object$sigmat %*%
                              (pars[seq(npmu + 1, length = npsc)]))
    xi <- object$shlink(object$shmat %*%
                          (pars[seq(npmu + npsc + 1, length = npsh)]))
  }
  # Calculate the loglikelihood contributions
  if (any(sigma <= 0)) {
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
nobs.ismev_gev <- function(object, ...) {
  return(length(object$data))
}

#' @export
coef.ismev_gev <- function(object, ...) {
  return(object$mle)
}

#' @export
vcov.ismev_gev <- function(object, ...) {
  return(object$cov)
}

#' @export
logLik.ismev_gev <- function(object, ...) {
  val <- -object$nllh
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(coef(object))
  class(val) <- "logLik"
  return(val)
}

# Methods for class gev.fit

#' @export
nobs.gev.fit <- function(object, ...) {
  return(length(object$data))
}

#' @export
coef.gev.fit <- function(object, ...) {
  return(object$mle)
}

#' @export
vcov.gev.fit <- function(object, ...) {
  return(object$cov)
}

#' @export
logLik.gev.fit <- function(object, ...) {
  val <- -object$nllh
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(coef(object))
  class(val) <- "logLik"
  return(val)
}
