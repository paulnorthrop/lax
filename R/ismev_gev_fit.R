# ================================ evd::fgev ================================ #

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
  if (object$trans & is.null(object$ydat)) {
    stop("Covariate data are needed.  Refit the model using oolax::oogev.fit")
  }
  if (!object$trans) {
    # If trans = FALSE then there are no covariates and object$data contains
    # the response data
    mu <- pars[1]
    sigma <- pars[2]
    xi <- pars[3]
    # Calculate the loglikelihood contributions
    if (sigma <= 0) {
      val <- -Inf
    } else {
      val <- revdbayes::dgev(object$data, loc = mu, scale = sigma,
                             shape = xi, log = TRUE)
    }
  } else {
    # If trans = TRUE then there are covariates, object$xdat contains the
    # response data and object$ydat the matrix of covariate data
    # object$model is a list containing the columns of ydat for mu, sigma, xi
    # Numbers of regression parameters for mu, sigma and xi
    reg_pars <- sapply(object$model, length)
    cum_reg_pars <- cumsum(reg_pars)
    mu <- pars[1]
    if (reg_pars[1] > 0) {
      mu_reg <- pars[2:(1 + reg_pars[1])]
      mu <- mu + object$ydat[, object$model[[1]], drop = FALSE] %*% mu_reg
    }
    sigma <- pars[2 + reg_pars[1]]
    if (reg_pars[2] > 0) {
      sigma_reg <- pars[(3 + reg_pars[1]):(2 + cum_reg_pars[1])]
      sigma <- sigma + object$ydat[, object$model[[2]], drop = FALSE] %*%
        sigma_reg
    }
    xi <- pars[3 + cum_reg_pars[2]]
    if (reg_pars[3] > 0) {
      xi_reg <- pars[(4 + cum_reg_pars[2]):(3 + cum_reg_pars[3])]
      xi <- xi + object$ydat[, object$model[[3]], drop = FALSE] %*% xi_reg
    }
    # Calculate the loglikelihood contributions
    if (any(sigma <= 0)) {
      val <- -Inf
    } else {
      val <- revdbayes::dgev(object$xdat, loc = mu, scale = sigma,
                             shape = xi, log = TRUE)
    }
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
