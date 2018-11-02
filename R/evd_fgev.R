# ================================ evd::fgev ================================ #

# Methods for class evd_fgev

#' @export
logLikVec.evd_fgev <- function(object, pars = NULL, ...) {
  if (!missing(...)) {
    warning("extra arguments discarded")
  }
  # If the parameter estimates have not been provided in pars then extract
  # them from the fitted object
  if (is.null(pars)) {
    pars <- coef(object)
  }
  n_pars <- length(pars)
  mu <- pars[1]
  sigma <- pars[n_pars - 1]
  xi <- pars[n_pars]
  if (n_pars > 3) {
    mu_reg <- pars[2:(n_pars - 2)]
    mu <- mu + as.matrix(object$nsloc) %*% mu_reg
  }
  # Calculate the loglikelihood contributions
  if (sigma <= 0) {
    val <- -Inf
  } else {
    val <- evd::dgev(object$data, loc = as.vector(mu), scale = sigma,
                     shape = xi, log = TRUE)
  }
  # Return the usual attributes for a "logLik" object
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- n_pars
  class(val) <- "logLikVec"
  return(val)
}

#' @export
nobs.evd_fgev <- function(object, ...) {
  return(object$n)
}

#' @export
coef.evd_fgev <- function(object, ...) {
  return(object$estimate)
}

#' @export
vcov.evd_fgev <- function(object, ...) {
  return(object$var.cov)
}

#' @export
logLik.evd_fgev <- function(object, ...) {
  return(logLik(logLikVec(object)))
}

# Methods for class gev (evd already has vcov and logLik methods)

#' @export
nobs.gev <- function(object, ...) {
  return(object$n)
}

#' @export
coef.gev <- function(object, ...) {
  return(object$estimate)
}

