# ================================ evd::fpot ================================ #

# Methods for class evd_fpot

#' @export
logLikVec.evd_fpot <- function(object, pars = NULL, ...) {
  if (!missing(...)) {
    warning("extra arguments discarded")
  }
  # If the parameter estimates have not been provided in pars then extract
  # them from the fitted object
  if (is.null(pars)) {
    pars <- coef(object)
  }
  n_pars <- length(pars)
  # If n_pars = 2 then model = "gp".
  # If n_pars = 3 then model = "pp".
  if (n_pars == 2) {
    sigma <- pars[1]
    xi <- pars[2]
    # Calculate the loglikelihood contributions
    if (sigma <= 0) {
      val <- -Inf
    } else {
      val <- revdbayes::dgp(object$exceedances, loc = object$threshold,
                            scale = sigma, shape = xi, log = TRUE)
    }
  } else {
    mu <- pars[1]
    sigma <- pars[2]
    xi <- pars[3]
    # Calculate the loglikelihood contributions
    if (sigma <= 0) {
      val <- -Inf
    } else {
      pp_loglik_vec <- function(x, u, mu, sigma, xi) {
        logFu <- revdbayes::pgev(q = u, loc = mu, scale = sigma, shape = xi,
                                 log.p = TRUE)
        logFx <- revdbayes::pgev(q = x, loc = mu, scale = sigma, shape = xi,
                                 log.p = TRUE)
        logfx <- revdbayes::dgev(x = x, loc = mu, scale = sigma,
                                 shape = xi, log = TRUE)
        rate_term <-  logFu / object$npp
        exc_term <- (x > u) * (logfx - logFx)
        return(rate_term + exc_term)
      }
      val <- pp_loglik_vec(x = object$data, u = object$threshold, mu = mu,
                           sigma = sigma, xi = xi)
    }
  }
  # Return the usual attributes for a "logLik" object
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- n_pars
  class(val) <- "logLikVec"
  return(val)
}

#' @export
nobs.evd_fpot <- function(object, ...) {
  return(object$nat)
}

#' @export
coef.evd_fpot <- function(object, ...) {
  return(object$estimate)
}

#' @export
vcov.evd_fpot <- function(object, ...) {
  return(object$var.cov)
}

#' @export
logLik.evd_fpot <- function(object, ...) {
  return(logLik(logLikVec(object)))
}

# Methods for class pot (evd already has vcov and logLik methods)

#' @export
nobs.pot <- function(object, ...) {
  return(object$nat)
}

#' @export
coef.pot <- function(object, ...) {
  return(object$estimate)
}

