# =============================== ismev::pp.fit ============================= #

# Methods for class ismev_pp

#' @export
logLikVec.ismev_pp <- function(object, pars = NULL, ...) {
  # object$data only contains the exceedances.  We need all the data.
  # We could pack it with -Inf for non-exceedances but if cluster is not NULL
  # then we need the dat to be in the correct order.
  if (is.null(object$xdat)) {
    stop("Please refit the model using lax::pp_refit")
  }
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
  if (!object$trans) {
    # If trans = FALSE then there are no covariates and object$data contains
    # the response data
    mu <- pars[1]
    sigma <- pars[2]
    xi <- pars[3]
  } else {
    # If trans = TRUE then there are covariates
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
    pp_loglik_vec <- function(x, u, mu, sigma, xi) {
      logFu <- revdbayes::pgev(q = u, loc = mu, scale = sigma, shape = xi,
                               log.p = TRUE)
      logFx <- revdbayes::pgev(q = x, loc = mu, scale = sigma, shape = xi,
                               log.p = TRUE)
      logfx <- revdbayes::dgev(x = x, loc = mu, scale = sigma,
                               shape = xi, log = TRUE)
      rate_term <-  logFu / object$npy
      exc_term <- ifelse(x > u, logfx - logFx, 0)
      return(rate_term + exc_term)
    }
    val <- pp_loglik_vec(x = response_data, u = object$threshold, mu = mu,
                         sigma = sigma, xi = xi)
  }
  # Return the usual attributes for a "logLik" object
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- n_pars
  class(val) <- "logLikVec"
  return(val)
}

#' @export
nobs.ismev_pp <- function(object, ...) {
  return(nrow(object$vals))
}

#' @export
coef.ismev_pp <- function(object, ...) {
  val <- object$mle
  names(val) <- ismev_gev_names(object)
  return(val)
}

#' @export
vcov.ismev_pp <- function(object, ...) {
  vc <- object$cov
  dimnames(vc) <- list(ismev_gev_names(object), ismev_gev_names(object))
  return(vc)
}

#' @export
logLik.ismev_pp <- function(object, ...) {
  return(logLik(logLikVec(object)))
}

# See ismev_methods.R for nobs, coef, vcov, logLik methods for class "pp.fit"
