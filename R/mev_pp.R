# ================================= mev::fit.pp ============================= #

# Methods for class laxmev_pp

#' @export
logLikVec.laxmev_pp <- function(object, pars = NULL, ...) {
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
  # If trans = FALSE then there are no covariates and object$data contains
  # the response data
  mu <- pars[1]
  sigma <- pars[2]
  xi <- pars[3]
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
      rate_term <-  logFu / object$npp
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
nobs.laxmev_pp <- function(object, ...) {
  return(object$nat / object$pat)
}

#' @export
coef.laxmev_pp <- function(object, ...) {
  return(object$estimate)
}

#' @export
vcov.laxmev_pp <- function(object, ...) {
  return(object$vc)
}

#' @export
logLik.laxmev_pp <- function(object, ...) {
  val <- -object$nllh
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(coef(object))
  class(val) <- "logLik"
  return(val)
}
