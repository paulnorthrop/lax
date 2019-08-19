# ================================== evir::pot ============================== #

# Methods for class evir_pot

#' @export
logLikVec.evir_pot <- function(object, pars = NULL, ...) {
  if (is.null(object$input_data)) {
    stop("Input data are needed.  Please refit the model using lax::pot_refit")
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
  response_data <- object$input_data
  npy <- length(response_data) / object$span
  mu <- pars["mu"]
  sigma <- pars["sigma"]
  xi <- pars["xi"]
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
      rate_term <-  logFu / npy
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
nobs.evir_pot <- function(object, ...) {
  return(object$n)
}

#' @export
coef.evir_pot <- function(object, complete = FALSE, ...) {
  if (complete) {
    val <- object$par.ests
  } else {
    val <- object$par.ests[c("xi", "sigma", "mu")]
  }
  return(val)
}

#' @export
vcov.evir_pot <- function(object, complete = FALSE, ...) {
  vc <- object$varcov
  par_names <- names(coef(object))
  if (complete) {
    vc <- rbind(cbind(vc, NA), NA)
    par_names <- c(names(coef(object)), "beta")
  }
  dimnames(vc) <- list(par_names, par_names)
  return(vc)
}

#' @export
logLik.evir_pot <- function(object, ...) {
  val <- -object$nllh.final
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(coef(object))
  class(val) <- "logLik"
  return(val)
}

# See evir_methods.R for nobs, coef, vcov, logLik methods for class "potd"
