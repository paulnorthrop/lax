# ============================== texmex::evmt =============================== #

# Methods for class texmex_evmOpt

#' @export
logLikVec.texmex_evmOpt <- function(object, pars = NULL, ...) {
  if (!missing(...)) {
    warning("extra arguments discarded")
  }
  # If the parameter estimates have not been provided in pars then extract
  # them from the fitted object
  if (is.null(pars)) {
    pars <- coef(object)
  }
  n_pars <- length(pars)
  response_data <- object$data$y
  # Determine the model and whether or not there are covariates
  if (object$family$name == "GEV") {
    if (n_pars == 3) {
      mu <- pars["mu: (Intercept)"]
      phi <- pars["phi: (Intercept)"]
      xi <- pars["xi: (Intercept)"]
    } else {
      mu_mat <- object$data$D$mu
      phi_mat <- object$data$D$phi
      xi_mat <- object$data$D$xi
      n_mu <- ncol(mu_mat)
      n_phi <- ncol(phi_mat)
      n_xi <- ncol(xi_mat)
      mu_pars <- pars[1:n_mu]
      phi_pars <- pars[(n_mu + 1):(n_mu + n_phi)]
      xi_pars <- pars[(n_mu + n_phi + 1):n_pars]
      mu <- as.vector(mu_mat %*% mu_pars)
      phi <- as.vector(phi_mat %*% phi_pars)
      xi <- as.vector(xi_mat %*% xi_pars)
    }
    # Calculate the loglikelihood contributions
    sigma <- exp(phi)
    if (any(sigma <= 0)) {
      val <- -Inf
    } else {
      val <- revdbayes::dgev(response_data, loc = mu, scale = sigma,
                             shape = xi, log = TRUE)
    }
  } else if (object$family$name == "GPD") {
    if (n_pars == 2) {
      phi <- pars["phi: "]
      xi <- pars["xi: "]
    } else {
      phi_mat <- object$data$D$phi
      xi_mat <- object$data$D$xi
      n_phi <- ncol(phi_mat)
      n_xi <- ncol(xi_mat)
      phi_pars <- pars[1:n_phi]
      xi_pars <- pars[(n_phi + 1):n_pars]
      phi <- as.vector(phi_mat %*% phi_pars)
      xi <- as.vector(xi_mat %*% xi_pars)
    }
    sigma <- exp(phi)
    # Calculate the loglikelihood contributions
    if (any(sigma <= 0)) {
      val <- -Inf
    } else {
      val <- revdbayes::dgp(response_data, loc = object$threshold,
                            scale = sigma, shape = xi, log = TRUE)
    }
  }
  # Return the usual attributes for a "logLik" object
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(pars)
  class(val) <- "logLikVec"
  return(val)
}

#' @export
nobs.texmex_evmOpt <- function(object, ...) {
  return(length(object$data$y))
}

#' @export
vcov.texmex_evmOpt <- function(object, ...) {
  vc <- object$cov
  par_names <- names(coef(object))
  dimnames(vc) <- list(par_names, par_names)
  return(vc)
}

#' @export
coef.texmex_evmOpt <- function(object, ...) {
  return(object$coefficients)
}

#' @export
logLik.texmex_evmOpt <- function(object, ...) {
  return(logLik(logLikVec(object)))
}

# See texmex_methods.R for nobs and vcov methods for class "evmOpt"
# (texmex already has coef and logLik methods)
