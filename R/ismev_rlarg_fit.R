# ============================= ismev::rlarg.fit ============================ #

# Methods for class ismev_rlarg

#' @export
logLikVec.ismev_rlarg <- function(object, pars = NULL, ...) {
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
    stop("Covariate data needed.  Refit the model using lax::rlarg_refit")
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
    rlarg_loglik_vec <- function(x, mu, sigma, xi) {
      logg <- apply(x, 2, revdbayes::dgev, loc = mu, scale = sigma,
                    shape = xi, log = TRUE)
      logG <- apply(x, 2, revdbayes::pgev, loc = mu, scale = sigma,
                    shape = xi, log.p = TRUE)
      logGmin <- revdbayes::pgev(min_response, loc = mu, scale = sigma,
                                 shape = xi, log.p = TRUE)
      loglik <- logGmin + rowSums(logg - logG, na.rm = TRUE)
      return(loglik)
    }
    min_response <- apply(response_data, 1, min, na.rm = TRUE)
    val <- rlarg_loglik_vec(x = response_data, mu = mu, sigma = sigma, xi = xi)
  }
  # Return the usual attributes for a "logLik" object
  attr(val, "names") <- NULL
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- n_pars
  class(val) <- "logLikVec"
  return(val)
}

#' @export
nobs.ismev_rlarg <- function(object, ...) {
  return(nrow(object$vals))
}

#' @export
coef.ismev_rlarg <- function(object, ...) {
  val <- object$mle
  names(val) <- ismev_gev_names(object)
  return(val)
}

#' @export
vcov.ismev_rlarg <- function(object, ...) {
  vc <- object$cov
  dimnames(vc) <- list(ismev_gev_names(object), ismev_gev_names(object))
  return(vc)
}

#' @export
logLik.ismev_rlarg <- function(object, ...) {
  return(logLik(logLikVec(object)))
}

# See ismev_methods.R for nobs, coef, vcov, logLik methods for class "rlarg.fit"
