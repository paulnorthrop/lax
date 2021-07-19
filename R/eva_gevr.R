# =============================== eva::gevrFit ============================== #

# Methods for class laxeva_gevr

#' @export
logLikVec.laxeva_rlarg <- function(object, pars = NULL, ...) {
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
  # Response data r-column matrix of r-largest order statistics
  response_data <- object$data
  # Design matrices for the location, scale and shape parameters
  mumat <- object$covars[[1]]
  sigmamat <- object$covars[[2]]
  ximat <- object$covars[[3]]
  # Link functions
  mulink <- object$links[[1]]
  sigmalink <- object$links[[2]]
  xilink <- object$links[[3]]
  # Numbers of parameter estimates for mu, sigma and xi
  npmu <- object$parnum[1]
  npsc <- object$parnum[2]
  npsh <- object$parnum[3]
  # object$mumat, object$sigmat, object$shmat contain design matrices
  # Values of mu, sigma, xi for each observation
  mu <- mulink(mumat %*% (pars[1:npmu]))
  sigma <- sigmalink(sigmamat %*% (pars[seq(npmu + 1, length = npsc)]))
  xi <- xilink(ximat %*% (pars[seq(npmu + npsc + 1, length = npsh)]))
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
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- n_pars
  class(val) <- "logLikVec"
  return(val)
}

#' @export
nobs.laxeva_rlarg <- function(object, ...) {
  return(object$n)
}

#' @export
coef.laxeva_rlarg <- function(object, ...) {
  return(object$par.ests)
}

#' @export
vcov.laxeva_rlarg <- function(object, ...) {
  vc <- object$varcov
  dimnames(vc) <- list(names(coef(object)), names(coef(object)))
  return(vc)
}

#' @export
logLik.laxeva_rlarg <- function(object, ...) {
  val <- -object$nllh.final
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(coef(object))
  class(val) <- "logLik"
  return(val)
}
