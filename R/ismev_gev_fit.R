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
  if (object$trans) {
    stop("Covariates are not allowed")
  }
  # If trans = FALSE then there are no covariates and object$data contains
  # the response data
  mu <- pars[1]
  sigma <- pars[2]
  xi <- pars[3]
  # Calculate the weighted loglikelihood contributions
  if (sigma <= 0) {
    val <- -Inf
  } else {
    val <- revdbayes::dgev(object$data, loc = mu, scale = sigma,
                           shape = xi, log = TRUE)
  }
  return(val)
}

#' @export
nobs.ismev_gev <- function(object) {
  return(length(object$data))
}

#' @export
coef.ismev_gev <- function(object) {
  return(object$mle)
}

#' @export
vcov.ismev_gev <- function(object, ...) {
  return(object$cov)
}
