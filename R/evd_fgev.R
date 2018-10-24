# ================================ evd::fgev ================================ #

#' @export
logLikVec.evd_fgev <- function(object, contrib = FALSE, ...) {
  if (!missing(...)) {
    warning("extra arguments discarded")
  }
  # Extract the parameter estimates
  pars <- object$estimate
  n_pars <- length(pars)
  mu <- pars[1]
  sigma <- pars[n_pars - 1]
  xi <- pars[n_pars]
  if (n_pars > 3) {
    mu_reg <- pars[2:(n_pars - 2)]
    mu <- mu + as.matrix(object$nsloc) %*% mu_reg
  }
  # Calculate the weighted loglikelihood contributions
  if (sigma <= 0) {
    val <- -Inf
  } else {
    val <- revdbayes::dgev(object$data, loc = as.vector(mu), scale = sigma,
                           shape = xi, log = TRUE)
  }
  # Sum them if we want the overall loglikelihood
  # ... and return the usual logLik object
  if (!contrib) {
    val <- sum(val)
    attr(val, "nobs") <- length(object$data)
    attr(val, "df") <- n_pars
    class(val) <- "logLik"
  }
  return(val)
}

##' @export
#logLik.evd_fgev <- function(object, contrib = FALSE, ...) {
#  if (!missing(...)) {
#    warning("extra arguments discarded")
#  }
#  # Extract the parameter estimates
#  pars <- object$estimate
#  n_pars <- length(pars)
#  mu <- pars[1]
#  sigma <- pars[n_pars - 1]
#  xi <- pars[n_pars]
#  if (n_pars > 3) {
#    mu_reg <- pars[2:(n_pars - 2)]
#    mu <- mu + as.matrix(object$nsloc) %*% mu_reg
#  }
#  # Calculate the loglikelihood contributions
#  if (sigma <= 0) {
#    val <- -Inf
#  } else {
#    val <- revdbayes::dgev(object$data, loc = as.vector(mu), scale = sigma,
#                           shape = xi, log = TRUE)
#  }
#  # Sum them if we want the overall loglikelihood
#  # ... and return the usual logLik object
#  if (!contrib) {
#    val <- sum(val)
#    attr(val, "nobs") <- length(object$data)
#    attr(val, "df") <- n_pars
#    class(val) <- "logLik"
#  }
#  return(val)
#}

##' @export
#logLikFn.evd_fgev <- function(pars, fitted_object, contrib = TRUE) {
#  new_object <- fitted_object
#  new_object$estimate <- pars
#  return(logLik(new_object, contrib = contrib))
#}

#' @export
nobs.evd_fgev <- function(object) {
  return(length(object$data))
}

#' @export
coef.evd_fgev <- function(object) {
  return(object$estimate)
}

#' @export
vcov.evd_fgev <- function(x, ...) {
  class(x) <- "evd"
  return(vcov(x))
}
