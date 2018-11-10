# ================================ evd::fgev ================================ #

# Methods for class evd_fgev
# The returned object has class c("gev", "uvdata", "evd")

#' @export
logLikVec.evd_fgev <- function(object, pars = NULL, ...) {
  if (!missing(...)) {
    warning("extra arguments discarded")
  }
  # Extract from object all the parameter estimates: free and fixed
  all_pars <- coef(object, complete = TRUE)
  free_pars <- coef(object, complete = FALSE)
  # If pars is supplied then overwrite the values of the free parameters
  if (!is.null(pars)) {
    names_all_pars <- names(all_pars)
    names_free_pars <- names(free_pars)
    which_free <- which(all_pars %in% free_pars)
    all_pars[which_free] <- pars
  }
  mu <- all_pars["loc"]
  sigma <- all_pars["scale"]
  xi <- all_pars["shape"]
  if (!is.null(object$nsloc)) {
    mu_reg <- all_pars[paste0("loc", names(object$nsloc))]
    mu <- mu + as.matrix(object$nsloc) %*% mu_reg
  }
  # Calculate the loglikelihood contributions
  if (sigma <= 0) {
    val <- -Inf
  } else {
    val <- revdbayes::dgev(object$data, loc = as.vector(mu), scale = sigma,
                           shape = xi, log = TRUE)
  }
  # Return the usual attributes for a "logLik" object
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(free_pars)
  class(val) <- "logLikVec"
  return(val)
}

#' @export
nobs.evd_fgev <- function(object, ...) {
  return(object$n)
}

#' @export
coef.evd_fgev <- function(object, complete = FALSE, ...) {
  if (complete) {
    val <- object$param
  } else {
    val <- object$estimate
  }
  return(val)
}

#' @export
vcov.evd_fgev <- function(object, complete = FALSE, ...) {
  # Variance-covariance matrix for the free parameters
  vc <- object$var.cov
  free_pars <- names(coef(object, complete = FALSE))
  dimnames(vc) <- list(free_pars, free_pars)
  if (complete) {
    all_pars <- names(coef(object, complete = TRUE))
    np <- length(all_pars)
    dummy <- matrix(NA, np, np)
    which_free <- which(all_pars %in% free_pars)
    dummy[which_free, which_free] <- vc
    vc <- dummy
    dimnames(vc) <- list(all_pars, all_pars)
  }
  return(vc)
}

#' @export
logLik.evd_fgev <- function(object, ...) {
  return(logLik(logLikVec(object)))
}

# See evd_methods.R for nobs and coef methods for class "evd"
# (evd already has vcov and logLik methods)
