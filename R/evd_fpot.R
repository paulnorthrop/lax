# ================================ evd::fpot ================================ #

# Methods for class evd_fpot

#' @export
logLikVec.evd_fpot <- function(object, pars = NULL, ...) {
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
  n_all_pars <- length(all_pars)
  # If n_all_pars = 2 then model = "gp".
  # If n_all_pars = 3 then model = "pp".
  if (n_all_pars == 2) {
    sigma <- all_pars["scale"]
    xi <- all_pars["shape"]
    # Calculate the loglikelihood contributions
    if (sigma <= 0) {
      val <- -Inf
    } else {
      val <- revdbayes::dgp(object$exceedances, loc = object$threshold,
                            scale = sigma, shape = xi, log = TRUE)
    }
  } else {
    mu <- all_pars["loc"]
    sigma <- all_pars["scale"]
    xi <- all_pars["shape"]
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
        exc_term <- ifelse(x > u, logfx - logFx, 0)
        return(rate_term + exc_term)
      }
      val <- pp_loglik_vec(x = object$data, u = object$threshold, mu = mu,
                           sigma = sigma, xi = xi)
    }
  }
  # Return the usual attributes for a "logLik" object
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(free_pars)
  class(val) <- "logLikVec"
  return(val)
}

#' @export
nobs.evd_fpot <- function(object, ...) {
  if (length(object$param) == 2) {
    val <- object$nat
  } else if (length(object$param) == 3) {
    val <- length(object$data)
  } else {
    val <- NA
  }
  return(val)
}

#' @export
coef.evd_fpot <- function(object, complete = FALSE, ...) {
  if (complete) {
    val <- object$param
  } else {
    val <- object$estimate
  }
  return(val)
}

#' @export
vcov.evd_fpot <- function(object, complete = FALSE, ...) {
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
logLik.evd_fpot <- function(object, ...) {
  return(logLik(logLikVec(object)))
}

# Methods for class pot (evd already has vcov and logLik methods)

#' @export
nobs.pot <- function(object, ...) {
  if (length(object$param) == 2) {
    val <- object$nat
  } else if (length(object$param) == 3) {
    val <- length(object$data)
  } else {
    val <- NA
  }
  return(val)
}

#' @export
coef.pot <- function(object, complete = FALSE, ...) {
  if (complete) {
    val <- object$param
  } else {
    val <- object$estimate
  }
  return(val)
}
