# ================================ evd::fgev ================================ #

# Methods for class POT_fitgpd

#' @export
logLikVec.POT_fitgpd <- function(object, pars = NULL, ...) {
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
  sigma <- all_pars["scale"]
  xi <- all_pars["shape"]
  # Calculate the loglikelihood contributions
  if (sigma <= 0) {
    val <- -Inf
  } else {
    val <- revdbayes::dgp(object$exceed, loc = object$threshold,
                          scale = sigma, shape = xi, log = TRUE)
  }
  # Return the usual attributes for a "logLik" object
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(free_pars)
  class(val) <- "logLikVec"
  return(val)
}

#' @export
nobs.POT_fitgpd <- function(object, ...) {
  return(object$nat)
}

#' @export
coef.POT_fitgpd <- function(object, complete = FALSE, ...) {
  if (complete) {
    val <- object$param
  } else {
    val <- object$fitted.values
  }
  return(val)
}

#' @export
vcov.POT_fitgpd <- function(object, complete = FALSE, ...) {
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
logLik.POT_fitgpd <- function(object, ...) {
  return(logLik(logLikVec(object)))
}
