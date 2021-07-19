# ================================ eva::gpdFit ============================== #

# Methods for class laxeva_gpd

#' @export
logLikVec.laxeva_gpd <- function(object, pars = NULL, ...) {
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
  # Threshold exceedances (values that lie above the threshold)
  response_data <- object$exceedances
  # Design matrices for the scale and shape parameters
  sigmamat <- object$covars[[1]]
  ximat <- object$covars[[2]]
  # Link functions
  sigmalink <- object$links[[1]]
  xilink <- object$links[[2]]
  # Parameter estimates
  sigmapars <- pars[1:object$parnum[1]]
  xipars <- pars[(object$parnum[1] + 1):(object$parnum[1] + object$parnum[2])]
  # object$sigmat, object$shmat contain design matrices
  # Values of sigma, xi for each observation
  sigma <- sigmalink(sigmamat %*% sigmapars)
  xi <- xilink(ximat %*% xipars)
  # Calculate the loglikelihood contributions
  if (any(sigma <= 0)) {
    val <- -Inf
  } else {
    val <- revdbayes::dgp(response_data, loc = object$threshold, scale = sigma,
                          shape = xi, log = TRUE)
  }
  # Return the usual attributes for a "logLik" object
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- n_pars
  class(val) <- "logLikVec"
  return(val)
}

#' @export
nobs.laxeva_gpd <- function(object, ...) {
  return(object$n.exceed)
}

#' @export
coef.laxeva_gpd <- function(object, ...) {
  return(object$par.ests)
}

#' @export
vcov.laxeva_gpd <- function(object, ...) {
  vc <- object$varcov
  dimnames(vc) <- list(names(coef(object)), names(coef(object)))
  return(vc)
}

#' @export
logLik.laxeva_gpd <- function(object, ...) {
  val <- -object$nllh.final
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(coef(object))
  class(val) <- "logLik"
  return(val)
}
