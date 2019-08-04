# ============================== ismev::gev.fit ============================= #

# Methods for class ismev_gev

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
  if (object$trans & is.null(object$xdat)) {
    stop("Covariate data needed.  Refit the model using laxx::gev_refit")
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
    val <- revdbayes::dgev(response_data, loc = mu, scale = sigma,
                           shape = xi, log = TRUE)
  }
  # Return the usual attributes for a "logLik" object
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- n_pars
  class(val) <- "logLikVec"
  return(val)
}

#' @export
nobs.ismev_gev <- function(object, ...) {
  return(length(object$data))
}

#' @export
coef.ismev_gev <- function(object, ...) {
  val <- object$mle
  names(val) <- ismev_gev_names(object)
  return(val)
}

#' @export
vcov.ismev_gev <- function(object, ...) {
  vc <- object$cov
  dimnames(vc) <- list(ismev_gev_names(object), ismev_gev_names(object))
  return(vc)
}

#' @export
logLik.ismev_gev <- function(object, ...) {
  return(logLik(logLikVec(object)))
}

ismev_gev_names <- function(x) {
  if (x$trans) {
    if (is.null(colnames(x$ydat))) {
      loc_names <- paste0("loc", c("", x$model[[1]]))
      scale_names <- paste0("scale", c("", x$model[[2]]))
      shape_names <- paste0("shape", c("", x$model[[3]]))
    } else {
      cov_names <- colnames(x$ydat)
      loc_names <- paste0("loc", c("", cov_names[x$model[[1]]]))
      scale_names <- paste0("scale", c("", cov_names[x$model[[2]]]))
      shape_names <- paste0("shape", c("", cov_names[x$model[[3]]]))
    }
    val <- c(loc_names, scale_names, shape_names)
  } else {
    val <- c("loc", "scale", "shape")
  }
  return(val)
}

# See ismev_methods.R for nobs, coef, vcov, logLik methods for class "gev.fit"
