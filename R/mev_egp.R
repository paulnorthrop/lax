# =============================== mev::fit.egp ============================== #

# Method for class laxmev_egp

#' @export
logLikVec.laxmev_egp <- function(object, pars = NULL, ...) {
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
  n <- length(response_data)
  # Calculate the loglikelihood contributions
  if (pars[2] <= 0) {
    val <- -Inf
  } else {
    fn <- function(i) {
      return(mev::egp.ll(xdat = response_data[i], thresh = object$thresh,
                         par = pars, model = object$model))
    }
    val <- vapply(1:n, fn, 0.0)
  }
  # Return the usual attributes for a "logLik" object
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- n_pars
  class(val) <- "logLikVec"
  return(val)
}

#' @export
nobs.laxmev_egp <- function(object, ...) {
  return(object$nat)
}

#' @export
coef.laxmev_egp <- function(object, ...) {
  return(object$estimate)
}

#' @export
vcov.laxmev_egp <- function(object, ...) {
  return(object$vcov)
}

#' @export
logLik.laxmev_egp <- function(object, ...) {
  val <- -object$deviance / 2
  attr(val, "names") <- NULL
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(coef(object))
  class(val) <- "logLik"
  return(val)
}
