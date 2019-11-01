# ============================= mev::rlarg.fit ============================ #

# logLikVec method for class laxmev_rlarg

#' @export
logLikVec.laxmev_rlarg <- function(object, pars = NULL, ...) {
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
  response_data <- object$xdat
  # If trans = FALSE then there are no covariates and object$data contains
  # the response data
  mu <- pars[1]
  sigma <- pars[2]
  xi <- pars[3]
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
nobs.laxmev_rlarg <- function(object, ...) {
  return(nrow(object$xdat))
}

#' @export
coef.laxmev_rlarg <- function(object, ...) {
  return(object$estimate)
}

#' @export
vcov.laxmev_rlarg <- function(object, ...) {
  vc <- object$vcov
  dimnames(vc) <- list(names(coef(object)), names(coef(object)))
  return(vc)
}

#' @export
logLik.laxmev_rlarg <- function(object, ...) {
  val <- -object$nllh
  attr(val, "names") <- NULL
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(coef(object))
  class(val) <- "logLik"
  return(val)
}
