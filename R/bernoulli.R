#' Inference for the Bernoulli distribution
#'
#' Functions involved in making inferences about the probability of success
#' in a Bernoulli distribution.
#'
#' @param data A numeric vector of outcomes from Bernoulli trials: 0 for a
#'    failure, 1 for a success.  Alternatively, a logical vector with FALSE
#'    for a failure and TRUE for a success.
#' @param x,object A fitted model object returned from \code{fit_bernoulli()}.
#' @param pars A numeric parameter vector of length 1 containing the value of
#'   the Bernoulli success probability.
#' @param cluster A vector or factor indicating from which cluster each
#'   observation in \code{data} originates.
#' @param use_vcov A logical scalar.  Should we use the \code{vcov} S3 method
#'   for \code{x} (if this exists) to estimate the Hessian of the independence
#'   loglikelihood to be passed as the argument \code{H} to
#'   \code{\link[chandwich]{adjust_loglik}}?
#'   Otherwise, \code{H} is estimated inside
#'   \code{\link[chandwich]{adjust_loglik}} using
#'   \code{\link[stats:optim]{optimHess}}.
#' @param ... Further arguments to be passed to the functions in the
#'   sandwich package \code{\link[sandwich]{meat}} (if \code{cluster = NULL}),
#'   or \code{\link[sandwich:vcovCL]{meatCL}} (if \code{cluster} is not
#'   \code{NULL}).
#' @details
#' \code{fit_bernoulli}: fit a Bernoulli distribution
#'
#' \code{logLikVec.bernoulli}: calculates contributions to a loglikelihood based
#' on the Bernoulli distribution.  The loglikelihood is calculated up to an
#' additive constant.
#'
#' \code{nobs, coef, vcov} and \code{logLik} methods are provided.
#' @return
#' \code{fit_bernoulli} returns an object of class \code{"bernoulli"}, a list
#' with components: \code{logLik, mle, nobs, vcov, data, obs_data}, where
#' \code{data} are the input data and \code{obs_data} are the input data after
#' any missing values have been removed, using
#' \code{\link[stats:na.fail]{na.omit}}.
#'
#' \code{logLikVec.bernoulli} returns an object of class \code{"logLikVec"}, a
#' vector length \code{length(data)} containing the likelihood contributions
#' from the individual observations in \code{data}.
#' @seealso \code{\link[stats]{Binomial}}.  The Bernoulli distribution is the
#'   special case where \code{size = 1}.
#' @examples
#' # Set up data
#' x <- exdex::newlyn
#' u <- quantile(x, probs = 0.9)
#' exc <- x > u
#'
#' # Fit a Bernoulli distribution
#' fit <- fit_bernoulli(exc)
#'
#' # Calculate the loglikelihood at the MLE
#' res <- logLikVec(fit)
#'
#' # The logLik method sums the individual loglikelihood contributions.
#' logLik(res)
#'
#' # nobs, coef, vcov, logLik methods for objects returned from fit_bernoulli()
#' nobs(fit)
#' coef(fit)
#' vcov(fit)
#' logLik(fit)
#'
#' # Adjusted loglikelihood
#' # Create 5 clusters each corresponding approximately to 1 year of data
#' cluster <- rep(1:5, each = 579)[-1]
#' afit <- alogLik(fit, cluster = cluster, cadjust = FALSE)
#' summary(afit)
#' @name bernoulli
NULL
## NULL

#' @rdname bernoulli
#' @export
fit_bernoulli <- function(data) {
  res <- list()
  res$data <- data
  obs_data <- stats::na.omit(data)
  res$obs_data <- obs_data
  res$mle <- mean(obs_data)
  res$nobs <- length(obs_data)
  res$vcov <- as.matrix(res$mle * (1 - res$mle) / res$nobs)
  n1 <- sum(obs_data)
  n0 <- res$nobs - n1
  res$logLik <- n1 * log(res$mle) + n0 * log(1 - res$mle)
  class(res) <- "bernoulli"
  return(res)
}

#' @rdname bernoulli
#' @export
logLikVec.bernoulli <- function(object, pars = NULL, ...) {
  if (!missing(...)) {
    warning("extra arguments discarded")
  }
  # If the parameter estimates have not been provided in pars then extract
  # them from the fitted object
  if (is.null(pars)) {
    pars <- coef(object)
  }
  n_pars <- length(pars)
  prob <- pars[1]
  if (prob < 0 || prob > 1) {
    val <- -Inf
  } else {
    val <- stats::dbinom(object$obs_data, 1, prob, log = TRUE)
  }
  # Return the usual attributes for a "logLik" object
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- n_pars
  class(val) <- "logLikVec"
  return(val)
}

#' @rdname bernoulli
#' @export
nobs.bernoulli <- function(object, ...) {
  return(object$nobs)
}

#' @rdname bernoulli
#' @export
coef.bernoulli <- function(object, ...) {
  val <- object$mle
  names(val) <- "prob"
  return(val)
}

#' @rdname bernoulli
#' @export
vcov.bernoulli <- function(object, ...) {
  vc <- object$vcov
  par_names <- "prob"
  dimnames(vc) <- list(par_names, par_names)
  return(vc)
}

#' @rdname bernoulli
#' @export
logLik.bernoulli <- function(object, ...) {
  val <- object$logLik
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(coef(object))
  class(val) <- "logLik"
  return(val)
}

#' @rdname bernoulli
#' @export
alogLik.bernoulli <- function(x, cluster = NULL, use_vcov = TRUE, ...) {
  # Call adj_object() to adjust the loglikelihood
  res <- adj_object(x, cluster = cluster, use_vcov = use_vcov, ...)
  class(res) <- c("lax", "chandwich", "bernoulli", "stat")
  return(res)
}
