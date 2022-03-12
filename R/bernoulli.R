#' Inference for the Bernoulli distribution
#'
#' Functions involved in making inferences about the probability of success
#' in a Bernoulli distribution.
#'
#' @param prob A numeric scalar. Probability of success on each trial.
#' @param data A numeric vector of outcomes from Bernoulli trials: 0 for a
#'    failure, 1 for a success.  Alternatively, a logical vector with FALSE
#'    for a failure and TRUE for a success.
#' @details
#' \code{fit_bernoulli}: fit a Bernoulli distribution
#'
#' \code{bernoulli_loglik}: calculates contributions to a log-likelihood based
#' on the Bernoulli distribution.  The log-likelihood is calculated up to an
#' additive constant.
#'
#' \code{nobs, coef, vcov} and \code{logLik} methods are provided.
#' @return
#' \code{fit_bernoulli} returns an object of class \code{"bernoulli"}, a list
#' with components: \code{data, logLik, mle, nobs, vcov}, where \code{data}
#' are the input data.
#'
#' \code{bernoulli_loglik} returns an object of class \code{"logLikVec"}, a
#' vector length \code{length(data)} containing the likelihood contributions
#' from the individual observations in \code{data}.
#' @seealso \code{\link[stats]{Binomial}}.  The Bernoulli distribution is the
#'   special case where \code{size = 1}.
#' @examples
#' # Set up data
#' x <- exdex::newlyn
#' u <- quantile(x, probs = 0.6)
#' exc <- x > u
#'
#' # Fit a Bernoulli distribution
#' fit <- fit_bernoulli(exc)
#'
#' # Calculate the log-likelihood at the MLE
#' res <- bernoulli_loglik(fit$mle, exc)
#'
#' # The logLik method sums the individual log-likelihood contributions.
#' logLik(res)
#'
#' # nobs, coef, vcov, logLik methods for objects returned from fit_bernoulli()
#' nobs(fit)
#' coef(fit)
#' vcov(fit)
#' logLik(fit)
#' @name bernoulli
NULL
## NULL

#' @rdname bernoulli
#' @export
bernoulli_loglik <- function(prob, data) {
  if (prob < 0 || prob > 1) {
    return(-Inf)
  }
  res <- stats::dbinom(data, 1, prob, log = TRUE)
  class(res) <- "logLikVec"
  return(res)
}

#' @rdname bernoulli
#' @export
fit_bernoulli <- function(data) {
  res <- list()
  res$mle <- mean(data)
  res$nobs <- length(data)
  res$vcov <- as.matrix(res$mle * (1 - res$mle) / res$nobs)
  n1 <- sum(data)
  n0 <- res$nobs - n1
  res$logLik <- n1 * log(res$mle) + n0 * log(1 - res$mle)
  res$data <- data
  class(res) <- "bernoulli"
  return(res)
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
