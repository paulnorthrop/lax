# =================================== evd =================================== #

#' Loglikelihood adjustment of evd fits
#'
#' Description
#'
#' @param x A fitted model object.
#' @param cluster A vector or factor indicating from which cluster the
#'   respective loglikelihood contributions from \code{loglik} originate.
#'   This must have the same length as the vector returned by \code{loglik}.
#'   If \code{cluster} is not supplied then it is set inside
#'   \code{\link[chandwich]{adjust_loglik}} under the assumption that
#'   each observation forms its own cluster.
#'
#'   If the sandwich package
#'   \href{http://dx.doi.org/10.18637/jss.v016.i09}{(Zeleis, 2006)}
#'   is used to estimate the quantities required to adjust the loglikelihood
#'   then \code{cluster} determines whether the variance matrix \code{V}
#'   of the score vector is estimated using \code{\link[sandwich]{meat}}
#'   (\code{cluster} is \code{NULL}) or
#'   \code{\link[sandwich:vcovCL]{meatCL}}
#'   (\code{cluster} is not \code{NULL}).
#' @param use_vcov A logical scalar.  Should we use the \code{vcov} S3 method
#'   for \code{x} (if this exists) to estimate the Hessian of the independence
#'   loglikelihood to be passed as the argument \code{H} to
#'   \code{\link[chandwich]{adjust_loglik}}?
#'   Otherwise, \code{H} is estimated inside
#'   \code{\link[chandwich]{adjust_loglik}} using
#'   \code{\link[stats:optim]{optimHess}}.
#' @param ... Further arguments to be passed to the functions in the
#'   sandwich package \code{\link[sandwich]{meat}}, if \code{cluster = NULL},
#'   or \code{\link[sandwich:vcovCL]{meatCL}}, otherwise.
#' @examples
#' # We need the evd package
#' got_evd <- requireNamespace("evd", quietly = TRUE)
#'
#' if (got_evd) {
#'   library(evd)
#'   # An example from the evd::fgev documentation
#'   uvdata <- evd::rgev(100, loc = 0.13, scale = 1.1, shape = 0.2)
#'   M1 <- evd::fgev(uvdata, nsloc = (-49:50)/100)
#'   adj_fgev <- alogLik(M1)
#'   summary(adj_fgev)
#'
#'   # An example from Chandler and Bate (2007)
#'   y <- c(chandwich::owtemps[, "Oxford"], chandwich::owtemps[, "Worthing"])
#'   x <- rep(c(-1, 1), each = length(y) / 2)
#'   owfit <- evd::fgev(y, nsloc = x)
#'   year <- rep(1:(length(y) / 2), 2)
#'   adj_owfit <- alogLik(owfit, cluster = year)
#'   summary(adj_owfit)
#'
#'   # An example from the evd::fpot documentation
#'   uvdata <- rgpd(100, loc = 0, scale = 1.1, shape = 0.2)
#'   M1 <- fpot(uvdata, 1)
#'   adj_fpot <- alogLik(M1)
#'   summary(adj_fpot)
#' }
#' @export
alogLik.evd <- function(x, cluster = NULL, use_vcov = TRUE, ...) {
  # List of evd objects supported
  supported_by_oolax <- list(evd_fgev = c("gev", "uvevd", "evd"),
                             evd_fpot = c("pot", "uvevd", "evd"))
  # Does x have a supported class?
  is_supported <- NULL
  for (i in 1:length(supported_by_oolax)) {
    is_supported[i] <- identical(class(x), unlist(supported_by_oolax[i],
                                                  use.names = FALSE))
  }
  if (!any(is_supported)) {
    stop(paste("x's class", deparse(class(x)), "is not supported"))
  }
  # Set the class
  name_of_class <- names(supported_by_oolax)[which(is_supported)]
  class(x) <- name_of_class
  # Call oola::adjust_object to adjust the loglikelihood
  res <- adj_object(x, cluster = cluster, use_vcov = use_vcov, ...)
  class(res) <- c("oolax", "chandwich")
  return(res)
}
